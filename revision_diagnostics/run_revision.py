"""Cacao marker provenance, genotype-correlation and gene-background analyses."""
import argparse
import csv
import hashlib
import importlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback
import uuid
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VERSION = '1.0.0'
STAGES = {
    'harmonization': ('provenance_audit', 'coordinate_audit', 'methods_audit'),
    'genotypes': ('genotype_audit', 'ld_audit'),
    'adjusted_ld': ('genotype_audit', 'within_adjusted_audit', 'cross_peak_audit'),
    'background': ('background_audit',),
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 ** 2), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def within(path, parent):
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


class Context:
    """Read the published table catalog with the original module interface."""
    def __init__(self, stage, supplement, sources, out):
        self.stage, self.supplement, self.raw, self.out = stage, supplement, sources, out
        self.module_root = ROOT / 'modules' / stage
        self.lock = json.loads((self.module_root / 'baseline_lock.json').read_text(encoding='utf-8'))
        self.baseline_config = json.loads((self.module_root / 'baseline_study_config.json').read_text(encoding='utf-8'))
        self.source_paths = {r['role']: sources / r['relative_path'] for r in self.lock['source_files']}
        self.summary, self.checks, self.notes, self.input_audit = {}, [], [], []
        self.book, self.table_map, self.cache = None, {}, {}

    def verify_file(self, path, expected, role):
        if not path.is_file():
            raise FileNotFoundError('Required input missing: ' + str(path))
        actual = digest(path)
        self.input_audit.append(dict(role=role, path=str(path), bytes=path.stat().st_size,
                                     sha256=actual, expected_sha256=expected,
                                     status='PASS' if actual == expected else 'FAIL'))
        if actual != expected:
            raise ValueError('Input differs from the published reference: ' + str(path))

    def validate_inputs(self):
        package = json.loads((ROOT / 'PACKAGE_SHA256.json').read_text(encoding='utf-8'))
        for name, expected in package.items():
            path = (ROOT / name).resolve()
            if not within(path, ROOT):
                raise ValueError('Unsafe package manifest path')
            self.verify_file(path, expected, 'code:' + name)
        self.verify_file(self.module_root / 'baseline_study_config.json',
                         self.lock['analysis_files']['study_config.json'], 'published_configuration')
        self.verify_file(self.supplement, self.lock['supplement_sha256'], 'published_supplement_S1')
        for row in self.lock['source_files']:
            self.verify_file(self.source_paths[row['role']], row['sha256'], 'original_workbook:' + row['role'])

    def read(self, workbook, sheet):
        key = (workbook, sheet)
        if key in self.cache:
            return self.cache[key]
        print('[READ]', workbook + ':' + sheet, flush=True)
        if self.book is None:
            Book = importlib.import_module('xlsx_read').Book
            self.book = Book(self.supplement)
            for _, row in self.book.records('Contents'):
                self.table_map[(row['source_workbook'], row['source_worksheet'])] = (row['worksheet'], row['data_rows'])
        public, count = self.table_map[key]
        rows = [row for _, row in self.book.records(public)]
        self.check('supplement_row_count:' + public, len(rows), count)
        self.cache[key] = rows
        return rows

    def check(self, name, observed, expected):
        passed = observed == expected
        self.checks.append(dict(check=name, observed=observed, expected=expected,
                                status='PASS' if passed else 'FAIL'))
        if not passed:
            print('[CHECK FAIL]', name, 'observed=', observed, 'expected=', expected, flush=True)
        return passed

    def note(self, text):
        self.notes.append(text)

    def write_csv(self, filename, rows, fieldnames=None):
        rows = list(rows)
        if fieldnames is None:
            fieldnames = list(dict.fromkeys(k for row in rows for k in row))
        if not fieldnames:
            fieldnames = ['no_records']
        path = (self.out / 'tables' / filename).resolve()
        if not within(path, (self.out / 'tables').resolve()):
            raise ValueError('Unsafe table output path')
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', newline='', encoding='utf-8-sig') as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction='raise')
            writer.writeheader()
            for row in rows:
                writer.writerow({k: json.dumps(sorted(v) if isinstance(v, set) else v, ensure_ascii=False)
                                 if isinstance(v, (list, dict, set, tuple)) else v for k, v in row.items()})


def run_one(args):
    started = time.monotonic()
    supplement = Path(args.supplement_s1).resolve()
    sources = Path(args.source_root).resolve() if args.source_root else ROOT / 'unused_sources'
    output = Path(args.output_root).resolve()
    for source in (ROOT, supplement, sources):
        if within(output, source) or within(source, output):
            raise ValueError('Output must be separate from code and source inputs: ' + str(source))
    stamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8]
    out = output / ('cacao_revision_' + args.stage + '_' + stamp)
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(ROOT / 'modules' / args.stage))
    ctx = Context(args.stage, supplement, sources, out)
    error = None
    try:
        if shutil.disk_usage(out).free < 2 * 1024 ** 3:
            raise ValueError('At least 2 GiB free output space is required.')
        ctx.validate_inputs()
        for name in STAGES[args.stage]:
            print('[MODULE]', name, flush=True)
            importlib.import_module(name).run(ctx)
        for item in ctx.input_audit:
            ctx.check('input_unchanged:' + item['role'] + ':' + item['path'], digest(Path(item['path'])), item['sha256'])
    except Exception as exc:
        error = type(exc).__name__ + ': ' + str(exc)
        (out / 'ERROR.txt').write_text(traceback.format_exc(), encoding='utf-8')
        print('[ERROR]', error, flush=True)
    finally:
        if ctx.book:
            ctx.book.z.close()
    failed = sum(row['status'] == 'FAIL' for row in ctx.checks)
    status = 'FAILED' if error else ('CHECK_FAILED' if failed else 'COMPLETED')
    report = dict(software='cacao_revision_diagnostics', version=VERSION, stage=args.stage,
                  status=status, mode='PUBLISHED_TABLE_AND_SOURCE_WORKBOOK_REPRODUCTION',
                  python_version=sys.version, platform=platform.platform(),
                  elapsed_seconds=round(time.monotonic() - started, 3),
                  check_count=len(ctx.checks), failed_checks=failed, error=error,
                  scientific_summaries=ctx.summary, interpretation_notes=ctx.notes,
                  scope='Marker provenance, coordinate, genotype-correlation and gene-background analyses from published tables and source workbooks.')
    write_json(out / 'RUN_SUMMARY.json', report)
    write_json(out / 'CHECKS.json', ctx.checks)
    write_json(out / 'INPUT_MANIFEST.json', ctx.input_audit)
    write_json(out / 'OUTPUT_SHA256.json', {p.relative_to(out).as_posix(): digest(p)
                                          for p in sorted(out.rglob('*')) if p.is_file()})
    archive = out.with_suffix('.zip')
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as writer:
        for path in sorted(out.rglob('*')):
            if path.is_file():
                writer.write(path, arcname=out.name + '/' + path.relative_to(out).as_posix())
    with zipfile.ZipFile(archive) as reader:
        if reader.testzip() is not None:
            raise RuntimeError('Result archive failed integrity check')
    print('[STATUS]', status, flush=True)
    print('[RESULT ZIP]', archive, flush=True)
    return 0 if status == 'COMPLETED' else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['all'] + list(STAGES), default='all')
    parser.add_argument('--supplement-s1', required=True, help='Complete published Supplementary Data S1 workbook')
    parser.add_argument('--source-root', help='Root containing the original unedited source workbook paths')
    parser.add_argument('--output-root', required=True)
    args = parser.parse_args()
    if args.stage != 'background' and not args.source_root:
        parser.error('--source-root is required except for the background stage')
    if args.stage == 'all':
        for stage in STAGES:
            command = [sys.executable, '-B', str(Path(__file__).resolve()), '--stage', stage,
                       '--supplement-s1', args.supplement_s1, '--source-root', args.source_root,
                       '--output-root', args.output_root]
            result = subprocess.run(command, check=False)
            if result.returncode:
                return result.returncode
        return 0
    return run_one(args)


if __name__ == '__main__':
    sys.exit(main())
