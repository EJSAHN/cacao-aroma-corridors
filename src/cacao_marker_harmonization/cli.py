"""Command-line access to the raw-input workflow and numerical validation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import unittest
from . import __version__
from .configuration import load_study
from .errors import DataError


def main(argv=None):
    p=argparse.ArgumentParser(description='Cacao marker harmonization: raw-input, values-only analysis; no figures.')
    p.add_argument('command',choices=['selftest','preflight','run','verify-results'])
    p.add_argument('--config',type=Path,default=(Path(__file__).resolve().parents[2]/'config/study.json' if (Path(__file__).resolve().parents[2]/'config/study.json').is_file() else Path(__file__).with_name('study.json')))
    for key in ['source-root','output-root','scratch','cache-root','tools-root','reference-cache','reference-dir','reference-fasta','reference-gff','bowtie2-dir','analysis-run','reference-coverage','reference-ld']:
        p.add_argument('--'+key,type=Path)
    p.add_argument('--offline',action='store_true',help='No automatic reference/tool download.')
    p.add_argument('--threads',type=int)
    args=p.parse_args(argv)
    try:
        if sys.version_info<(3,10):raise DataError('Python 3.10 or newer is required.')
        if args.command=='selftest':
            tests=Path(__file__).resolve().parents[2]/'tests'
            if not tests.is_dir():raise DataError('Use the source distribution for the bundled tests.')
            result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover(str(tests)))
            return 0 if result.wasSuccessful() else 1
        if args.command=='verify-results':
            if args.analysis_run is None:raise DataError('--analysis-run is required.')
            from .validation import verify_results
            from .workflow import finalize
            result=verify_results(args.analysis_run,args.reference_coverage,args.reference_ld)
            finalize(args.analysis_run.resolve())
            print('[AGREEMENT]',json.dumps(result,ensure_ascii=False),flush=True)
            return 0 if result['status']=='PASS' else 1
        if args.source_root is None:raise DataError('--source-root is required.')
        args.config=args.config.resolve();cfg=load_study(args.config)
        if args.threads is not None and not 1<=args.threads<=8:raise DataError('--threads must be in 1..8.')
        args.source_root=args.source_root.expanduser().resolve()
        from .workflow import source_paths, run
        raw,paths=source_paths(args.source_root,cfg)
        if args.command=='preflight':
            print('Cacao marker harmonization',__version__)
            print('Read-only source root:',raw)
            for role,path in paths.items():print('[INPUT]',role,':',path.relative_to(raw))
            print('Configuration and source paths resolved. Scientific content is checked during run.')
            return 0
        if args.output_root is None:raise DataError('--output-root is required; no implicit writes to the current directory.')
        args.output_root=args.output_root.expanduser().resolve();base=args.output_root.parent
        for key,default in [('scratch',base/'scratch'),('cache_root',base/'alignment_cache'),('tools_root',base/'tools'),('reference_cache',base/'reference_cache')]:
            setattr(args,key,(getattr(args,key) or default).expanduser().resolve())
        for key in ['reference_dir','reference_fasta','reference_gff','bowtie2_dir']:
            if getattr(args,key) is not None:setattr(args,key,getattr(args,key).expanduser().resolve())
        out=run(args,cfg)
        return 0
    except KeyboardInterrupt:
        print('[STOPPED] Interrupted. Partial outputs are not completed analyses.',file=sys.stderr);return 130
    except Exception as exc:
        print('[ERROR]',type(exc).__name__+':',str(exc),file=sys.stderr);return 1

if __name__=='__main__':raise SystemExit(main())
