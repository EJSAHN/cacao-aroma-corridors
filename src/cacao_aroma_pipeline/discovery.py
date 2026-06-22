from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd

from cacao_aroma_pipeline.constants import (
    SUPPORTED_ARCHIVE_EXTENSIONS,
    SUPPORTED_DOCUMENT_EXTENSIONS,
    SUPPORTED_SPREADSHEET_EXTENSIONS,
)
from cacao_aroma_pipeline.models import RawFileRecord
from cacao_aroma_pipeline.utils import ensure_dir, file_md5


def scan_raw_files(raw_dir: Path, allowed_extensions: set[str] | None = None) -> list[RawFileRecord]:
    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw directory does not exist: {raw_dir}")
    allowed = allowed_extensions
    records: list[RawFileRecord] = []
    for path in sorted(raw_dir.rglob("*")):
        if not path.is_file():
            continue
        ext = path.suffix.lower()
        if allowed is not None and ext not in allowed:
            continue
        records.append(
            RawFileRecord(
                path=path,
                relative_path=str(path.relative_to(raw_dir)),
                extension=ext,
                size_bytes=path.stat().st_size,
                md5=file_md5(path),
            )
        )
    return records


def build_inventory_frame(records: list[RawFileRecord]) -> pd.DataFrame:
    rows = []
    hash_counts: dict[str, int] = {}
    for record in records:
        hash_counts[record.md5] = hash_counts.get(record.md5, 0) + 1
    for record in records:
        category = "other"
        if record.extension in SUPPORTED_SPREADSHEET_EXTENSIONS:
            category = "spreadsheet"
        elif record.extension in SUPPORTED_ARCHIVE_EXTENSIONS:
            category = "archive"
        elif record.extension in SUPPORTED_DOCUMENT_EXTENSIONS:
            category = "document"
        rows.append(
            {
                "relative_path": record.relative_path,
                "extension": record.extension,
                "size_bytes": record.size_bytes,
                "md5": record.md5,
                "duplicate_hash_count": hash_counts[record.md5],
                "category": category,
                "source_archive": record.source_archive,
                "extracted_member": record.extracted_member,
            }
        )
    return pd.DataFrame(rows)


def extract_supported_archive_members(
    archive_record: RawFileRecord,
    *,
    raw_dir: Path,
    staging_dir: Path,
    allowed_extensions: set[str],
) -> list[RawFileRecord]:
    ensure_dir(staging_dir)
    extracted: list[RawFileRecord] = []
    with zipfile.ZipFile(archive_record.path) as archive:
        for member in archive.infolist():
            member_path = Path(member.filename)
            if member.is_dir():
                continue
            if member_path.suffix.lower() not in allowed_extensions:
                continue
            destination = staging_dir / archive_record.md5 / member.filename
            ensure_dir(destination.parent)
            with archive.open(member, "r") as source, destination.open("wb") as sink:
                sink.write(source.read())
            extracted.append(
                RawFileRecord(
                    path=destination,
                    relative_path=str(destination.relative_to(staging_dir.parent)),
                    extension=destination.suffix.lower(),
                    size_bytes=destination.stat().st_size,
                    md5=file_md5(destination),
                    source_archive=str(archive_record.path.relative_to(raw_dir)),
                    extracted_member=member.filename,
                )
            )
    return extracted


def collect_candidate_spreadsheets(
    records: list[RawFileRecord],
    *,
    raw_dir: Path,
    staging_dir: Path,
    include_archives: bool,
    archive_member_extensions: set[str],
) -> list[RawFileRecord]:
    candidates: list[RawFileRecord] = []
    seen_hashes: set[str] = set()

    for record in records:
        if record.extension in SUPPORTED_SPREADSHEET_EXTENSIONS:
            if record.md5 not in seen_hashes:
                seen_hashes.add(record.md5)
                candidates.append(record)

    if include_archives:
        for record in records:
            if record.extension not in SUPPORTED_ARCHIVE_EXTENSIONS:
                continue
            for extracted in extract_supported_archive_members(
                record,
                raw_dir=raw_dir,
                staging_dir=staging_dir,
                allowed_extensions=archive_member_extensions,
            ):
                if extracted.md5 in seen_hashes:
                    continue
                seen_hashes.add(extracted.md5)
                candidates.append(extracted)

    return sorted(candidates, key=lambda x: (x.source_archive or "", x.relative_path))
