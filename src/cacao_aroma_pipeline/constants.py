from __future__ import annotations

DEFAULT_WINDOWS_BP = [0, 10_000, 50_000, 250_000]
SUPPORTED_SPREADSHEET_EXTENSIONS = {".xlsx", ".xls"}
SUPPORTED_ARCHIVE_EXTENSIONS = {".zip"}
SUPPORTED_DOCUMENT_EXTENSIONS = {".pdf"}

STANDARD_MARKER_COLUMNS = [
    "dataset_id",
    "dataset_name",
    "parser_type",
    "source_path",
    "source_archive",
    "marker_id",
    "chromosome",
    "position",
    "alleles",
    "strand",
]

STANDARD_SAMPLE_COLUMNS = [
    "dataset_id",
    "dataset_name",
    "parser_type",
    "source_path",
    "source_archive",
    "sample_id",
    "germplasm_name",
    "dna_sample_id",
    "origin",
    "accession_number",
]

MAX_EXCEL_SHEETNAME_LEN = 31
