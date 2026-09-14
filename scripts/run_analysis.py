#!/usr/bin/env python3
"""Run from a source checkout without an editable installation."""
import sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from cacao_marker_harmonization.cli import main
if __name__=='__main__':raise SystemExit(main())
