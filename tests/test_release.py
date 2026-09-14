"""Version, configuration and package boundary checks."""
from pathlib import Path
import ast,json,re,unittest
from cacao_marker_harmonization import __version__
from cacao_marker_harmonization.configuration import load_study

ROOT=Path(__file__).resolve().parents[1]
class ReleaseTests(unittest.TestCase):
    def test_versions_agree(self):
        project=(ROOT/'pyproject.toml').read_text(encoding='utf-8');citation=(ROOT/'CITATION.cff').read_text(encoding='utf-8')
        self.assertIn('version = "'+__version__+'"',project)
        self.assertIn('version: '+__version__,citation)
    def test_configuration_copy(self):
        self.assertEqual(json.loads((ROOT/'config/study.json').read_text(encoding='utf-8')),json.loads((ROOT/'src/cacao_marker_harmonization/study.json').read_text(encoding='utf-8')))
    def test_python310_syntax(self):
        for p in (ROOT/'src').rglob('*.py'):
            ast.parse(p.read_text(encoding='utf-8'),feature_version=(3,10))
    def test_no_figure_dependencies(self):
        text=(ROOT/'pyproject.toml').read_text(encoding='utf-8').lower()
        for word in ['matplotlib','seaborn','plotly']:self.assertNotIn(word,text)
    def test_no_history_directories_required(self):
        cfg=load_study(ROOT/'config/study.json');keys=json.dumps(cfg)
        for word in ['required_mapping_version','required_impact_version','required_refinement_version','baseline_run','evidence_run','run_paths']:
            self.assertNotIn(word,keys)
