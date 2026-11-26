import json
import sys
import tempfile
import unittest
from pathlib import Path
from subprocess import run


class ReceiptReplay(unittest.TestCase):
    def test_tamper_and_resource_limit(self):
        with tempfile.TemporaryDirectory() as d:
            spec, record = Path(d) / "spec.json", Path(d) / "receipt.json"
            leaf = {"op": "family", "sites": ["p"], "minimal": [["p"]]}
            spec.write_text(json.dumps({"kind": "family", "expression": leaf}))
            cmd = [sys.executable, "-m", "ftqc_compose"]
            run([*cmd, "check", spec, "-o", record], check=True, capture_output=True)
            report = json.loads(record.read_text())
            report["result"]["enumerator"] = {"0": 1}
            record.write_text(json.dumps(report))
            result = run(
                [*cmd, "verify", spec, record],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("differs", result.stderr)
            leaf["sites"] += [str(i) for i in range(16)]
            spec.write_text(json.dumps({"kind": "family", "expression": leaf}))
            result = run(
                [*cmd, "check", spec], check=False, capture_output=True, text=True
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("16 physical sites", result.stderr)
