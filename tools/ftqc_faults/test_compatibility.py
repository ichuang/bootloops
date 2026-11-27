"""Compatibility boundaries between the former packages and the two fault engines."""

import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class CompatibilityTests(unittest.TestCase):
    def test_same_module_objects(self):
        for old, new in [
            ("ftqc_audit.circuit", "ftqc_faults.audit.circuit"),
            ("ftqc_audit.certify", "ftqc_faults.audit.certify"),
            ("ftqc_zx.diagram", "ftqc_faults.zx.diagram"),
            ("ftqc_zx.adapters.stim_bridge", "ftqc_faults.zx.adapters.stim_bridge"),
            ("ftqc_audit.injection_cli", "ftqc_extractors.injection_cli"),
            ("ftqc_audit.injection_sampling", "ftqc_extractors.injection_sampling"),
        ]:
            with self.subTest(old=old):
                self.assertIs(
                    importlib.import_module(old), importlib.import_module(new)
                )

    def test_existing_entry_points(self):
        def run(package, *args):
            return subprocess.check_output([sys.executable, "-m", package, *args])

        with tempfile.TemporaryDirectory() as tmp:
            r = str(Path(tmp) / "audit.json")
            run("ftqc_audit", "analyze", "cat:4,verify", "-o", r)
            self.assertEqual(
                json.loads(Path(r).read_text())["schema"], "ftqc_audit-receipt-1"
            )
            run("ftqc_faults", "audit", "verify", r, "--full")
            opts = [
                "contract",
                "steane_5anc",
                "--code",
                "steane_7_1_3",
                "--noise-model",
                "circuit_css",
                "--p",
                "1e-4",
                "--json",
            ]
            legacy = run("ftqc_zx", *opts)
            self.assertEqual(json.loads(legacy)["schema"], "ftqc_zx-contract-v1")
            self.assertEqual(legacy, run("ftqc_faults", "zx", *opts))
