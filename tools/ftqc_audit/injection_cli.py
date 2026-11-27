"""Compatibility import for the extractor injection study."""
import sys
from ftqc_extractors import injection_cli
sys.modules[__name__] = injection_cli
