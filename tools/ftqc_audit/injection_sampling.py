"""Compatibility import for the extractor injection study."""
if __name__ == "__main__":
    from ftqc_extractors.injection_sampling import main
    main()
else:
    import sys
    from ftqc_extractors import injection_sampling
    sys.modules[__name__] = injection_sampling
