"""Decoder parameters in Liu et al.'s released injection study."""


def build_decoder(circuit):
    from stimbposd.bp_osd import BPOSD

    return BPOSD(
        circuit.detector_error_model(),
        max_bp_iters=1000,
        bp_method="ms",
        osd_order=10,
        osd_method="osd_cs",
        ms_scaling_factor=0,
    )
