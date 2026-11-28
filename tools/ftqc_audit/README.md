# ftqc_audit

The location-fault engine of [ftqc_faults](../ftqc_faults/GUIDE.md). Its code
is in `tools/ftqc_faults/audit/`; this directory keeps the old import name and
command line working.

For one scheduled Clifford circuit it gives the fault model (detectors fired,
observables flipped and output damage for every Pauli at every location),
circuit and fault distance as intervals with replayable witnesses, strict
fault tolerance at low order, and a location census, by exact enumeration and
by SAT with a re-checked LRAT proof, optionally cross-checked against Stim.

    export PYTHONPATH="$PWD/tools"                       # from the BootLoops root
    python3 -m ftqc_faults audit analyze surface:3,2 --legs A,B,S
    python3 tools/ftqc_audit/selftest.py                 # battery

See [`GUIDE.md`](GUIDE.md).
