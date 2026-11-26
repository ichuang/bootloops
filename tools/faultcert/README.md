# faultcert

Circuit-level fault properties of Clifford gadgets as receipts: the fault model,
circuit distance, fault distance, strict fault tolerance at low order and a
location census, by exact enumeration, by SAT with a checked LRAT proof, and
cross-checked against Stim.

    python3 -m faultcert analyze surface:3,2 --legs A,B,S   # from tools/
    python3 tools/faultcert/selftest.py                     # battery, from the repository root

See `GUIDE.md`.
