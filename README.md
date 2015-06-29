# Coarsegraining-Flow

Code to take DFAX series from atomistic to coarsegrained (and possibly run some analysis and data collection)

Flowchart:

(1) Create a topology file by using DFAG & martinize.py and merging the two results

(2) Transform topology file into VOTCA-usable form

(3) Use VOTCA to get data on bonds, angles, and dihedrals from atomistic system

(4) Create tabulated potentials and run a short (~30 ns) CG run

(5) Do analysis of resultant CG run

(6) Assuming it's fine, move to doing long production runs
