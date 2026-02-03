# md_simulations/rdf

Example setup for computing the radial distribution function $g(r)$ between O and O* atoms.

## Contents

- `sys.psf`, `sys.crd`  
  Topology and coordinates of the system.

- `toppar.xml`  
  Force-field parameters for the system.

- `run_rdf.py`, `run_rdf.sh`  
  MD script and submission script used to run a 50 ns simulation suitable for RDF analysis.

- `get_RDF.py`  
  Post-processing script that reads the trajectory and computes the O–O* RDF, writing `RDF_data.txt`.

- `RDF_data.txt`  
  Example RDF data file produced by `get_RDF.py`.

## Workflow

1. Run the trajectory (if needed):

   ```bash
   ./run_rdf.sh
   ```
2. Compute the RDF:
   ```bash
   python get_RDF.py
   ```
3. The resulting `RDF_data.txt` (which is the processed data in `data/rdf/rdf_O_Ostar.dat`) can be compared with the corresponding figure in the paper.
