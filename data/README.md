# data

This directory collects processed numerical data used directly in the figures of the paper.

```text
data/
├── qm_pes/           # QM reference PES along the reaction path
├── rdf/              # Radial distribution functions
├── reactive_pmfs/    # PMFs for the reactive EVB simulations (bulk and interface)
└── solvation_pmfs/   # PMFs for solvation and surface affinity (droplets and slabs)
```

## Contents
`qm_pes/`
- `pes_irc_dlpno.dat`\
  DLPNO single-point energies along the intrinsic reaction coordinate (IRC).
  Used as the reference quantum-mechanical PES for EVB parametrization and validation.

`rdf/`
- `rdf_O_Ostar.dat`\
  Radial distribution function $g(r)$ between O and O* atoms (as defined in the main text).

`reactive_pmfs/`
- `pmf_bulk_128w.dat`, `pmf_bulk_500w.dat`, `pmf_bulk_1000w.dat`\
  1D PMFs for the SN2 reaction in bulk water with different system sizes.
- `pmf_slab_z_*.dat`\
  1D PMFs for the reaction at the air–water interface, with the reaction center constrained at different interfacial depths (values of $z$).
- `ecn_z_*` subdirectories (`ecn_z_4`, `ecn_z_7`, `ecn_z_10`, `ecn_z_15`)\
  PMFs resolved by effective coordination number (ECN) for configurations at selected $z$-positions.
  Each directory contains several `pmf_slab_ecn_*.dat` files corresponding to different ECN windows.

`solvation_pmfs/`
PMFs related to solvation and surface affinity:
- `ch3cl_droplet_*w.dat`\
  Solvation/free-energy profiles for CH<sub>3</sub>Cl in droplets of varying sizes (e.g., 500, 1000, 2000 waters).
- `ch3cl_slab_*w.dat`\
  Solvation free energies for CH<sub>3</sub>Cl in slab geometries (varying lateral system sizes).
- `ohr_droplet_*w.dat`, `ohr_slab_*w.dat`
  Analogous data for OH radical species.
- `ohr_coating_droplet_1000w.dat`, `ohr_coating_slab_1000w.dat`
  PMFs describing the behavior of OHR when there is CH<sub>3</sub>Cl surface coating.

## Notes
- All PMFs in this directory are processed and ready to be plotted. The raw WHAM or histogram data and trajectories are located in the corresponding directories in `md_simulations/`.
- Units and zero references (e.g. $k_\mathrm{B}T$, kcal/mol) are documented in the main text and in comments in the processing scripts in `analysis/`.
