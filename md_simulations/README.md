# md_simulations

This directory contains representative molecular dynamics (MD) setups used in the paper. 
They are organized by the physical quantity being computed (RDF, reactive PMFs, solvation PMFs) and by system type (bulk, slab, droplet).

```text
md_simulations/
├── rdf/              # RDF calculation for O–O* distances
├── reactive_bulk/    # EVB reaction in bulk water (128w, 500w, 1000w)
├── reactive_slab/    # EVB reaction at the air–water interface (various z, ECN)
└── solvation_droplet/# Solvation and surface affinity of CH3Cl and OHR in droplets
```

Each subdirectory contains:

- **Initialization files** (`inits/`): structures (*.psf, *.crd, *.pdb), EVB/FF parameters (para.json, toppar.*, off.xml), and small helper scripts.
- **Equilibration scripts** (`eqm/`): input and job submission scripts to equilibrate the system.
- **Production runs** (`pull-1/`, `runfolder`, `US.sh`, etc.): steered MD or umbrella sampling setups.
- **Analysis scripts** (`plot_cv_PMF.py`, `get_RDF.py`, etc.): to compute PMFs or RDFs from trajectories.
- **Representative trajectories and CV files**: selected `*.dcd`, `*.cv`, and restart files illustrating typical output.

Many subfolders (e.g. specific $z$-positions or droplet sizes) include their own `README` files describing that particular system and workflow.

For a high-level overview of what each subdirectory corresponds to in the paper, see:

- `md_simulations/reactive_bulk/README`
- `md_simulations/reactive_slab/README`
- `md_simulations/solvation_droplet/README`
- `md_simulations/rdf/README`
