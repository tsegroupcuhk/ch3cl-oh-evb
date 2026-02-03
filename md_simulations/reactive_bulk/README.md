# md_simulations/reactive_bulk

This directory contains EVB-based simulations of the SN2 reaction in **bulk water** for different system sizes:

```text
reactive_bulk/
├── 128w/
├── 500w/
└── 1000w/
```

Each subdirectory corresponds to a different number of water molecules and has an identical internal layout.

## Common directory structure (e.g. `128w/`)
```text
128w/
├── config_bank/        # Pre-generated restart configurations for SMD windows
├── CV/                 # Collective variable (CV) file(s)
├── eqm/                # Equilibration input scripts and EVB module
├── pull-1/             # Example SMD run and analysis
├── template/           # Template SMD job folder used to generate pull directories
└── generate_running_folder_1-150.sh  # Helper script to create many SMD windows
```

### Key components
- `config_bank/`\
Restart files `*.rst` used as starting points for multiple SMD realizations (up to 200 configurations per system).

- `CV/1.cv`\
Example collective variable (CV) time series files obtained from SMD simulations.

- `eqm/`\
Equilibration scripts:
  - `eqm.sh`, `eqm.py`: run _NVT_ equilibration with the EVB model.
  - `evb/`: local copy of the EVB code.
  - `inits/`: topology, coordinates, parameter files (`*.psf`, `*.crd`, `para.json`, `toppar.*`, etc.).
  - `restart/`: example restart file(s) after equilibration (`restart.rst`).

- `pull-1/`\
  Example SMD window:
  - `runSMD.py`, `SMD-SN2.sh`: scripts to run a single SMD trajectory.
  - `SMD.cv`: resulting CV time series (same as `CV/1.cv`).
  - `nvtrun.dcd`, `nvtrun.log`: example trajectory and log file.
  - `plot_cv_PMF.py`: local script to convert CV/work data into an approximate PMF.
  - `evb/`, `inits/`: as in `eqm/`.

- `template/`\
Template directory used together with `generate_running_folder_1-150.sh` to create multiple SMD windows (`pull-2`,
`pull-3`, …).

### Typical workflow
1. Equilibrate (if not already done):
  ```bash
  cd 128w/eqm
  ./eqm.sh
  ```
2. Generate multiple SMD folders (if desired):
  ```bash
   cd ..
   ./generate_running_folder_1-150.sh
  ```

3. Run SMD simulations in each `pull-*` folder (e.g., via `SMD-SN2.sh`), collect CV and work data.
4. Analyze work distributions using `analysis/smd_jarzynski_analysis.py` or the local `plot_cv_PMF.py`.
5. Combine into a PMF and compare with the processed PMFs in `data/reactive_pmfs/pmf_bulk_*.dat`.

See the main text for details on the SMD protocol (pulling speed, force constant, number of trajectories
used in Jarzynski averaging, etc.).


