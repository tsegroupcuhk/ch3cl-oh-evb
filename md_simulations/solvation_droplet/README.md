# md_simulations/solvation_droplet

This directory contains umbrella sampling (US) and steered MD setups for solvation and surface affinity of CH<sub>3</sub>Cl and OH radical in water droplets of various sizes.

The naming convention is:

```text
solvation_droplet/
├── ch3cl_1000w/         # Full US setup for CH3Cl in a 1000-water droplet
├── _ch3cl_500w/         # Documentation folder for 500-water droplet runs (README only)
├── _ch3cl_2000w/        # Documentation folder for 2000-water droplet runs (README only)
├── _ohr_500w/           # Documentation folder for OHR in 500-water droplet
├── _ohr_1000w/          # Documentation folder for OHR in 1000-water droplet
├── _ohr_2000w/          # Documentation folder for OHR in 2000-water droplet
└── _ohr_coating_1000w/  # Documentation folder for OHR with CH3Cl coating simulations
```

Folders starting with `_` typically contain a `README` explaining which data files under `data/solvation_pmfs/` correspond to those simulations, but do not necessarily include full input sets.

## Example: `ch3cl_1000w/`
This folder contains a complete umbrella sampling workflow for CH<sub>3</sub>Cl in a 1000-water droplet.

### Layout
```text
ch3cl_1000w/
├── SMD_to_determine_US_windows/  # Initial SMD to choose US windows and force constants
├── runfolder/                    # Initial US run (generation of restarts files)
├── runfolder-1/                  # Long production US run + analysis
├── submitscript/                 # Batch submission scripts for each US window
├── sys.psf, sys.crd, toppar.xml  # System topology, coordinates, parameters
├── nvt.py, nvt2.py               # Scripts for NVT runs
├── submit.sh, US.sh              # Convenience drivers for running all US windows
└── US_settings.txt               # Example US settings (may be duplicated in subfolders)
```

### SMD to determine US windows
`SMD_to_determine_US_windows/`:
- `eqm.py`, `eqm.sh`: Equilibration of the droplet with CH<sub>3</sub>Cl.
- `pull.py`, `pull.sh`: SMD along the chosen reaction coordinate (i.e., distance from the droplet center).
- `pull.cv`, `pull.dcd`: SMD CV and trajectory.
- `smooth_PMF.py`: Script to smooth the PMF.
- `US_settings_optimizer_kd.py`, `US_settings_func.py`: Scripts/functions to determine US window centers and force constants based on SMD PMF, giving `US_settings.txt`
- `get_rst_for_each_CV_from_cvfile_dcdfile.py`: Scripts/functions to get rst structure files for each US window from SMD .dcd file.
- `nPMF.txt`, `PMF.txt`: SMD-derived PMFs used to design the US windows.
- `US_settings.txt`: Resulting US window settings.

### Umbrella sampling runs
`runfolder/`:
- `CV/0.cv`–`19.cv`: Resulting US CV files for each window.
- `rst/0.rst`–`19.rst`: Starting structures for each window.
- `US_settings.txt`: US parameters (window centers, force constants).
- `prepare_window.py`: Script calling `get_rst_for_each_CV_from_cvfile_dcdfile.py` to generate `metafile` and `rst/0.rst`–`19.rst`.
- `metafile`: WHAM metafile.

`runfolder-1/` (production):
- `CV/`, `rst/`: As above, for the longer production run.
- `example_traj/`: Selected `nvt-*.dcd` trajectories illustrating sampling in each window.
- `metafile`, `metafile-trun`: WHAM metafile.
- `check_overlap.py`: Script to assess histogram overlap between neighboring windows (produces overlapping.png).
- `wham.sh`: Script to run WHAM and produce `WHAM_output`.
- `WHAM_output`: WHAM output files.
- `US_settings.txt`: Copy of the US setup used for this production run.

### Submission scripts
`submitscript/US*.sh`:
- Per-window batch scripts to submit each US window to a cluster scheduler.

Top-level files:
- `submit.sh`: Example high-level script to submit all US windows.
- `US.sh`: Convenience script that may run all windows.

### Relation to processed data
The processed solvation/surface PMFs in `data/solvation_pmfs/ch3cl_droplet_1000w.dat` were obtained by running and analyzing the US setup in `ch3cl_1000w/` using WHAM.

Other droplet sizes and solutes follow analogous protocols; see the `_ch3cl_*w` and `_ohr_*w` READMEs for details.
