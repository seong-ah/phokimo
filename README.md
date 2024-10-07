# PhoKiMo

## Required conditions

# Directory structure

* Respect to calculation path, folder structure should follow the example below.

directory (calculation_path in a toml)/ ├── 0000_s0min/ │ ├── excited_sp/ │ │ └── tc.out │ └── ground_sp/ │ └── tc.out ├── 0001_s1min/ │ └── ground_sp/ │ └── tc.out └── 0002_s1s0_ci/ └── excited_sp/ └── tc.out

## How to run the main module

Run at the directory of input.toml file

```
python3 -m phokimo input.toml
```

### Inputs 
input.toml: toml file that includes all required information about running mechanism

### Outputs
All output files are stored at the same directory of input.toml file

phokimo.toml: toml file with output data values

phokimo_state_energy.png: scattered graph with energy values of states

phokimo_kinetics_state.png: concentration change of each state vs time

phokimo_kinetics_spin.png: concentration change of each spin vs time

phokimo_expfitting.png: exponential fitting graph on the concentration change of each spin vs time


## How to draw the mechanism graph

Run at the directory of input.toml file

```
python3 -m phokimo.independent_graphing.py
```

## TOML file guideline

Check /phokimo/used_tomls/ for the example tomle files and its guideline

## Some limitations and features to add

