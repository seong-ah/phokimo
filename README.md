# PhoKiMo

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

