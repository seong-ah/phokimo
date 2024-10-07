# PhoKiMo

PhoKiMo is a general kinetic model based on rate law. By utilizing given mechanism on toml file, you can calculate the kinetic scheme of the reaction. This model is using calculation result of `terachem`.

## Prerequisites

Directory from the calculation path should look like the structure below. Numbering does not matter, but folder name should **ENDS** with the `folder_name` at the toml file.

```
calculation_path/
├── 0000_aaa/
│   └── sp/
│       └── tc.out
├── 0001_bbb/
│   ├── sp/
│   │   └── tc.out
│   └── ground_sp/
│       └── tc.out
└── 0002_ccc/
    └── sp/
        └── tc.out
```

## How to run the kinetic model

Kinetic model should be run on the cluster that has terachem output file.

```
python3 -m phokimo input.toml
```

### Inputs

`input.toml`: toml file that includes information of the target mechanism

### Outputs

`phokimo.toml`: 

## How to dramw the mechansim diagram

Drawing mechanism on the cluster is not supported, so it should be run independently on your local environment with the toml file.

```
python3 -m phokimo.mechanism input.toml
```

## TOML file guideline

Check /phokimo/used_tomls/ for the example tomle files and its guideline.

## Some limitations and features to add

