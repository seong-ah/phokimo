import toml
import os

from phokimo.src.io.terachem import TeraChemOutputReader
from phokimo.src.toml_reader import TomlReader
from phokimo.src.rate_constants import RateCalculator, AdhocRelaxation, RelaxationTheory

# Get the directory of the currently executing Python script
current_dir = os.path.dirname(os.path.abspath(__file__))

# Assume the TOML file is in the same directory as the Python script
toml_file_path = os.path.join(current_dir, "sample.toml")

# Read the TOML file
with open(toml_file_path, "r") as f:
    toml_data = toml.load(f)

toml_preparation = TomlReader(toml_file_path)

value = toml_preparation._normal_modes()
#print(value)

state = toml_preparation._ts_existence(1, 4)
#print(state)

constant = RateCalculator()
rate = constant.relaxation_theory.compute_rate(1, 63, 20)
print(rate)