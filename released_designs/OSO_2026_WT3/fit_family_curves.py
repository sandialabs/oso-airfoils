# load the jsons of polynomial data
import json
from pathlib import Path
import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()
path_to_oso  = path_to_here.parent.parent   # released_designs/pareto_data → oso-airfoils/

json_path = path_to_oso / "released_designs/pareto_data/unconstrained/tau_21/shape_functions.json"
with open(json_path, "r") as f:
    shape_functions = json.load(f)
print(shape_functions.keys())