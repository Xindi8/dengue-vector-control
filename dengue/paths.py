"""File locations, relative to the project folder, so code runs from anywhere."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AB = ROOT / "data_and_literature"
AB_DATA = AB / "data"
PARAMETERS_JSON = AB / "parameters" / "parameters.json"
PARAMETERS_CSV = AB / "parameters" / "parameters.csv"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
