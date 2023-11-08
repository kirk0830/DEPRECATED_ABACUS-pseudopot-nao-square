# import two modules
import generation.generation as tg
import analysis.analysis as ta
# for input file parsing
import json

"""
Sub-driver for test generation

Parameters
----------
file : str
    input file, use input_generation.json
"""
def prepare_test(file: str):

    with open(file, "r") as f:
        work_status = json.load(f)

    work_status = tg.convert_to_absolutepath(work_status)
    test_status = tg.generate_test_status(work_status)

    tg.generate_test(test_status)

"""
Sub-driver for test analysis

Parameters
----------
file : str
    input file, use input_analysis.json
"""
def analyze_test(file: str):

    with open(file, "r") as f:
        work_status = json.load(f)
    
    test_status_file = work_status["test_status"]
    with open(test_status_file, "r") as f:
        test_status = json.load(f)

    ta.analyze_test(work_status, test_status)

"""
Pseudopotential test project for ABACUS

Parameters
----------
mode : str
    generation or analysis
file : str
    input file, 
    for generation, use input_generation.json, 
    for analysis, use both input_analysis.json and generated test_status.json
"""
def main(mode: str, file: str):

    if mode == "generation":
        prepare_test(file)
    elif mode == "analysis":
        analyze_test(file)

if __name__ == "__main__":

    main("generation", "input.json")

