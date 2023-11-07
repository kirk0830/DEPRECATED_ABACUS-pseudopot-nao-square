import generation as tg
import analysis as ta
import json

def prepare_test(file: str):

    with open(file, "r") as f:
        work_status = json.load(f)

    work_status = tg.convert_to_absolutepath(work_status)
    test_status = tg.generate_test_status(work_status)

    tg.generate_test(test_status)

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
        pass

if __name__ == "__main__":

    main("generation", "input.json")

