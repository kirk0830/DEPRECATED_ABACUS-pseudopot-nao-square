import abg
import astr
def analyze_test(work_status: dict, test_status: dict) -> None:

    analyze_items = []
    for system in work_status["systems"].keys():
        analyze_items_system = work_status["systems"][system]["experimental_values"].keys()
        if len(analyze_items) == 0:
            analyze_items = list(analyze_items_system)
        else:
            if len(analyze_items) != len(analyze_items_system):
                print("Inconsistent experimental values for system " + system)
                exit(1)

    for analyze_item in analyze_items:
        if analyze_item == "band_gap":
            analyze_band_gap(work_status, test_status)
        elif analyze_item == "stress":
            analyze_stress(work_status, test_status)

def analyze_band_gap(work_status: dict, test_status: dict):

    if test_status["general"]["software"] == "Quantum ESPRESSO":
        pass
    elif test_status["general"]["software"] == "ABACUS":
        pass
    else:
        print("Error: software not recognized.")
        exit(1)

def analyze_stress(work_status: dict, test_status: dict):

    if test_status["general"]["software"] == "Quantum ESPRESSO":
        pass
    elif test_status["general"]["software"] == "ABACUS":
        pass
    else:
        print("Error: software not recognized.")
        exit(1)