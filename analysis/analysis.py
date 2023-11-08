import analysis.abg as abg
import analysis.astr as astr
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
    if len(analyze_items) != len(test_status["global"]["analysis_items"]):
        print("Inconsistent analysis items defined in input file.")
        exit(1)
        
    for analyze_item in analyze_items:
        if analyze_item == "band_gap":
            analyze_band_gap(work_status, test_status)
        elif analyze_item == "stress":
            analyze_stress(work_status, test_status)

def analyze_band_gap(work_status: dict, test_status: dict):

    if test_status["global"]["software"] == "Quantum ESPRESSO":
        raise NotImplementedError("Band gap analysis for Quantum ESPRESSO is not implemented yet.")
    elif test_status["global"]["software"] == "ABACUS":
        raise NotImplementedError("Band gap analysis for ABACUS is not implemented here. Instead, use ABACUS Test in-build function.")
    else:
        print("Error: software not recognized.")
        exit(1)

def analyze_stress(work_status: dict, test_status: dict):

    if test_status["global"]["software"] == "Quantum ESPRESSO":
        raise NotImplementedError("Stress analysis for Quantum ESPRESSO is not implemented yet.")
    elif test_status["global"]["software"] == "ABACUS":
        raise NotImplementedError("Stress analysis for ABACUS is not implemented yet.")
    else:
        print("Error: software not recognized.")
        exit(1)