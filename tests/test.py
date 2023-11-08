from ..base import file_io as fio
import json

if __name__ == "__main__":
    
    # load test_status_test.json
    with open("test_status_test.json", 'r') as f:
        test_status = json.load(f)
    pseudpotentials_to_test = test_status["systems"]["InAs"]["sg1510sg1510_D7D10"]["pseudopotentials"]
    numerical_orbitals_to_test = test_status["systems"]["InAs"]["sg1510sg1510_D7D10"]["numerical_orbitals"]
    cif_file_to_test = "./tests/support/InAs_cub.cif"
    fio.cif_to_STRU(cif_file_to_test, pseudpotentials_to_test, numerical_orbitals_to_test)