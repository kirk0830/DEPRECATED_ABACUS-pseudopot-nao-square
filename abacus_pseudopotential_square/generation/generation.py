import os
import json
import time

import abacus_pseudopotential_square.base.scans as sc
import abacus_pseudopotential_square.base.file_io as fio
import itertools as it
import abacus_pseudopotential_square.generation.gppt as gppt
import abacus_pseudopotential_square.generation.glpt as glpt

PROGRAM_PATH = "./abacus_pseudopotential_square/"

def convert_to_absolutepath(work_status: dict) -> dict:

    """ convert all paths to absolute paths """
    work_status["work_folder"] = os.path.abspath(work_status["work_folder"])
    fio.check_folder_availability(work_status["work_folder"])
    work_status["pseudopotentials"]["pseudopot_paths"]["work"] = os.path.abspath(PROGRAM_PATH + work_status["pseudopotentials"]["pseudopot_paths"]["work"])
    fio.check_folder_availability(work_status["pseudopotentials"]["pseudopot_paths"]["work"])
    work_status["pseudopotentials"]["pseudopot_paths"]["resources"] = os.path.abspath(PROGRAM_PATH + work_status["pseudopotentials"]["pseudopot_paths"]["resources"])
    fio.check_folder_availability(work_status["pseudopotentials"]["pseudopot_paths"]["resources"])
    work_status["numerical_orbitals"]["nao_paths"]["work"] = os.path.abspath(PROGRAM_PATH + work_status["numerical_orbitals"]["nao_paths"]["work"])
    fio.check_folder_availability(work_status["numerical_orbitals"]["nao_paths"]["work"])
    work_status["numerical_orbitals"]["nao_paths"]["resources"] = os.path.abspath(PROGRAM_PATH + work_status["numerical_orbitals"]["nao_paths"]["resources"])
    fio.check_folder_availability(work_status["numerical_orbitals"]["nao_paths"]["resources"])

    return work_status

def import_info_to_test_status(
        test_status: dict, 
        valid_data: dict, 
        system: str, 
        test_name: str, 
        element: str, 
        pseudo_element: str, 
        nao_type_element = "") -> dict:

    if nao_type_element == "":
        test_status["systems"][system][test_name]["pseudopotentials"]["files"][element] = valid_data[element][pseudo_element]["file"]
        test_status["systems"][system][test_name]["pseudopotentials"]["info"][element] = {
            "kind": valid_data[element][pseudo_element]["kind"],
            "version": valid_data[element][pseudo_element]["version"],
            "appendix": valid_data[element][pseudo_element]["appendix"]
        }
    else:
        test_status["systems"][system][test_name]["numerical_orbitals"]["files"][element] = valid_data[element][pseudo_element][nao_type_element]["file"]
        test_status["systems"][system][test_name]["numerical_orbitals"]["info"][element] = {
            "type": valid_data[element][pseudo_element][nao_type_element]["type"],
            "rcut": valid_data[element][pseudo_element][nao_type_element]["rcut"],
            "appendix": valid_data[element][pseudo_element][nao_type_element]["appendix"]
        }
    return test_status

def generate_test_status(work_status: dict) -> dict:
    
    """ generate test_status for passing test information """
    test_status = {
        "global": {
            "test_mode": work_status["global"]["test_mode"],
            "analysis_items": work_status["global"]["analysis_items"],
            "software": work_status["global"]["software"]
        },
        "calculation": {
            "basis_type": work_status["calculation"]["basis_type"],
            "functionals": work_status["calculation"]["functionals"],
            "ecutwfc": work_status["calculation"]["ecutwfc"],
            "stress_deformation_ratios": work_status["calculation"]["stress_deformation_ratios"]
        },
        "systems": {},
        "paths": {
            "root": os.getcwd(),
            "work_folder": work_status["work_folder"],
            "pseudopotential": {
                "work": work_status["pseudopotentials"]["pseudopot_paths"]["work"],
                "resources": work_status["pseudopotentials"]["pseudopot_paths"]["resources"]
            },
            "numerical_orbital": {
                "work": work_status["numerical_orbitals"]["nao_paths"]["work"],
                "resources": work_status["numerical_orbitals"]["nao_paths"]["resources"]
            }
        }
    }
    # get valid pseudopotentials
    valid_pseudopotentials = sc.scan_valid_pseudopotentials(work_status)
    if work_status["calculation"]["basis_type"] == "lcao":
        valid_numerical_orbitals = sc.scan_valid_numerical_orbitals(work_status, valid_pseudopotentials)
        if len(list(valid_numerical_orbitals.keys())) != len(list(valid_pseudopotentials.keys())):
            print("Error: number of elements in valid_numerical_orbitals is not equal to number of elements in valid_pseudopotentials.")
            exit(1)
    # loop over systems, for every system, find all possible combinations of pseudopotentials
    # and corresponding numerical orbitals
    # note that it may be different for different systems

    
    for system in work_status["systems"].keys():
        # goal is to store all possible tests to test_status system by system
        test_status["systems"][system] = {}
        # allocate memory for storing valid pseudopotentials list element by element
        valid_pseudopots_to_combine = []
        # get elements in present system
        elements = sc.scan_elements(system)

        for element in elements:
            if work_status["calculation"]["basis_type"] == "pw":
                # in this case... because the valid_pseudopotentials has format:
                # {
                #     "As": {
                #         "pd_04_spd": {
                #             "file": ...,
                #             "kind": ...,
                #             "version": ...,
                #             "appendix": ...,
                #         },
                #         "sg15_10": {
                #             "file": ...,
                #             "kind": ...,
                #             "version": ...,
                #             "appendix": ...,
                #         }
                #         ...
                #     },
                # }
                # so we need to get the keys of the dictionary
                valid_pseudopots_to_combine.append([[valid_pseudopot] for valid_pseudopot in list(valid_pseudopotentials[element].keys())])
                # thus for every element, valid_pseudopots_to_combine will be like:
                # [
                #     [["pd_04_spd"], ["sg15_10"], ...],
                #     [["pd_04_spd"], ["sg15_10"], ...],
                #     ...
                # ]
            elif work_status["calculation"]["basis_type"] == "lcao":
                # in this case... valid_numerical_orbitals has format:
                # {
                #     "As": {
                #         "pd_04_spd": {
                #             "D6": {
                #                "type": "DZP",
                #                "rcut": 6,
                #                "appendix": "",
                #                "file": "As_gga_6au_100Ry_2s2p1d.orb",
                #             },
                #             ...
                #             "T7": {
                #                 "type": "TZDP",
                #                 "rcut": 7,
                #                 "appendix": "",
                #                 "file": "As_gga_7au_100Ry_3s3p2d.orb",
                #             },
                #             ...
                #         }
                #         "sg15_10": {
                #             "D6": {
                #                 "type": "DZP",
                #                 "rcut": 6,
                #                 "appendix": "",
                #                 "file": "As_gga_6au_100Ry_2s2p1d.orb",
                #             },
                #             ...
                #             "D6_osci_rm": {
                #                 "type": "DZP",
                #                 "rcut": 6,
                #                 "appendix": "osci_rm",
                #                 "file": "As_gga_6au_100Ry_2s2p1d_osci_rm.orb",
                #             },
                #             ...
                #         }
                #     },
                #     ...
                # }
                # so we need to get the keys of the dictionary
                valid_pseudopot_nao_type_element = []
                for pseudopotential in valid_numerical_orbitals[element].keys():
                    for nao_type in valid_numerical_orbitals[element][pseudopotential].keys():
                        valid_pseudopot_nao_type_element.append([pseudopotential, nao_type])
                valid_pseudopots_to_combine.append(valid_pseudopot_nao_type_element)
                # thus for every element, valid_pseudopots_to_combine will be like:
                # [
                #     pseudopot and corresponding numerical orbital in each []
                #     [["pd_04_spd", "D6"], ["sg15_10", "D6"], ...], # element 1
                #     [["pd_04_spd", "D6"], ["sg15_10", "D6"], ...], # element 2
                #     ...
                # ]
        # make Cartesian product of all valid pseudopotentials
        # test combinations (per) system
        test_combinations_system = list(it.product(*valid_pseudopots_to_combine))
        # such that each combination will be a list of tuples of list:
        # for lcao, it will be like:
        # [
        #             element1             element2           element3
        #     test1: (["pd_04_spd", "D6"], ["sg15_10", "D6"], ...),
        #     test2: (["pd_04_spd", "D6"], ["sg15_10", "D6_osci_rm"]),
        #     ...
        # ]
        # for pw, it will be like:
        # [
        #             element1       element2     element3
        #     test1: (["pd_04_spd"], ["sg15_10"], ...),
        #     ...
        # ]
        # loop over all combinations
        for test in test_combinations_system:
            # convert test name from tuple (["pd_04_spd", "D6"], ["sg15_10", "D6"], ...) to pd04spdsg1510..._D6D6...., each test is equivalent with one tuple
            test_name = ""
            for item in range(len(test[0])): # item can only be 1 or 2, 1 is for pseudopotential without ecutwfc tests, 
                                             #                          2 is for numerical orbitals or pseudopotential with ecutwfc tests
                for axis_element in test: # axis can be many, for instance InP, it is 2, BaTiO3, it is 3
                    test_name += axis_element[item].replace("_", "")
                test_name += "_" # add _ between pseudopotential and numerical orbitals
            if work_status["calculation"]["basis_type"] == "pw":
                # add ecutwfc to test name if there are more than one ecutwfc
                for ecutwfc in test_status["calculation"]["ecutwfc"]:
                    if len(test_status["calculation"]["ecutwfc"]) > 1:
                        test_name += "ecut" + str(ecutwfc)
                    else:
                        test_name = test_name[:-1] # remove the last _
                    print("prepare test: " + test_name + " for system: " + system)
                    # initialize test_status for this test
                    test_status["systems"][system][test_name] = {
                        "elements": elements,
                        "pseudopotentials": {
                            "files": {},
                            "info": {}
                        },
                        "numerical_orbitals": {
                            "files": {},
                            "info": {}
                        },
                        "ecutwfc": ecutwfc
                    }
                    for ie, element in enumerate(elements): # so we can get the index of element
                        pseudo_element = test[ie][0]
                        # get pseudopotential information and copy to test_status
                        test_status = import_info_to_test_status(test_status, valid_pseudopotentials, system, test_name, element, pseudo_element)
                    test_name = test_name.split("ecut")[0] # remove ecutwfc from test name

            elif work_status["calculation"]["basis_type"] == "lcao":
                test_name = test_name[:-1] # remove the last _
                print("prepare test: " + test_name + " for system: " + system)
                # initialize test_status for this test
                test_status["systems"][system][test_name] = {
                    "elements": elements,
                    "pseudopotentials": {
                        "files": {},
                        "info": {}
                    },
                    "numerical_orbitals": {
                        "files": {},
                        "info": {}
                    }
                }
                # loop over elements
                for ie, element in enumerate(elements): # so we can get the index of element
                    pseudo_element = test[ie][0]
                    nao_type_element = test[ie][1]
                    # get pseudopotential information and copy to test_status
                    test_status = import_info_to_test_status(test_status, valid_pseudopotentials, system, test_name, element, pseudo_element)
                    # get numerical orbital information and copy to test_status
                    test_status = import_info_to_test_status(test_status, valid_numerical_orbitals, system, test_name, element, pseudo_element, nao_type_element)
    
    if work_status["global"]["save_test_status"] or "save_test_status" not in work_status["global"].keys():
        # save test_status to json file, name marked with time stamp
        # note that the test_status is saved in the work folder
        print("save test_status to json file")
        os.chdir(test_status["paths"]["work_folder"])
        test_status_json = "test_status_" + time.strftime("%Y%m%d_%H%M%S", time.localtime()) + ".json"
        with open(test_status_json, "w") as f:
            json.dump(test_status, f, indent=4)
        print("test_status saved to " + test_status_json)
    else:
        print("warning: test_status is necessary for analyzing test results, but it is not saved according to input.")
    # change back to root folder
    os.chdir(test_status["paths"]["root"])
    return test_status

def generate_test(test_status: dict) -> None:

    if test_status["calculation"]["basis_type"] == "pw":
        gppt.generate_pw_pseudopotential_test(test_status)
    elif test_status["calculation"]["basis_type"] == "lcao":
        glpt.generate_lcao_pseudopotential_test(test_status)
    else:
        print("Error: basis_type not recognized.")
        exit(1)

