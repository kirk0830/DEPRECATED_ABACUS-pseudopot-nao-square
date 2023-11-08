import sys
import sys
import os

import base.scans as sc
import itertools as it
import generation.gppt as gppt
import generation.glpt as glpt

def convert_to_absolutepath(work_status: dict) -> dict:

    """ convert all paths to absolute paths """
    work_status["work_folder"] = os.path.abspath(work_status["work_folder"])
    work_status["pseudopotentials"]["pseudopot_paths"]["work"] = os.path.abspath(work_status["pseudopotentials"]["pseudopot_paths"]["work"])
    work_status["pseudopotentials"]["pseudopot_paths"]["resources"] = os.path.abspath(work_status["pseudopotentials"]["pseudopot_paths"]["resources"])
    work_status["numerical_orbitals"]["nao_paths"]["work"] = os.path.abspath(work_status["numerical_orbitals"]["nao_paths"]["work"])
    work_status["numerical_orbitals"]["nao_paths"]["resources"] = os.path.abspath(work_status["numerical_orbitals"]["nao_paths"]["resources"])
    return work_status

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
        # make Cartesian product of all valid pseudopotentials
        # test combinations (per) system
        test_combinations_system = list(it.product(*valid_pseudopots_to_combine))
        # such that each combination will be a list of tuples of list:
        # [
        #     (["pd_04_spd", "D6"], ["sg15_10", "D6"]),
        #     (["pd_04_spd", "D6"], ["sg15_10", "D6_osci_rm"]),
        #     ...
        # ]
        # loop over all combinations
        for test in test_combinations_system:
            # convert test name from (["pd_04_spd", "D6"], ["sg15_10", "D6"]) to pd04spdsg1510_D6D6
            test_name = ""
            for item in range(len(test[0])): # item can only be 1 or 2, 1 is for pseudopotential, 2 is for numerical orbitals
                for axis in test: # axis can be many, for instance InP, it is 2, BaTiO3, it is 3
                    test_name += axis[item].replace("_", "")
                test_name += "_" # add _ between pseudopotential and numerical orbitals
            test_name = test_name[:-1] # remove the last _
            print("prepare test: " + test_name)
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
                if work_status["calculation"]["basis_type"] == "lcao":
                    nao_type_element = test[ie][1]
                # get pseudopotential information and copy to test_status
                test_status["systems"][system][test_name]["pseudopotentials"]["files"][element] = valid_pseudopotentials[element][pseudo_element]["file"]
                test_status["systems"][system][test_name]["pseudopotentials"]["info"][element] = {
                    "kind": valid_pseudopotentials[element][pseudo_element]["kind"],
                    "version": valid_pseudopotentials[element][pseudo_element]["version"],
                    "appendix": valid_pseudopotentials[element][pseudo_element]["appendix"]
                }
                # get numerical orbital information and copy to test_status
                if work_status["calculation"]["basis_type"] == "lcao":
                    test_status["systems"][system][test_name]["numerical_orbitals"]["files"][element] = valid_numerical_orbitals[element][pseudo_element][nao_type_element]["file"]
                    test_status["systems"][system][test_name]["numerical_orbitals"]["info"][element] = {
                        "type": valid_numerical_orbitals[element][pseudo_element][nao_type_element]["type"],
                        "rcut": valid_numerical_orbitals[element][pseudo_element][nao_type_element]["rcut"],
                        "appendix": valid_numerical_orbitals[element][pseudo_element][nao_type_element]["appendix"]
                    }
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

