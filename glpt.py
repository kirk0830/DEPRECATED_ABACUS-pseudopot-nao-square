import os
import atom_in as ai

def generate_lcao_pseudopotential_test(test_status: dict):

    if test_status["general"]["software"] == "ABACUS":
        glpt_abacus(test_status)
    else:
        print("Error: software not recognized.")
        exit(1)

def glpt_abacus(test_status: dict):

    os.chdir(test_status["paths"]["work_folder"])
    for system in test_status["systems"].keys():
        for test in test_status["systems"][system].keys():
            folder = "t_" + test + "_" + system
            if os.path.exists(folder):
                print("Warning: " + folder + " already exists, clean.")
                os.system("rm -rf " + folder)
            os.mkdir(folder)
            input_file_path = test_status["paths"]["numerical_orbital"]["work"] + "\\templates\\ABACUS\\INPUT"
            stru_file_path = test_status["paths"]["numerical_orbital"]["work"] + "\\templates\\ABACUS\\STRU_" + system
            kpt_file_path = test_status["paths"]["numerical_orbital"]["work"] + "\\templates\\ABACUS\\KPT"
            # copy INPUT, STRU, KPT to test folder, use unified name STRU on STRU_*
            os.system("cp " + input_file_path + " " + folder + "\\INPUT")
            os.system("cp " + stru_file_path + " " + folder + "\\STRU")
            os.system("cp " + kpt_file_path + " " + folder + "\\KPT")
            # copy pseudopotential and numerical orbital to test folder
            for element in test_status["systems"][system][test]["elements"]:
                pseudopot_path = test_status["paths"]["pseudopotential"]["resources"] + "\\"
                pseudopot_kva = test_status["systems"][system][test]["pseudopotentials"]["info"][element]["kind"] + "_"
                pseudopot_kva += test_status["systems"][system][test]["pseudopotentials"]["info"][element]["version"]
                if test_status["systems"][system][test]["pseudopotentials"]["info"][element]["appendix"] != "":
                    pseudopot_kva += "_" + test_status["systems"][system][test]["pseudopotentials"]["info"][element]["appendix"]
                pseudopot_path += pseudopot_kva + "\\"
                pseudopot_file = test_status["systems"][system][test]["pseudopotentials"]["files"][element]
                print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                os.system("cp " + pseudopot_path + pseudopot_file + " " + folder)
                # swap pseudopotential file name in input file
                os.system("sed -i 's/{}_pseudopot/".format(element) + pseudopot_file + "/g' " + folder + "\\STRU")

                nao_path = test_status["paths"]["numerical_orbital"]["resources"] + "\\"
                nao_path += pseudopot_kva + "\\"
                nao_et = str(ai.get_element_index(element)) + "_" + element + "_"
                nao_et += test_status["systems"][system][test]["numerical_orbitals"]["info"][element]["type"]
                nao_path += nao_et + "\\"
                nao_file = test_status["systems"][system][test]["numerical_orbitals"]["files"][element]
                print("copy " + nao_path + nao_file + " to " + folder)
                os.system("cp " + nao_path + nao_file + " " + folder)
                # swap numerical orbital file name in input file
                os.system("sed -i 's/{}_numerical_orbital/".format(element) + nao_file + "/g' " + folder + "\\STRU")

    os.chdir(test_status["paths"]["root"])

def packup_tests(test_status: dict):

    pass