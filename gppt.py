import os

def generate_pw_pseudopotential_test(test_status: dict):

    if test_status["general"]["software"] == "Quantum ESPRESSO":
        gppt_quantum_espresso(test_status)
    elif test_status["general"]["software"] == "ABACUS":
        gppt_abacus(test_status)
    else:
        print("Error: software not recognized.")
        exit(1)

def gppt_quantum_espresso(test_status: dict):

    os.chdir(test_status["paths"]["work_folder"])
    for functional in test_status["dft_settings"]["functionals"]:
        for system in test_status["systems"].keys():
            for test in test_status["systems"][system].keys():
                folder = "t_" + functional + "_" + test + "_" + system
                if os.path.exists(folder):
                    print("Warning: " + folder + " already exists, clean.")
                    os.system("rm -rf " + folder)
                os.mkdir(folder)
                input_file_path = test_status["paths"]["pseudopotential"]["work"] + "\\templates\\Quantum_ESPRESSO\\template_" + system + ".in"
                # copy input file to test folder, use unified name scf.in
                os.system("cp " + input_file_path + " " + folder + "\\scf.in")
                os.system("sed -i 's/functional_to_test/" + functional + "/g' " + folder + "\\scf.in")
                # copy pseudopotential to test folder
                for element in test_status["systems"][system][test]["elements"]:
                    pseudopot_path = test_status["paths"]["pseudopotential"]["resources"] + "\\"
                    pseudopot_path += test_status["systems"][system][test]["pseudopotentials"]["info"][element]["kind"] + "_"
                    pseudopot_path += test_status["systems"][system][test]["pseudopotentials"]["info"][element]["version"]
                    if test_status["systems"][system][test]["pseudopotentials"]["info"][element]["appendix"] != "":
                        pseudopot_path += "_" + test_status["systems"][system][test]["pseudopotentials"]["info"][element]["appendix"]
                    pseudopot_path += "\\"
                    pseudopot_file = test_status["systems"][system][test]["pseudopotentials"]["files"][element]
                    print("copy " + pseudopot_path + pseudopot_file + " to " + folder)
                    os.system("cp " + pseudopot_path + pseudopot_file + " " + folder)
                    # swap pseudopotential file name in input file
                    os.system("sed -i 's/{}_pseudopot/".format(element) + pseudopot_file + "/g' " + folder + "\\scf.in")
    os.chdir(test_status["paths"]["root"])
    print_str = """Generation Done.\n
To run Quantum ESPRESSO tests on ABACUS Test, use the following parameters:\n
rundft: mpirun -np 16 pw.x -i scf.in | tee scf.log; rm -rf pwscf.save\n
Bohrium image: registry.dp.tech/dptech/prod-471/abacus-vasp-qe:20230116\n
Bohrium_machine_type: c32_m128_cpu\n
Bohrium_platform: ali
"""
    print(print_str)
def gppt_abacus(test_status: dict):

    raise NotImplementedError("ABACUS pseudopotential test is not implemented yet.")
