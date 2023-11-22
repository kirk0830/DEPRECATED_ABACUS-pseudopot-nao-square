import abacus_pseudopotential_square.analysis.fold_band as fb
import os
import matplotlib.pyplot as plt
import numpy as np
import time
import json

def analysis_band_gap(work_status: dict, test_status: dict) -> None:

    if test_status["global"]["software"] == "Quantum ESPRESSO":
        log_status = read_energies_quantum_espresso(work_status, test_status)
        if test_status["draw_statistics"]:
            statistics(work_status, test_status, log_status)
        if test_status["dos"]["draw"]:
            log_status = calculate_dos(log_status, delta_e = test_status["dos"]["delta_e"], smear = test_status["dos"]["smear"])
            draw_dos(log_status)
    elif test_status["global"]["software"] == "ABACUS":
        pass
    else:
        print("abg>> Error: software not recognized.")
        exit(1)

def read_energies_quantum_espresso(work_status: dict, test_status: dict):
    results_path = test_status["paths"]["pseudopotential"]["work"]
    results_path += "\\results"
    os.chdir(results_path)

    log_status = {}
    for system in test_status["systems"].keys():
        log_status[system] = {}
        for test in test_status["systems"][system].keys():
            folder = "t_" + test
            if os.path.isdir(folder):
                os.chdir(folder)
                print("abg>> Analyzing " + folder + "...")
                if os.path.exists("CRASH"):
                    print("abg>> CRASH detected in folder ", folder)
                    # for we use absolute path, we do not need to change directory back
                    continue
                log_status[system][test] = {
                    "pseudopotential": {
                        "pseudopot_kind": "",
                        "version": "",
                        "appendix": ""
                    },
                    "dos": {
                        "window_width": work_status["analysis"]["dos"]["window_width"],
                        "emin": -1, "emax": -1,
                        "delta_e": work_status["analysis"]["dos"]["delta_e"],
                        "smear": work_status["analysis"]["dos"]["smear"]
                    },
                    "energies": {
                        "band_gap": -1, "band_energies": [],
                        "efermi": -1, "homo": -1, "lumo": -1,
                    },
                    "nelec": -1, "nband": -1
                }
                concatenated_pseudopot_kind = ""
                concatenated_pseudopot_version = ""
                concatenated_pseudopot_appendix = ""
                pseudopot_informations = test_status["systems"][system]["pseudopotentials"]["info"]
                for element in pseudopot_informations.keys():
                    concatenated_pseudopot_kind += pseudopot_informations[element]["kind"]
                    concatenated_pseudopot_version += pseudopot_informations[element]["version"]
                    concatenated_pseudopot_appendix += pseudopot_informations[element]["appendix"]
                log_status[system][test]["pseudopotential"]["pseudopot_kind"] = concatenated_pseudopot_kind
                log_status[system][test]["pseudopotential"]["version"] = concatenated_pseudopot_version
                log_status[system][test]["pseudopotential"]["appendix"] = concatenated_pseudopot_appendix
                
                if os.path.exists("scf.log"):
                    print("abg>> parsing ", folder)
                    band_energies, nelec, nband, efermi = fb.parse_qe_output("scf.log")
                    log_status[system][test]["energies"]["band_energies"] = band_energies
                    log_status[system][test]["nelec"] = nelec
                    log_status[system][test]["nband"] = nband
                    log_status[system][test]["energies"]["efermi"] = efermi
                    print("abg>> nelec = ", nelec, " nband = ", nband, " efermi = ", efermi)
                    band_gap, band_energies, homo, lumo = fb.calculate_band_gap(band_energies, nelec, efermi, True)
                    log_status[system][test]["energies"]["band_gap"] = band_gap
                    log_status[system][test]["energies"]["band_energies"] = band_energies
                    log_status[system][test]["energies"]["homo"] = homo
                    log_status[system][test]["energies"]["lumo"] = lumo

                    print("abg>> band gap calculated to be: ", round(band_gap, 4), " ev")
                os.chdir("..")
                print("abg>> back to folder ", os.getcwd())
    os.chdir(test_status["paths"]["root"])
    return log_status

def count_test_kinds(log_status: dict) -> int:
    
    return len(log_status[list(log_status.keys())[0]].keys())

def calculate_dos(log_status: dict, delta_e = 0.01, smear = 0.02) -> dict:
    for system in log_status.keys():
        for test in log_status[system].keys():
            print("abg>> calculating density of states for system {0}, test {1}".format(system, test))
            emin = log_status[system][test]["energies"]["efermi"] - log_status[system][test]["dos"]["window_width"]/2
            emax = log_status[system][test]["energies"]["efermi"] + log_status[system][test]["dos"]["window_width"]/2
            log_status[system][test]["dos"]["emin"] = emin
            log_status[system][test]["dos"]["emax"] = emax
            log_status[system][test]["dos"]["data"] = fb.calculate_dos(band_energies=log_status[system][test]["energies"]["band_energies"], 
                                                                              emin = emin, emax = emax,
                                                                              delta_e = delta_e, smear = smear)
    return log_status

def draw_dos(log_status: dict) -> None:

    print("abg>> drawing density of states")
    plt_nrow = count_test_kinds(log_status)
    plt_ncol = len(log_status.keys())
    fig, axs = plt.subplots(plt_nrow, plt_ncol, figsize=(8, 15), squeeze=False)

    for icol, system in enumerate(log_status.keys()):
        for irow, test in enumerate(log_status[system].keys()):
            x = np.linspace(
                log_status[system][test]["dos"]["emin"],
                log_status[system][test]["dos"]["emax"],
                len(log_status[system][test]["dos"]["data"])
                )
            axs[irow, icol].plot(x, log_status[system][test]["dos"]["data"])
            axs[irow, icol].axvline(x = log_status[system][test]["energies"]["efermi"], color = 'b', linestyle='-.')
            axs[irow, icol].axvline(x = log_status[system][test]["energies"]["homo"], color = 'r')
            axs[irow, icol].axvline(x = log_status[system][test]["energies"]["lumo"], color = 'r')
            axs[irow, icol].set_xlabel("energy (eV)")
            axs[irow, icol].set_ylabel("DOS")

            title = "System: " + system 
            title += ", Pseudopotential: " + log_status[system][test]["pseudopotential"]["pseudopot_kind"]
            title += ", Version: " + log_status[system][test]["pseudopotential"]["version"]
            if log_status[system][test]["pseudopotential"]["appendix"] != "":
                title += log_status[system][test]["pseudopotential"]["appendix"]
            
            axs[irow, icol].set_xlim(
                log_status[system][test]["dos"]["emin"],
                log_status[system][test]["dos"]["emax"]
            )
            axs[irow, icol].set_ylim(0, 2)
            axs[irow, icol].fill_between(
                x,
                log_status[system][test]["dos"]["data"].flatten(),
                y2 = -1,
                where=((x >= log_status[system][test]["dos"]["emin"]) & (x <= log_status[system][test]["energies"]["efermi"])),
                alpha = 0.5
            )
            axs[irow, icol].text(
                0.95,
                0.95, 
                'Band gap = {} eV'.format(
                    round(
                        log_status[system][test]["energies"]["band_gap"],
                        4)
                    ), 
                    transform=axs[irow, icol].transAxes, 
                    ha='right', va='top'
            )
            axs[irow, icol].set_title(title)
                
    plt.subplots_adjust(hspace=0.5, wspace=0.5, top=0.9, bottom=0.1, left=0.1, right=0.9)
    plt.tight_layout()
    # draw
    plt.show()

def statistics(work_status: dict, test_status: dict, log_status: dict) -> None:

    print("abg>> To do precision statistics, experimental values are needed for all systems.")
    ncol = 2
    fig, axs = plt.subplots(int(np.ceil(len(work_status["systems"].keys())/ncol)), ncol, figsize=(8, 8), squeeze=False)
    
    add_reference_line = False
    num_reference_values = 0

    for system in work_status["systems"].keys():
        if work_status["systems"][system]["experimental_value"]["band_gap"] >= 0:
            num_reference_values += 1
    if num_reference_values == len(log_status.keys()):
        add_reference_line = True
    
    band_gaps = []
    band_gaps_abs_error = {}
    band_gaps_rank = []
    for icol, system in enumerate(log_status.keys()):
        band_gaps.append([])
        pseudopotentials = []
        for irow, test in enumerate(log_status[system].keys()):
            band_gap = log_status[system][test]["energies"]["band_gap"]
            band_gaps[-1].append(band_gap)
            
            title = log_status[system][test]["pseudopotential"]["pseudopot_kind"]
            if log_status[system][test]["pseudopotential"]["version"] != "":
                title += "_" + log_status[system][test]["pseudopotential"]["version"]
            if log_status[system][test]["pseudopotential"]["appendix"] != "":
                title += log_status[system][test]["pseudopotential"]["appendix"]
            pseudopotentials.append(title)

            if add_reference_line:
                if title not in band_gaps_abs_error.keys():
                    band_gaps_abs_error[title] = 0
                band_gaps_abs_error[title] += abs(band_gap - test_status["systems"][system]["experimental_value"])
        
        _irow = int(icol/ncol)
        _icol = icol%ncol
        if add_reference_line:
            axs[_irow][_icol].hlines(
                y=test_status["systems"][system]["experimental_value"]["band_gap"], 
                xmin=-0.5, xmax=len(band_gaps[-1])+0.5, 
                linestyle='--', color = 'red')
            axs[_irow][_icol].text(
                len(band_gaps[-1])+1.0, 
                test_status["systems"][system]["experimental_value"]["band_gap"], 
                "Experimental value", 
                ha='left', va='center', 
                color='red', fontsize=12)
            # readjust to ensure reference line always echo
            axs[_irow][_icol].set_ylim(
                bottom = 0, top = max(max(band_gaps[-1]), test_status["systems"][system]["experimental_value"]["band_gap"]) + 0.5)
        axs[_irow][_icol].bar(pseudopotentials, band_gaps[-1])
        title = "System: " + system
        if title != "":
            axs[_irow][_icol].set_title(title)
        axs[_irow][_icol].set_xlabel("Pseudopotentials")
        axs[_irow][_icol].set_ylabel("Band gap (eV)")
        axs[_irow][_icol].tick_params(axis='x', rotation=45)
    # resize subplots
    plt.subplots_adjust(hspace=0.5, wspace=0.5, top=0.9, bottom=0.1, left=0.1, right=0.9)
    # draw
    plt.show()

    if add_reference_line:
        _data = []
        for key in band_gaps_abs_error.keys():
            _data.append((band_gaps_abs_error[key], key))
        dtype = [('error', float), ('pseudopotential name', 'S20')]
        _sorted = np.sort(np.array(_data, dtype), order = 'error')
        _score = [1/sorted[0] for sorted in _sorted]
        _name = [sorted[1] for sorted in _sorted]
        plt.bar(_name, _score)

        plt.xticks(rotation=45)
        plt.xlabel("Pseudopotentials", fontsize=12)
        plt.ylabel("Score (1/eV)", fontsize=12)
        plt.title("Band gap accuracy score", fontsize=12)
        plt.show()

def capitalize_folder_names():

    """ this is a temperary function whose function can be anything """
    # check in folder ./results and change all 'insb' or 'inas' in folders' name to 'InSb' or 'InAs'
    os.chdir("./results")
    for folder in os.listdir():
        if folder.startswith('t_'):
            if folder.endswith('inas'):
                os.rename(folder, folder.replace('inas', 'InAs'))
            elif folder.endswith('insb'):
                os.rename(folder, folder.replace('insb', 'InSb'))
    os.chdir('../')

def preprocess_2d(tests: dict):

    exist_tests_1d_startswith = {}
    folder_names = os.listdir(tests["work_folder"])
    folder_names = [folder for folder in folder_names if folder.startswith('t_')]
    for pseudopot_kind in tests["pseudopot_kinds"]:
        for version in tests["versions"]:
            for appendix in tests["appendices"]:
                folder_startswith = "t_"+pseudopot_kind+"_"+version
                if appendix != "":
                    folder_startswith += "_" + appendix
                folder_startswith += "_" # to ensure the system name is not included
                for folder_name in folder_names:
                    if folder_name.startswith(folder_startswith):
                        if pseudopot_kind not in exist_tests_1d_startswith.keys():
                            exist_tests_1d_startswith[pseudopot_kind] = {}
                        if version not in exist_tests_1d_startswith[pseudopot_kind].keys():
                            exist_tests_1d_startswith[pseudopot_kind][version] = []
                        if appendix not in exist_tests_1d_startswith[pseudopot_kind][version]:
                            exist_tests_1d_startswith[pseudopot_kind][version].append(appendix)
                        print("abg>> " + "-"*50)
                        print("abg>> Preprocess_2d\nFor pseudopot_kind = ", pseudopot_kind, ", \nversion = ", version, ", \nappendix = ", appendix, ", \nfind system = ", folder_name.split("_")[-1])
    if exist_tests_1d_startswith != {}:
        print("abg>> " + "-"*50)

    # then combine all possible combinations
    tests_2d = {}
    # and record the mapping from 2d to 1d
    mapping_tests_2d_to_1d = {}

    pseudopot_kinds_2d = []
    versions_2d = []
    appendices_2d = []
    for pseudopot_kind_i in exist_tests_1d_startswith.keys():
        for pseudopot_kind_j in exist_tests_1d_startswith.keys():
            # concatenate pseudopot_kind_i and pseudopot_kind_j if they are different
            pseudopot_kind = pseudopot_kind_i
            if pseudopot_kind_i != pseudopot_kind_j:
                pseudopot_kind += pseudopot_kind_j
            if pseudopot_kind not in pseudopot_kinds_2d:
                pseudopot_kinds_2d.append(pseudopot_kind)

            for version_i in exist_tests_1d_startswith[pseudopot_kind_i].keys():
                for version_j in exist_tests_1d_startswith[pseudopot_kind_j].keys():
                    # concatenate version_i and version_j if they are different
                    version = version_i
                    if version_i != version_j:
                        version += version_j
                    if version not in versions_2d:
                        versions_2d.append(version)

                    for appendix_i in exist_tests_1d_startswith[pseudopot_kind_i][version_i]:
                        for appendix_j in exist_tests_1d_startswith[pseudopot_kind_j][version_j]:
                            # concatenate appendix_i and appendix_j if they are different
                            appendix = appendix_i
                            if appendix_i != appendix_j:
                                appendix += appendix_j
                            if appendix not in appendices_2d:
                                appendices_2d.append(appendix)

                            concatenated_test = pseudopot_kind+"_"+version
                            if appendix != "":
                                concatenated_test += "_"+appendix
                            test_i = pseudopot_kind_i+"_"+version_i
                            if appendix_i != "":
                                test_i += "_"+appendix_i
                            test_j = pseudopot_kind_j+"_"+version_j
                            if appendix_j != "":
                                test_j += "_"+appendix_j
                            if concatenated_test not in mapping_tests_2d_to_1d.keys():
                                # exclude the case like sg15_11_fr being decompsed into (sg15_11_fr, sg15_11), instead, 
                                # it is (sg15_11_fr, sg15_11_fr)
                                if test_i.startswith(test_j) or test_j.startswith(test_i):
                                    test = len(test_i) > len(test_j) and test_i or test_j
                                    mapping_tests_2d_to_1d[concatenated_test] = [test, test]
                                else:
                                    mapping_tests_2d_to_1d[concatenated_test] = [test_i, test_j]

    tests_2d["pseudopot_kinds"] = pseudopot_kinds_2d
    tests_2d["versions"] = versions_2d
    tests_2d["appendices"] = appendices_2d
    # copy all other information from tests
    tests_2d["work_folder"] = tests["work_folder"]
    tests_2d["dos"] = tests["dos"]
    tests_2d["draw_statistics"] = tests["draw_statistics"]
    tests_2d["systems"] = tests["systems"]

    # export the test_2d to json and name the file as "test_status" with time stamp
    os.chdir(tests["work_folder"])
    with open("test_status_2d_"+time.strftime("%Y%m%d-%H%M%S")+".json", "w") as f:
        json.dump(tests_2d, f, indent=4)
    with open("test_status_exists_startswith_1d_"+time.strftime("%Y%m%d-%H%M%S")+".json", "w") as f:
        json.dump(exist_tests_1d_startswith, f, indent=4)
    os.chdir("..")

    return exist_tests_1d_startswith, tests_2d, mapping_tests_2d_to_1d

def postprocess_2d(exist_tests_1d_startswith: dict, tests_2d: dict, mapping_tests_2d_to_1d: dict, log_status: dict) -> dict:

    """ log_status cannot distinguish between different combinations of pseudopotentials diagonal and off-diagonal, use tests dict to help it """
    # because the property is of interest is only band gap, so only band gap is needed
    # the return dict(s) will organize information in the following way:
    # {
    #     axis_i: [pseudopot_kind_version_appendix, ...],
    #     axis_j: [pseudopot_kind_version_appendix, ...],
    #     data: {
    #     "system_1": [[], [], ...]],
    #     "system_2": [[], [], ...]],
    #     ...
    # }
    # }
    
    # count number of diagonal elements, create axis_i and axis_j
    num_diagonal = 0
    axis_i = []
    axis_j = []
    for pseudopot_kind in exist_tests_1d_startswith.keys():
        for version in exist_tests_1d_startswith[pseudopot_kind].keys():
            for appendix in exist_tests_1d_startswith[pseudopot_kind][version]:
                num_diagonal += 1
                test = pseudopot_kind+"_"+version
                if appendix != "":
                    test += "_"+appendix
                axis_i.append(test)
                axis_j.append(test)
    # create return dict
    return_dict = {
        "axis_i": axis_i,
        "axis_j": axis_j,
        "data": {
            "values": {},
            "availabilities": {}
        }
    }
    # fill in return dict
    for system in log_status.keys():
        data_system = np.zeros(shape = (num_diagonal, num_diagonal))
        availability = np.zeros(shape = (num_diagonal, num_diagonal))
        for test in log_status[system].keys():
            pseudopot_kind_version_appendix = log_status[system][test]["pseudopotential"]["pseudopot_kind"]
            if log_status[system][test]["pseudopotential"]["version"] != "":
                pseudopot_kind_version_appendix += "_"+log_status[system][test]["pseudopotential"]["version"]
            if log_status[system][test]["pseudopotential"]["appendix"] != "":
                pseudopot_kind_version_appendix += "_"+log_status[system][test]["pseudopotential"]["appendix"]
            test_i = mapping_tests_2d_to_1d[pseudopot_kind_version_appendix][0]
            test_j = mapping_tests_2d_to_1d[pseudopot_kind_version_appendix][1]
            i = axis_i.index(test_i)
            j = axis_j.index(test_j)
            data_system[i][j] = log_status[system][test]["energies"]["band_gap"] - tests_2d["systems"][system]["experimental_value"]
            availability[i][j] = 1

        return_dict["data"]["values"][system] = data_system.tolist()
        return_dict["data"]["availabilities"][system] = availability.tolist()

    return return_dict

def plot_band_gap_2d(mode: str, data: dict, ncol: int, color_map = "winter") -> None:

    availabibities = data["data"]["availabilities"]
    values = data["data"]["values"]
    
    summary = np.zeros_like(data["data"]["values"][list(data["data"]["values"].keys())[0]])

    nrows = int(np.ceil(len(values.keys())/ncol))
    fig, axs = plt.subplots(nrows, ncol, figsize=(8, 15), squeeze=False)
    # for each system/subplot, for each data point [i][j], if availability[i][j] == 0, means no data, then set values[i][j] = 0 but not display
    for icol, system in enumerate(values.keys()):
        for irow, test_i in enumerate(data["axis_i"]):
            for jcol, test_j in enumerate(data["axis_j"]):
                # for data not available, set value to empty such that imshow will not display it
                if availabibities[system][irow][jcol] == 0:
                    values[system][irow][jcol] = -1000
                    summary[irow][jcol] = -1000
                else:
                    summary[irow][jcol] += np.abs(values[system][irow][jcol])
        if mode == "imshow":
            # use imshow
            axs[int(icol/ncol)][icol%ncol].imshow(values[system], cmap = color_map, vmin = -1, vmax = 1)
            # add colorbar
            cbar = axs[int(icol/ncol)][icol%ncol].figure.colorbar(axs[int(icol/ncol)][icol%ncol].images[0], ax=axs[int(icol/ncol)][icol%ncol])
            cbar.ax.set_ylabel('Band gap error (eV)', rotation=-90, va="bottom")
            axs[int(icol/ncol)][icol%ncol].set_xlabel("Pseudopotential kinds of atom {} ->".format(system[2:]))
            axs[int(icol/ncol)][icol%ncol].set_ylabel("Pseudopotential kinds of atom In ->")
        elif mode == "contour_f":
            vlevel = np.linspace(-1, 1, 50)
            axs[int(icol/ncol)][icol%ncol].contourf(values[system], vlevel, cmap = color_map)
        # set xticks and yticks
        axs[int(icol/ncol)][icol%ncol].set_xticks(np.arange(len(data["axis_j"])))
        axs[int(icol/ncol)][icol%ncol].set_yticks(np.arange(len(data["axis_i"])))
        # set xlabel and ylabel
        axs[int(icol/ncol)][icol%ncol].set_xticklabels(data["axis_j"])
        axs[int(icol/ncol)][icol%ncol].set_yticklabels(data["axis_i"])
        # set title of subplot
        axs[int(icol/ncol)][icol%ncol].set_title(system)
        # rotate xticks
        plt.setp(axs[int(icol/ncol)][icol%ncol].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        # add grid
        #axs[int(icol/ncol)][icol%ncol].grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.5)
    
    # resize subplots
    plt.subplots_adjust(hspace=0.5, wspace=0.5, top=0.9, bottom=0.1, left=0.1, right=0.9)
    # draw
    plt.show()

    # plot summary
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(summary, cmap = color_map, vmin = 0, vmax = 1)
    # add colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Band gap absolute error sum (eV)', rotation=-90, va="bottom")
    # set xticks and yticks
    ax.set_xticks(np.arange(len(data["axis_j"])))
    ax.set_yticks(np.arange(len(data["axis_i"])))
    # set xlabel and ylabel
    ax.set_xticklabels(data["axis_j"])
    ax.set_yticklabels(data["axis_i"])
    # rotate xticks
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    # add grid
    #ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.5)
    # set title of subplot
    plt.title("Summary")
    # draw
    plt.show()
