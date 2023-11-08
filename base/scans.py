import os
import base.atom_in as ai

def scan_elements(system: str) -> list:

    # read letter in system one by one, from the beginning capitalized letter, if meet another, then cut the string
    # and save the previous string to the list, until the end of the string
    elements = []
    element = ""
    for letter in system:
        if letter.isupper():
            if element != "":
                elements.append(element)
            element = letter
        else:
            element += letter
    elements.append(element)
    return elements

def scan_valid_pseudopotentials(work_status: dict):
    """
    Imagin the pseudpotentials are stored in this way:
    in [work_folder]/psedupotentials/resources/

    pd_04_spd/
        As.PD04.PBE.UPF
        ...
    sg15_10/
        As_ONCV_PBE-1.0.upf
        ...

    Output is organized as:
    {
        "As": {
            "pd_04_spd": {
                file: "As.PD04.PBE.UPF",
                kind: "pd",
                version: "04",
                appendix: "",
            },
            "sg15_10": {
                file: "As_ONCV_PBE-1.0.upf",
                kind: "sg15",
                version: "10",
                appendix: "",
            }
            ...
        },
    }
    """
    elements = []
    valid_pseudopotentials = {}
    # collect all elements needed for test
    for system in work_status["systems"].keys():
        for element in scan_elements(system):
            if element not in elements:
                elements.append(element)
    # check-in
    os.chdir(work_status["pseudopotentials"]["pseudopot_paths"]["resources"])
    # then for each element, scan available pseudopotentials
    for element in elements:
        for pseudopot_kind in work_status["pseudopotentials"]["kinds"][element]:
            for version in work_status["pseudopotentials"]["versions"][element]:
                for appendix in work_status["pseudopotentials"]["appendices"][element]:
                    pseudopotential = pseudopot_kind + "_" + version
                    if appendix != "":
                        pseudopotential += "_" + appendix
                    # if find the folder named as pseudopot_kind + "_" + version + "_" + appendix, then it is valid
                    if os.path.isdir(pseudopotential):
                        # check if the pseudopotential file is in the folder, sometimes the file startswith lowercase
                        files = os.listdir(pseudopotential)
                        pseudopot_valid = False
                        for file in files:
                            if file.startswith(element) or file.startswith(element.lower()):
                                pseudopot_valid = True
                                break
                        if not pseudopot_valid:
                            continue
                        # if really find the pseudopotential, then save it to the list
                        if element not in valid_pseudopotentials.keys():
                            valid_pseudopotentials[element] = {}
                        valid_pseudopotentials[element][pseudopotential] = {
                            "file": file,
                            "kind": pseudopot_kind,
                            "version": version,
                            "appendix": appendix,
                        }
                        
    
    return valid_pseudopotentials

def scan_valid_numerical_orbitals(work_status: dict, valid_pseudopotentials: dict):
    """
    Imagine the numerical orbitals are stored in this way:
    in [work_folder]/numerical_orbitals/resources/

    pd_04_spd/
        33_As_DZP/
            As_gga_6au_100Ry_2s2p1d.orb
            As_gga_7au_100Ry_2s2p1d.orb
            ...
        33_As_TZDP/
            As_gga_6au_100Ry_2s2p1d.orb
            ...
    sg15_10/
        33_As_DZP/
            As_gga_6au_100Ry_2s2p1d.orb
            ...
            As_gga_6au_100Ry_2s2p1d_osci_rm.orb
        ...
    ...

    Output is organized as:
    {
        "As": {
            "pd_04_spd: {
                "D6": {
                    "type": "DZP",
                    "rcut": 6,
                    "appendix": "",
                    "file": "As_gga_6au_100Ry_2s2p1d.orb",
                },
                ...
                "T7": {
                    "type": "TZDP",
                    "rcut": 7,
                    "appendix": "",
                    "file": "As_gga_7au_100Ry_3s3p2d.orb",
                },
                ...
            }
            "sg15_10": {
                "D6": {
                    "type": "DZP",
                    "rcut": 6,
                    "appendix": "",
                    "file": "As_gga_6au_100Ry_2s2p1d.orb",
                },
                ...
                "D6_osci_rm": {
                    "type": "DZP",
                    "rcut": 6,
                    "appendix": "osci_rm",
                    "file": "As_gga_6au_100Ry_2s2p1d_osci_rm.orb",
                },
                ...
            }
        },
        ...
    }
    """
    elements = []
    valid_numerical_orbitals = {}
    # collect all elements needed for test
    for system in work_status["systems"].keys():
        for element in scan_elements(system):
            if element not in elements:
                elements.append(element)
    # check-in
    os.chdir(work_status["numerical_orbitals"]["nao_paths"]["resources"])
    print("enter folder: %s" % os.getcwd())
    # after getting valid pseudopotentials, we scan the numerical orbitals according to it
    for element in elements:
        for pseudopotential in valid_pseudopotentials[element].keys():
            if os.path.isdir(pseudopotential):
                os.chdir(pseudopotential)
                folder_header = str(ai.get_element_index(element)) + "_" + element
                for nao_type in work_status["numerical_orbitals"]["types"][element]:
                    folder = folder_header + "_" + nao_type
                    if os.path.isdir(folder):
                        files = os.listdir(folder)
                        for rcut in work_status["numerical_orbitals"]["rcuts"][element]:
                            for appendix in work_status["numerical_orbitals"]["appendices"][element]:
                                nao_startswith = "%s_gga_%sau_100Ry" % (element, rcut)
                                nao_endswith = "%s.orb" % appendix
                                nao_valid = False
                                for file in files:
                                    if file.startswith(nao_startswith) and file.endswith(nao_endswith):
                                        nao_valid = True
                                        break
                                if not nao_valid:
                                    continue
                                if element not in valid_numerical_orbitals.keys():
                                    valid_numerical_orbitals[element] = {}
                                if pseudopotential not in valid_numerical_orbitals[element].keys():
                                    valid_numerical_orbitals[element][pseudopotential] = {}
                                key_nao = nao_type[0] + str(rcut)
                                if appendix != "":
                                    key_nao += "_" + appendix
                                valid_numerical_orbitals[element][pseudopotential][key_nao] = {
                                    "type": nao_type,
                                    "rcut": rcut,
                                    "appendix": appendix,
                                    "file": file,
                                }
                os.chdir("..")
    return valid_numerical_orbitals