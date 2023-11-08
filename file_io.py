import atom_in as ai
import numpy as np

def cif_to_STRU(
        cif_file: str, 
        pseudopotentials_in: dict, 
        numerical_orbitals_in: dict) -> None:
    """
    combine with pseudpotentials and numerical orbitals information (if available)
    to generate STRU file

    pseudopotentials dict should be in the form of:
    "pseudopotentials_in": {
        "files": {
            "In": "In_ONCV_PBE-1.0.upf",
            "N": "N_ONCV_PBE-1.0.upf"
        },
        "info": {
            "In": {
                "kind": "sg15",
                "version": "10",
                "appendix": ""
            },
            "N": {
                "kind": "sg15",
                "version": "10",
                "appendix": ""
            }
        }
    },

    numerical_orbitals dict should be in the form of:
    "numerical_orbitals_in": {
        "files": {
            "In": "In_gga_6au_100Ry_2s2p2d1f.orb",
            "N": "N_gga_6au_100Ry_2s2p1d.orb"
        },
        "info": {
            "In": {
                "type": "DZP",
                "rcut": "6",
                "appendix": ""
            },
            "N": {
                "type": "DZP",
                "rcut": "6",
                "appendix": ""
            }
        }
    }
    """
    STRU_filename = "STRU"
    # convert cif to STRU-organizing dict
    atomic_positions, lattice = parse_cif(cif_file)

    atomic_species = {
        "elements": [],
        "mass": [],
        "pseudopotentials": [],
    }
    for atom in atomic_positions:
        atomic_species["elements"].append(atom["element"])
        atomic_species["mass"].append(ai.get_element_mass(atom["element"]))
        atomic_species["pseudopotentials"].append(
            pseudopotentials_in["files"][atom["element"]],
        )
    
    numerical_orbital = {
        "elements": [],
        "numerical_orbital": [],
    }
    for atom in atomic_positions:
        numerical_orbital["elements"].append(atom["element"])
        numerical_orbital["numerical_orbital"].append(
            numerical_orbitals_in["files"][atom["element"]],
        )

    lattice = {
        "lattice_constant": 0,
        "lattice_vectors": [],
    }

    stru = ""
    stru += write_atomic_species(atomic_species)
    stru += write_numerical_orbital(numerical_orbital)
    stru += write_lattice(lattice)
    stru += write_atomic_positions(atomic_positions)
    with open(STRU_filename, 'w') as f:
        f.write(stru)

def vasp_to_STRU(
        vasp_file: str,
        pseudpotentials: dict,
        numerical_orbitals: dict) -> None:
    raise NotImplementedError("vasp_to_STRU not implemented yet")

def STRU_to_cif(stru: dict):
    raise NotImplementedError("STRU_to_cif not implemented yet")

def parse_cif(cif_file: str) -> dict:
    cif = {
        "lattice_parameters": {
            "a": 0,
            "b": 0,
            "c": 0,
            "alpha": 0,
            "beta": 0,
            "gamma": 0,
        },
        "atomic_position": {},
    }
    with open(cif_file, 'r') as f:
        cif_contents = f.readlines()
    read_atomic_position = False
    for line in cif_contents:
        clean_line = line.strip()
        # import lattice parameters information
        if clean_line.startswith("_cell_length_a"):
            cif["lattice_parameters"]["a"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_length_b"):
            cif["lattice_parameters"]["b"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_length_c"):
            cif["lattice_parameters"]["c"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_alpha"):
            cif["lattice_parameters"]["alpha"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_beta"):
            cif["lattice_parameters"]["beta"] = float(clean_line.split()[1])
        elif clean_line.startswith("_cell_angle_gamma"):
            cif["lattice_parameters"]["gamma"] = float(clean_line.split()[1])
        # import atomic position information
        elif clean_line.startswith("_atom_site_occupancy"):
            read_atomic_position = True
            cif["atomic_position"] = {}
        elif read_atomic_position:
            if clean_line.startswith("_"):
                read_atomic_position = False
            else:
                atom_name = clean_line.split()[1]
                cif["atomic_position"][atom_name] = {
                    "element": clean_line.split()[0],
                    "x": float(clean_line.split()[3]),
                    "y": float(clean_line.split()[4]),
                    "z": float(clean_line.split()[5]),
                    "occupancy": float(clean_line.split()[6]),
                }

    atomic_positions = {}
    for atom in cif["atomic_position"].keys(): # O1, Al1 or something like that
        atom_info = cif["atomic_position"][atom]
        if atom_info["element"] not in atomic_positions.keys():
            atomic_positions[atom_info["element"]] = {
                "n": 1,
                "mag": 0.0,
                "positions": [
                    [atom_info["x"], atom_info["y"], atom_info["z"]],
                ]
            }
        else:
            atomic_positions[atom_info["element"]]["n"] += 1
            atomic_positions[atom_info["element"]]["positions"].append(
                [atom_info["x"], atom_info["y"], atom_info["z"]],
            )
    lattice = {}
    vector_a = cif["lattice_parameters"]["a"] * np.array([1, 0, 0])
    vector_b = cif["lattice_parameters"]["b"] * np.array([np.cos(np.deg2rad(cif["lattice_parameters"]["gamma"])), np.sin(np.deg2rad(cif["lattice_parameters"]["gamma"])), 0])
    vector_c_1 = np.cos(np.deg2rad(cif["lattice_parameters"]["beta"]))
    vector_c_2 = (np.cos(np.deg2rad(cif["lattice_parameters"]["alpha"])) - np.cos(np.deg2rad(cif["lattice_parameters"]["gamma"])) * np.cos(np.deg2rad(cif["lattice_parameters"]["beta"]))) / np.sin(np.deg2rad(cif["lattice_parameters"]["gamma"]))
    vector_c_3 = np.sqrt(1 - vector_c_1**2 - vector_c_2**2)
    vector_c = cif["lattice_parameters"]["c"] * np.array([vector_c_1, vector_c_2, vector_c_3])
    lattice["lattice_constant"] = 1.889716
    lattice["lattice_vectors"] = [
        vector_a,
        vector_b,
        vector_c,
    ]
    return atomic_positions, lattice

def write_atomic_species(atomic_species: dict) -> str:
    
    print(atomic_species)
    return_str = "ATOMIC_SPECIES\n"
    for atom_symbol in atomic_species.keys():
        atom = atomic_species[atom_symbol]
        return_str += atom["element"] + " " + atom["mass"] + " " + atom["pseudopotential"] + "\n"
    return_str += "\n"
    return return_str

def write_numerical_orbital(numerical_orbital: dict) -> str:
    
    return_str = "NUMERICAL_ORBITAL\n"
    for atom in numerical_orbital:
        return_str += atom["element"] + atom["numerical_orbital"] + "\n"
    return_str += "\n"
    return return_str

def write_lattice(lattice: dict) -> str:
    
    return_str = "LATTICE_CONSTANT\n%10.20f\n\nLATTICE_VECTORS\n"%lattice["lattice_constant"]
    for vector in lattice["lattice_vectors"]:
        return_str += "%10.20f %10.20f %10.20f\n"%(vector[0], vector[1], vector[2])
    return return_str

def write_atomic_positions(atomic_positions: dict) -> str:
    
    return_str = "ATOMIC_POSITIONS\nDirect\n"
    for atom in atomic_positions:
        return_str += atom["element"] + "\n"
        if "mag" in atom:
            return_str += str(atom["mag"]) + "\n"
        else:
            return_str += "0.0\n"
        for position in atom["positions"]:
            return_str += "%10.20f %10.20f %10.20f 1 1 1\n"%(position[0], position[1], position[2])
    return_str += "\n"
    return return_str

