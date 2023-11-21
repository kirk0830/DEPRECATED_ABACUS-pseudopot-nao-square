from abacus_pseudopotential_square.base.file_io import parse_cif
import abacus_pseudopotential_square.base.atom_in as ai

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
    atoms, lattice = parse_cif(cif_file)

    atomic_species = {
        "elements": [],
        "mass": [],
        "pseudopotentials": [],
    }
    for element in atoms.keys():
        atomic_species["elements"].append(element)
        atomic_species["mass"].append(ai.get_element_mass(element))
        atomic_species["pseudopotentials"].append(
            pseudopotentials_in["files"][element],
        )
    
    numerical_orbital = {
        "elements": [],
        "numerical_orbitals": [],
    }
    for element in atoms.keys():
        numerical_orbital["elements"].append(element)
        numerical_orbital["numerical_orbitals"].append(
            numerical_orbitals_in["files"][element],
        )

    lattice = {
        "lattice_constant": 1.889716,
        "lattice_vectors": [
            lattice["lattice_vectors"][0].tolist(),
            lattice["lattice_vectors"][1].tolist(),
            lattice["lattice_vectors"][2].tolist()
        ],
    }

    stru = ""
    stru += write_atomic_species(atomic_species)
    stru += write_numerical_orbital(numerical_orbital)
    stru += write_lattice(lattice)
    stru += write_atomic_positions(atoms)
    with open(STRU_filename, 'w') as f:
        f.write(stru)

def STRU_to_cif(stru: dict):
    raise NotImplementedError("STRU_to_cif not implemented yet")

def write_atomic_species(atomic_species: dict) -> str:
    
    print(atomic_species)
    return_str = "ATOMIC_SPECIES\n"
    for iatom in range(len(atomic_species["elements"])):
        return_str += atomic_species["elements"][iatom] + " " + str(atomic_species["mass"][iatom]) + " " + atomic_species["pseudopotentials"][iatom] + "\n"
    return_str += "\n"
    return return_str

def write_numerical_orbital(numerical_orbital: dict) -> str:
    
    print(numerical_orbital)
    return_str = "NUMERICAL_ORBITAL\n"
    for iatom in range(len(numerical_orbital["elements"])):
        return_str += numerical_orbital["elements"][iatom] + " " + numerical_orbital["numerical_orbitals"][iatom] + "\n"
    return_str += "\n"
    return return_str

def write_lattice(lattice: dict) -> str:
    
    return_str = "LATTICE_CONSTANT\n%20.10f\n\nLATTICE_VECTORS\n"%lattice["lattice_constant"]
    for vector in lattice["lattice_vectors"]:
        return_str += "%20.10f %20.10f %20.10f\n"%(vector[0], vector[1], vector[2])
    return_str += "\n"
    return return_str

def write_atomic_positions(atoms: dict) -> str:
    
    return_str = "ATOMIC_POSITIONS\nDirect\n"
    for element in atoms.keys():
        return_str += element + "\n"
        if "mag" in atoms[element].keys():
            return_str += str(atoms[element]["mag"]) + "\n"
        else:
            return_str += "0.0\n"
        for position in atoms[element]["positions"]:
            return_str += "%20.10f %20.10f %20.10f 1 1 1\n"%(position[0], position[1], position[2])
    return_str += "\n"
    return return_str
