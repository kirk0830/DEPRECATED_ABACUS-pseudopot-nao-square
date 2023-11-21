from abacus_pseudopotential_square.base.file_io import parse_cif
import abacus_pseudopotential_square.base.atom_in as ai

basic_parameters = {
    "control": {
        "calculation": "scf",
        "restart_mode": "from_scratch",
        "prefix": "abacus",
        "pseudo_dir": "./",
        "outdir": "./",
        "tprnfor": ".true.",
        "tstress": ".true.",
        "wf_collect": ".true.",
        "nstep": 100,
        "verbosity": "high",
    },
    "system": {
        "ibrav": 0,
        "nat": 0,
        "ntyp": 0,
        "ecutwfc": 100,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.02,
        "nspin": 1,
        "starting_magnetization": 0.0,
        "nosym": ".true.",
        "noinv": ".true.",
        "noncolin": ".false.",
        "lspinorb": ".false.",
    },
    "electrons": {
        "electron_maxstep": 100,
        "conv_thr": 1e-7,
        "mixing_beta": 0.7,
        "mixing_mode": "plain",
        "mixing_ndim": 8,
        "diagonalization": "david"
    },
    "ions": {
        "ion_dynamics": "bfgs"
    },
    "cell": {
        "cell_dynamics": "bfgs"
    },
    "k_points": {
        "k_points": [1, 1, 1],
        "k_points_shift": [0, 0, 0],
    }
}

def write_control(other_parameters: dict) -> str:
    """
    write control parameters
    """
    control = basic_parameters["control"]
    if "control" in other_parameters.keys():
        for key, value in other_parameters["control"].items():
            control[key] = value
    control_str = "&CONTROL\n"
    for key, value in control.items():
        control_str += "\t%s = %s,\n"%(key, str(value))
    control_str += "/\n"
    return control_str

def write_system(other_parameters: dict) -> str:
    """
    write system parameters
    """
    system = basic_parameters["system"]
    if "system" in other_parameters.keys():
        for key, value in other_parameters["system"].items():
            system[key] = value
    system_str = "&SYSTEM\n"
    for key, value in system.items():
        system_str += "\t%s = %s,\n"%(key, str(value))
    system_str += "/\n"
    return system_str

def write_electrons(other_parameters: dict) -> str:
    """
    write electrons parameters
    """
    electrons = basic_parameters["electrons"]
    if "electrons" in other_parameters.keys():
        for key, value in other_parameters["electrons"].items():
            electrons[key] = value
    electrons_str = "&ELECTRONS\n"
    for key, value in electrons.items():
        electrons_str += "\t%s = %s,\n"%(key, str(value))
    electrons_str += "/\n"
    return electrons_str

def write_ions(other_parameters: dict) -> str:
    """
    write ions parameters
    """
    ions = basic_parameters["ions"]
    if "ions" in other_parameters.keys():
        for key, value in other_parameters["ions"].items():
            ions[key] = value
    ions_str = "&IONS\n"
    for key, value in ions.items():
        ions_str += "\t%s = %s,\n"%(key, str(value))
    ions_str += "/\n"
    return ions_str

def write_cell(other_parameters: dict) -> str:
    """
    write cell parameters
    """
    cell = basic_parameters["cell"]
    if "cell" in other_parameters.keys():
        for key, value in other_parameters["cell"].items():
            cell[key] = value
    cell_str = "&CELL\n"
    for key, value in cell.items():
        cell_str += "\t%s = %s,\n"%(key, str(value))
    cell_str += "/\n"
    return cell_str

def write_atomic_species(atomic_species: dict) -> str:
    """
    write atomic species
    """
    atomic_species_str = "ATOMIC_SPECIES\n"
    for iatom in range(len(atomic_species["elements"])):
        atomic_species_str += atomic_species["elements"][iatom] + " " + str(atomic_species["mass"][iatom]) + " " + atomic_species["pseudopotentials"][iatom] + "\n"
    atomic_species_str += "\n"
    return atomic_species_str

def write_cell_parameters(lattice: dict) -> str:
    """
    write cell parameters
    """
    cell_parameters_str = "CELL_PARAMETERS (angstrom)\n"
    for i in range(3):
        cell_parameters_str += "\t%12.8f %12.8f %12.8f\n"%(lattice["lattice_vectors"][i][0], lattice["lattice_vectors"][i][1], lattice["lattice_vectors"][i][2])
    cell_parameters_str += "\n"
    return cell_parameters_str

def write_atomic_positions(atoms: dict) -> str:
    """
    write atomic positions
    """
    atomic_positions_str = "ATOMIC_POSITIONS (angstrom)\n"
    for atom in atoms.keys():
        natom = len(atoms[atom]["positions"])
        for iatom in range(natom):
            atomic_positions_str += "\t%s %12.8f %12.8f %12.8f\n"%(atom, atoms[atom]["positions"][iatom][0], atoms[atom]["positions"][iatom][1], atoms[atom]["positions"][iatom][2])
    atomic_positions_str += "\n"
    return atomic_positions_str

def write_k_points(other_parameters: dict) -> str:
    """
    write k points
    """
    k_points = basic_parameters["k_points"]
    if "k_points" in other_parameters.keys():
        for key, value in other_parameters["k_points"].items():
            k_points[key] = value
    k_points_str = "K_POINTS automatic\n"
    k_points_str += "\t%d %d %d %d %d %d\n"%(k_points["k_points"][0], k_points["k_points"][1], k_points["k_points"][2], k_points["k_points_shift"][0], k_points["k_points_shift"][1], k_points["k_points_shift"][2])
    k_points_str += "\n"
    return k_points_str

def cif_to_quantum_espresso(
        cif_file: str,
        pseudopotentials_in: dict,
        other_parameters: dict) -> None:
    qe_filename = cif_file.split(".")[0] + ".in"
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
    return_str = ""
    return_str += write_control(other_parameters)
    return_str += write_system(other_parameters)
    return_str += write_electrons(other_parameters)
    return_str += write_ions(other_parameters)
    return_str += write_cell(other_parameters)
    return_str += write_atomic_species(atomic_species)
    return_str += write_cell_parameters(lattice)
    return_str += write_atomic_positions(atoms)
    return_str += write_k_points(other_parameters)
    with open(qe_filename, 'w') as f:
        f.write(return_str)

if __name__ == "__main__":

    pseudopotentials_in = {
        "files": {
            "La": "La_pseudopot",
            "O": "O_pseudopot",
        },
    }

    cif_to_quantum_espresso("mp-26.cif", pseudopotentials_in, {})