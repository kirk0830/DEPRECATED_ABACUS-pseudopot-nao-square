"""
Interface to external files

Author: Kirk0830
Date: 2023-11-22

Description:
    This file contains functions that can read and write files.

Functions:
    check_folder_availability(folder: str) -> None
    parse_cif(cif_file: str, save_json = False) -> tuple[dict, dict]
"""
import abacus_pseudopotential_square.base.atom_in as ai
import numpy as np
import os
import json

def check_folder_availability(folder: str) -> None:
    """
    check if the folder exists, if not, create one
    """
    if os.path.exists(folder):
        print("Folder check: %s, pass."%folder)
    else:
        print("Folder check: %s, not found, exit."%folder)
        exit(1)

def parse_cif(cif_file: str, save_json = False) -> tuple[dict, dict]:
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
    lattice["lattice_parameters"] = cif["lattice_parameters"]

    if save_json:
        cif_json = cif_file.split(".")[0] + "_origin.json"
        with open(cif_json, 'w') as f:
            json.dump(cif, f, indent = 4)
        print("CIF file original information %s has been converted to JSON file %s."%(cif_file, cif_json))
        lattice_json = cif_file.split(".")[0] + "_lattice.json"
        lattice["lattice_vectors"] = [
            vector_a.tolist(),
            vector_b.tolist(),
            vector_c.tolist(),
        ]
        with open(lattice_json, 'w') as f:
            json.dump(lattice, f, indent = 4)
        print("Lattice information has been saved to JSON file %s."%lattice_json)
        atomic_positions_json = cif_file.split(".")[0] + "_atoms.json"
        for atom in atomic_positions.keys():
            atomic_positions[atom]["positions"] = np.array(atomic_positions[atom]["positions"]).tolist()
        with open(atomic_positions_json, 'w') as f:
            json.dump(atomic_positions, f, indent = 4)
        print("Atomic positions information has been saved to JSON file %s."%atomic_positions_json)
        total_json = cif_file.split(".")[0] + "_total.json"
        total = {
            "cif": cif,
            "lattice": lattice,
            "atomic_positions": atomic_positions,
        }
        with open(total_json, 'w') as f:
            json.dump(total, f, indent = 4)
        print("Total information has been saved to JSON file %s."%total_json)
        
    return atomic_positions, lattice


