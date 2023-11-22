"""
Online interface to Materials Project database

Author: Kirk0830
Date: 2023-11-22

Description:
    This file contains functions that can download data from Materials Project database.
    website: https://materialsproject.org/

Functions:
    single_element_structures(api_key: str, element: str, num_cif = 1) -> list[str]
"""

from mp_api.client import MPRester

def single_element_structures(api_key: str, element: str, num_cif = 1):
    """
    down load single element structures from materials project

    Args:

    :param api_key: Materials Project API key
    :param element: element name
    :param num_cif: number of cif files to download, default is 1

    Returns:

    :return: cif_filenames: list of cif filenames
    """
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(elements=[element], num_elements=1)
        
    mpi_ids = [doc.material_id for doc in docs]

    num_cif_downloaded = 0
    cif_filenames = []

    for mpi_id in mpi_ids:
        with MPRester(api_key) as mpr:
            structure = mpr.get_structure_by_material_id(mpi_id)
        structure.to(filename=f"{mpi_id}.cif", fmt="cif")
        num_cif_downloaded += 1
        cif_filenames.append(f"{mpi_id}.cif")
        if num_cif_downloaded == num_cif:
            break

    return cif_filenames

if __name__ == "__main__":

    api_key = "wV1HUdmgESPVgSmQj5cc8WvttCO8NTHp"
    elements = ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']
    for element in elements:
        single_element_structures(api_key, element, num_cif = 1)
    