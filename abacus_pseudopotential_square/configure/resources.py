from mp_api.client import MPRester

def single_element_structures(api_key: str, element: str, num_cif = 1):
    
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(elements=[element], num_elements=1)
        
    mpi_ids = [doc.material_id for doc in docs]

    num_cif_downloaded = 0
    for mpi_id in mpi_ids:
        with MPRester(api_key) as mpr:
            structure = mpr.get_structure_by_material_id(mpi_id)
        structure.to(filename=f"{mpi_id}.cif", fmt="cif")
        num_cif_downloaded += 1
        if num_cif_downloaded == num_cif:
            break

if __name__ == "__main__":

    api_key = "wV1HUdmgESPVgSmQj5cc8WvttCO8NTHp"
    elements = ['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']
    for element in elements:
        single_element_structures(api_key, element, num_cif = 1)
    