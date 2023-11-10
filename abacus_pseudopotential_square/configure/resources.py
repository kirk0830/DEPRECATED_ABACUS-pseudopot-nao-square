from mp_api.client import MPRester

def download_cif_from_materials_project(api_key: str):

    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(elements=["Si", "O"], 
                                band_gap=(0.5, 1.0))
    for doc in docs:
        print(doc)

if __name__ == "__main__":

    download_cif_from_materials_project("wV1HUdmgESPVgSmQj5cc8WvttCO8NTHp")