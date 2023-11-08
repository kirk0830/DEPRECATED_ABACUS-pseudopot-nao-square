from mp_api.client import MPRester

def download_cif_from_materials_project(api_key: str):

    with MPRester(api_key) as mpr:
        docs = mpr.summary.search(material_ids=["mp-149", "mp-13", "mp-22526"])

    for doc in docs:
        print(doc)

if __name__ == "__main__":

    download_cif_from_materials_project("5feX86roqGMRzW0kr9G")