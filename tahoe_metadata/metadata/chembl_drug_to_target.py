import requests
import pandas as pd
from tqdm import tqdm

BASE = "https://www.ebi.ac.uk/chembl/api/data"

def get_chembl_id(drug_name):
    url = f"{BASE}/molecule/search.json?q={drug_name}"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return None
    data = r.json()
    if data["page_meta"]["total_count"] == 0:
        return None
    return data["molecules"][0]["molecule_chembl_id"]

def get_drug_mechanisms(chembl_id):
    url = f"{BASE}/drug_mechanism.json?molecule_chembl_id={chembl_id}"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return []
    data = r.json()
    return data["drug_mechanisms"]

def get_target_gene(target_chembl_id):
    url = f"{BASE}/target/{target_chembl_id}.json"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return None
    data = r.json()
    comps = data.get("target_components", [])
    if len(comps) == 0:
        return None
    synonyms = comps[0].get("target_component_synonyms", [])
    for s in synonyms:
        if s["syn_type"] == "GENE_SYMBOL":
            return s["component_synonym"]
    return None

# ========= 主流程 =========
drug_list = pd.read_csv("tahoe_drugs.csv")["drug"].unique()

results = []

for drug in tqdm(drug_list):
    chembl_id = get_chembl_id(drug)
    if chembl_id is None:
        continue

    mechanisms = get_drug_mechanisms(chembl_id)
    for m in mechanisms:
        target_id = m.get("target_chembl_id")
        gene = get_target_gene(target_id)
        if gene:
            results.append({
                "drug": drug,
                "chembl_id": chembl_id,
                "target_chembl_id": target_id,
                "target_gene": gene,
                "mechanism": m.get("mechanism_of_action")
            })

df = pd.DataFrame(results).drop_duplicates()
df.to_csv("drug_target_map_chembl.csv", index=False)

print("Saved drug_target_map_chembl.csv")
