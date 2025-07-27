import requests

BASE_URL = "https://api.gdc.cancer.gov"

matched_case_ids = [
    "TCGA-A7-A26I", "TCGA-A7-A0DB", "TCGA-LD-A9QF", "TCGA-AC-A3OD",
    "TCGA-A7-A0CG", "TCGA-A7-A26J", "TCGA-A7-A13G", "TCGA-A7-A13E",
    "TCGA-AC-A3QQ", "TCGA-LD-A74U", "TCGA-A7-A0DC", "TCGA-AC-A2QH",
    "TCGA-E2-A572", "TCGA-A7-A26F", "TCGA-AC-A6IX", "TCGA-A7-A26E",
    "TCGA-AN-A041", "TCGA-A7-A13D"
]

file_uuids = []

for case_id in matched_case_ids:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.submitter_id", "value": [case_id]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Simple Nucleotide Variation"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["VCF"]}},
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["WXS"]}}
        ]
    }
    params = {
        "filters": filters,
        "fields": "file_id,file_name,cases.submitter_id",
        "format": "JSON",
        "size": 100
    }
    response = requests.post(f"{BASE_URL}/files", json=params)
    response.raise_for_status()
    hits = response.json()['data']['hits']
    for hit in hits:
        file_uuids.append((hit['file_id'], hit['file_name'], case_id))

# Print files found:
with open("manifest.txt", "w") as f:
    f.write("id\n") 
    for uuid, fname, case in file_uuids:
        print(f"Case: {case}, File: {fname}, UUID: {uuid}")
        f.write(f"{uuid}\n")