
import requests
import pandas as pd
from collections import defaultdict

# Define the API endpoint
BASE_URL = "https://api.gdc.cancer.gov"

# Step 1: Query case metadata
params = {
    "filters": {
        "op": "in",
        "content": {
            "field": "project.project_id",
            "value": ["TCGA-BRCA"]
        }
    },
    "fields": "case_id,submitter_id,samples.sample_type,samples.portions.portion_number,samples.portions.analytes.analyte_type",
    "format": "JSON",
    "size": 10000
}

response = requests.post(f"{BASE_URL}/cases", json=params)
response.raise_for_status()
data = response.json()

# Step 2: Parse case-level metadata
case_list = []

for hit in data["data"]["hits"]:
    case_id = hit.get("submitter_id", hit.get("case_id"))
    for sample in hit.get("samples", []):
        sample_type = sample.get("sample_type")
        if sample_type != "Primary Tumor":
            continue
        for portion in sample.get("portions", []):
            portion_number = portion.get("portion_number")
            if portion_number is None:
                continue
            try:
                pnum = int(portion_number)
                case_list.append({
                    "case_id": case_id,
                    "portion_number": pnum
                })
            except ValueError:
                continue

# Step 3: Convert to DataFrame
df = pd.DataFrame(case_list)

# Step 4: Identify matched FFPE and frozen
# Frozen if portion_number < 10, FFPE if portion_number >= 10
case_flags = defaultdict(lambda: {"frozen": False, "ffpe": False})

for _, row in df.iterrows():
    cid = row["case_id"]
    pnum = row["portion_number"]
    if pnum < 10:
        case_flags[cid]["frozen"] = True
    else:
        case_flags[cid]["ffpe"] = True

# Step 5: Collect matched case IDs
matched_case_ids = [cid for cid, flags in case_flags.items() if flags["frozen"] and flags["ffpe"]]
print(f"\nMatched cases with both frozen and FFPE samples: ({len(matched_case_ids)})")
for cid in matched_case_ids:
    print(cid)
