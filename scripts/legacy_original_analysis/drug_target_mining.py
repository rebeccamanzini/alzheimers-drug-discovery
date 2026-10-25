# ============================================================
# Drug-target mining via DGIdb GraphQL API
# Input: gene_list_for_drug_mining.csv (from R pipeline)
# Output: drug_gene_interactions.csv
# ============================================================

import requests
import pandas as pd

genes_df = pd.read_csv("gene_list_for_drug_mining.csv")
genes = genes_df['gene'].tolist()
genes_top50 = genes[:50]  # top 50 by adjusted p-value

url = "https://dgidb.org/api/graphql"

query = """
query genes($names: [String!]) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        drug {
          name
        }
        interactionScore
        interactionTypes {
          type
        }
      }
    }
  }
}
"""

variables = {"names": genes_top50}
response = requests.post(url, json={"query": query, "variables": variables}, timeout=30)
data = response.json()

rows = []
for gene_node in data["data"]["genes"]["nodes"]:
    gene_name = gene_node["name"]
    for interaction in gene_node["interactions"]:
        drug_name = interaction["drug"]["name"]
        score = interaction["interactionScore"]
        types = ", ".join([t["type"] for t in interaction["interactionTypes"]]) \
            if interaction["interactionTypes"] else "unspecified"
        rows.append({"gene": gene_name, "drug": drug_name,
                     "interaction_type": types, "score": score})

df = pd.DataFrame(rows)
df = df.sort_values(by="score", ascending=False)
df.to_csv("drug_gene_interactions.csv", index=False)

print(f"Total interactions found: {len(df)}")
print(f"Genes with at least one known drug interaction: {df['gene'].nunique()}")
print(df.head(20))
