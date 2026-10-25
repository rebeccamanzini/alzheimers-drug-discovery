
import datetime
import pandas as pd
import requests

genes_df = pd.read_csv("results/gene_list_for_drug_mining.csv")
top = genes_df.head(50)
genes_top50 = top["gene"].tolist()

url = "https://dgidb.org/api/graphql"
query = """
query genes($names: [String!]) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        drug { name }
        interactionScore
        interactionTypes { type }
      }
    }
  }
}
"""
response = requests.post(url, json={"query": query, "variables": {"names": genes_top50}}, timeout=30)
response.raise_for_status()
data = response.json()
if "errors" in data:
    raise RuntimeError(data["errors"])

nodes = data["data"]["genes"]["nodes"]
rows = []
for node in nodes:
    for inter in node["interactions"]:
        types = ", ".join(t["type"] for t in inter["interactionTypes"]) if inter["interactionTypes"] else "unspecified"
        rows.append({"gene": node["name"], "drug": inter["drug"]["name"],
                     "interaction_type": types, "score": inter["interactionScore"]})

df = pd.DataFrame(rows)
df = df.merge(top[["gene", "logFC", "adj.P.Val"]], on="gene", how="left")
df["direction_in_AD"] = df["logFC"].apply(lambda x: "higher in AD" if x > 0 else "lower in AD")
df = df.sort_values(by="score", ascending=False)
df.to_csv("results/drug_gene_interactions.csv", index=False)

not_found = sorted(set(genes_top50) - {n["name"] for n in nodes})
print("Queried DGIdb on", datetime.date.today())
print("Total interactions found:", len(df))
print("Genes with at least one known drug interaction:", df["gene"].nunique())
print("Genes not found in DGIdb:", not_found)
print(df.head(20).to_string())

