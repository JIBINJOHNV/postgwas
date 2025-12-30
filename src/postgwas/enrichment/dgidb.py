
import requests
import pandas as pd

def query_dgidb(gene_list):
    """
    Queries DGIdb v5 (GraphQL) for drug-gene interactions.
    """
    url = "https://dgidb.org/api/graphql"
    # GraphQL Query
    # We ask for: Gene Name -> Interactions -> Drug Name, Approval, Score, & Sources
    query = """
    query getInteractions($genes: [String!]!) {
      genes(names: $genes) {
        nodes {
          name
          interactions {
            interactionScore
            interactionTypes {
              type
              directionality
            }
            drug {
              name
              conceptId
              approved
              drugApprovalRatings {
                rating
              }
            }
            sources {
              sourceDbName
              pmid
            }
          }
        }
      }
    }
    """
    variables = {"genes": gene_list}
    try:
        print(f"Querying DGIdb for {len(gene_list)} genes...")
        response = requests.post(url, json={'query': query, 'variables': variables})
        response.raise_for_status()
        data = response.json()
        # Parse the nested JSON response
        all_rows = []
        # Check if data exists
        if data.get('data') and data['data'].get('genes'):
            for gene_node in data['data']['genes']['nodes']:
                gene_symbol = gene_node['name']
                for interaction in gene_node['interactions']:
                    # Extract Drug Info
                    drug_name = interaction['drug']['name']
                    is_approved = interaction['drug']['approved']
                    score = interaction['interactionScore']
                    # Extract Interaction Types (e.g., "inhibitor", "antagonist")
                    types = [t['type'] for t in interaction['interactionTypes']]
                    type_str = ", ".join(types) if types else "n/a"
                    # Extract Sources (e.g., "DrugBank", "ChEMBL")
                    sources = [s['sourceDbName'] for s in interaction['sources']]
                    source_str = ", ".join(set(sources)) # Remove duplicates
                    # Extract PMIDs (PubMed IDs)
                    pmids = [str(s['pmid']) for s in interaction['sources'] if s['pmid']]
                    pmid_str = ", ".join(set(pmids))
                    all_rows.append({
                        "Gene": gene_symbol,
                        "Drug": drug_name,
                        "Interaction Type": type_str,
                        "Approved": is_approved,
                        "Score": score,
                        "Sources": source_str,
                        "PMIDs": pmid_str
                    })
        return pd.DataFrame(all_rows)
    except Exception as e:
        print(f"DGIdb Request Failed: {e}")
        return pd.DataFrame()

# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Genes
    my_genes = ["EGFR", "BRAF", "TP53", "TNF", "BRCA1"]
    # 2. Run Query
    df = query_dgidb(my_genes)
    if not df.empty:
        # 3. Filter Results
        # Example: Show only FDA Approved drugs that are "inhibitors"
        approved_inhibitors = df[
            (df['Approved'] == True) & 
            (df['Interaction Type'].str.contains("inhibitor", case=False))
        ]
        print(f"\nTotal Interactions Found: {len(df)}")
        print(f"FDA Approved Inhibitors: {len(approved_inhibitors)}\n")
        # Display Top 10
        display_cols = ['Gene', 'Drug', 'Interaction Type', 'Score', 'Sources']
        print(approved_inhibitors[display_cols].head(10).to_string(index=False))
        # 4. Save to CSV
        # df.to_csv("dgidb_results.csv", index=False)
        # print("\nSaved to dgidb_results.csv")
    else:
        print("No interactions found or API error.")