import requests
import pandas as pd
import io
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def get_string_ids(gene_list, species=9606):
    """
    Maps gene symbols to STRING protein IDs.
    Species 9606 is Human.
    """
    url = "https://string-db.org/api/json/get_string_ids"
    params = {
        "identifiers": "\r".join(gene_list),
        "species": species,
        "limit": 1,  # Map to the best single hit per gene
        "echo_query": 1
    }
    try:
        response = requests.post(url, data=params)
        response.raise_for_status()
        data = response.json()
        # STRING returns a list of dicts. We want the 'stringId'
        mapped_ids = []
        for item in data:
            mapped_ids.append(item['stringId'])
        return list(set(mapped_ids)) # Remove duplicates
    except Exception as e:
        print(f"Mapping Error: {e}")
        return []

def run_string_enrichment(string_ids, species=9606):
    """
    Performs functional enrichment on the mapped STRING IDs.
    """
    url = "https://string-db.org/api/json/enrichment"
    params = {
        "identifiers": "\r".join(string_ids),
        "species": species
    }
    try:
        response = requests.post(url, data=params)
        response.raise_for_status()
        data = response.json()
        return pd.DataFrame(data)
    except Exception as e:
        print(f"Enrichment Error: {e}")
        return pd.DataFrame()

def get_network_image(string_ids, species=9606, filename="string_network.png"):
    """
    Downloads the static network image for the provided genes.
    """
    url = "https://string-db.org/api/image/network"
    params = {
        "identifiers": "\r".join(string_ids),
        "species": species,
        "add_white_nodes": 0,      # Don't add extra proteins
        "network_flavor": "evidence", # 'evidence' or 'confidence' lines
        "hide_disconnected_nodes": 1  # Clean up the plot
    }
    try:
        print("Downloading network image...")
        response = requests.post(url, data=params)
        response.raise_for_status()
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f"Network image saved to: {filename}")
    except Exception as e:
        print(f"Image Download Error: {e}")


def get_string_ppi_network(gene_list, species=9606, score_threshold=400):
    """
    Fetches protein-protein interactions from STRING-DB.
    Args:
        gene_list (list): List of gene symbols.
        species (int): 9606 for Human, 10090 for Mouse.
        score_threshold (int): Minimum confidence score (0-1000). 
                               400 = Medium, 700 = High, 900 = Highest.
    """
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list),  # Join with newline for API
        "species": species,
        "required_score": score_threshold
    }
    try:
        print(f"Fetching PPI network for {len(gene_list)} genes...")
        response = requests.post(url, data=params)
        response.raise_for_status()
        data = response.json()
        if not data:
            print("No interactions found.")
            return pd.DataFrame()
        # Convert to DataFrame
        # STRING returns: 'preferredName_A', 'preferredName_B', 'score', etc.
        df = pd.DataFrame(data)
        # Select relevant columns: Protein A, Protein B, Confidence Score
        # preferredName_A = Source, preferredName_B = Target
        df_clean = df[['preferredName_A', 'preferredName_B', 'score', 'nscore', 'fscore', 'pscore', 'ascore', 'escore', 'dscore', 'tscore']]
        # Rename for clarity
        df_clean.columns = ['Source', 'Target', 'Total_Score', 'Neighborhood', 'Fusion', 'Cooccurrence', 'Coexpression', 'Experimental', 'Database', 'TextMining']
        return df_clean
    except Exception as e:
        print(f"Error fetching network: {e}")
        return pd.DataFrame()

def analyze_hubs(df_interactions):
    """
    Uses NetworkX to find the most connected proteins (Hubs).
    """
    if df_interactions.empty:
        return
    # 1. Build Graph
    G = nx.from_pandas_edgelist(df_interactions, 'Source', 'Target', ['Total_Score'])
    # 2. Calculate Degree Centrality (Number of connections)
    degrees = dict(G.degree())
    # Sort by connections
    sorted_degrees = sorted(degrees.items(), key=lambda item: item[1], reverse=True)
    print("\n=== TOP 5 HUB GENES (Most Connected) ===")
    for gene, degree in sorted_degrees[:5]:
        print(f"{gene}: {degree} interactions")
    return G


# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Genes
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "AKT1", "MYC"]
    print(f"Input: {my_genes}")

    # 2. Map to STRING IDs
    mapped_ids = get_string_ids(my_genes)
    
    if mapped_ids:
        print(f"Successfully mapped {len(mapped_ids)} proteins.")
        
        # 3. Run Enrichment
        df = run_string_enrichment(mapped_ids)
        
        if not df.empty:
            # Clean and Sort Results
            # Columns usually: category, term, number_of_genes, fdr, description
            df_clean = df[['category', 'term', 'description', 'fdr', 'number_of_genes']]
            df_clean = df_clean.sort_values(by='fdr')
            
            print("\nTop 5 Enriched Terms:")
            print(df_clean.head().to_string(index=False))
            
            # Save CSV
            # df_clean.to_csv("string_enrichment_results.csv", index=False)
        
        # 4. Get Network Image (The cool part of STRING)
        get_network_image(mapped_ids)
        
    else:
        print("No genes could be mapped to STRING IDs.")

    # 2. Get Interactions
    # Score 400 is "Medium Confidence". Increase to 700 for "High Confidence".
    df_ppi = get_string_ppi_network(my_genes, score_threshold=400)
    
    if not df_ppi.empty:
        print(f"\nFound {len(df_ppi)} interactions.")
        
        # Show data
        print(df_ppi[['Source', 'Target', 'Total_Score']].head().to_string(index=False))
        
        # 3. Find Hubs
        G = analyze_hubs(df_ppi)
        
        # 4. Save to CSV for Cytoscape (Optional)
        # df_ppi.to_csv("string_interactions.csv", index=False)
        # print("\nSaved interactions to 'string_interactions.csv' (Load this into Cytoscape!)")
        
        # 5. Simple Visualization (Optional)
        plt.figure(figsize=(8, 8))
        pos = nx.spring_layout(G, seed=42) # Consistent layout
        nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=1500, font_weight='bold')
        plt.title("Protein-Protein Interaction Network")
        plt.show()
    else:
        print("No interactions returned.")