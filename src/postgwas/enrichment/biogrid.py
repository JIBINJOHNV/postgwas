import requests
import pandas as pd
import io

def fetch_biogrid_interactions(gene_list, access_key, tax_id=9606):
    """
    Fetches physical interactions from BioGRID for a list of genes.
    Args:
        gene_list (list): List of gene symbols (e.g., ['TP53', 'BRCA1'])
        access_key (str): Your BioGRID API Key.
        tax_id (int): 9606 = Human, 10090 = Mouse.
    """
    base_url = "https://webservice.thebiogrid.org/interactions/"
    # BioGRID expects genes separated by pipes (|)
    genes_query = "|".join(gene_list)
    params = {
        "accesskey": access_key,
        "geneList": genes_query,
        "taxId": tax_id,
        "searchNames": "true",      # Allow using Symbols instead of IDs
        "includeHeader": "true",    # Needed for pandas to read columns
        "throughputTag": "any",     # 'high', 'low', or 'any'
        "format": "tab2"            # Tab-delimited format is easiest for Pandas
    }
    try:
        print(f"Querying BioGRID for {len(gene_list)} genes...")
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        # BioGRID returns a TSV string
        # We wrap it in StringIO so pandas can read it like a file
        data_io = io.StringIO(response.text)
        df = pd.read_csv(data_io, sep="\t")
        return df
    except Exception as e:
        print(f"BioGRID Request Failed: {e}")
        return pd.DataFrame()

# --- Main Execution ---
if __name__ == "__main__":
    # 1. SETUP: You MUST paste your key here
    # Register at: https://webservice.thebiogrid.org/
    MY_ACCESS_KEY = "86cf1c2c7c8b2972cc5e6b31c45e8584"
    
    my_genes = ["BRCA1", "TP53", "EGFR", "MYC"]

    if MY_ACCESS_KEY == "PASTE_YOUR_ACCESS_KEY_HERE":
        print("[!] Error: You must replace MY_ACCESS_KEY with a valid BioGRID key.")
    else:
        # 2. Fetch Data
        df_biogrid = fetch_biogrid_interactions(my_genes, MY_ACCESS_KEY)

        if not df_biogrid.empty:
            # 3. Filter for PHYSICAL interactions only
            # BioGRID contains 'genetic' interactions (e.g. synthetic lethality), 
            # which you might want to remove if you only care about binding.
            df_phys = df_biogrid[df_biogrid['Experimental System Type'] == 'physical']
            
            print(f"\nTotal Interactions: {len(df_biogrid)}")
            print(f"Physical Interactions: {len(df_phys)}\n")
            
            # 4. Clean up columns for display
            # BioGRID returns A LOT of columns. Let's select the useful ones.
            cols = [
                'Official Symbol Interactor A', 
                'Official Symbol Interactor B', 
                'Experimental System Name', 
                'Publication Source',  # Pubmed ID or Author
                'Throughput'
            ]
            
            # Display Top 10
            print(df_phys[cols].head(10).to_string(index=False))
            
            # OPTIONAL: Save to CSV
            # df_phys.to_csv("biogrid_physical_interactions.csv", index=False)
        else:
            print("No interactions found.")