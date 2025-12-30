
import omnipath as op
import pandas as pd

# https://omnipathdb.org/ ; https://omnipath.readthedocs.io/en/latest/api.html 
# 
def get_interactions_official(gene_list):
    """
    Fetches interactions using the official 'omnipath' library.
    This handles all parameter naming automatically.
    """
    try:
        # 1. Fetch the Human Network
        # This automatically maps organism=9606 and genesymbols=True
        print("Fetching OmniPath data via official client...")
        df = op.interactions.OmniPath.get(genesymbols=True, organism='human')
        # 2. Filter for your genes
        if gene_list:
            # Columns are usually 'source_genesymbol' and 'target_genesymbol'
            # OR 'source' and 'target' depending on the client version
            # Find the correct column names dynamically
            src_col = 'source_genesymbol' if 'source_genesymbol' in df.columns else 'source'
            tgt_col = 'target_genesymbol' if 'target_genesymbol' in df.columns else 'target'
            df_filtered = df[
                (df[src_col].isin(gene_list)) | 
                (df[tgt_col].isin(gene_list))
            ]
            return df_filtered
        return df
    except Exception as e:
        print(f"OmniPath Client Error: {e}")
        return pd.DataFrame()

# --- Usage ---
if __name__ == "__main__":
    my_genes = ["TP53", "EGFR", "MYC", "BRCA1"]
    
    df_result = get_interactions_official(my_genes)
    
    if not df_result.empty:
        print(f"\nSuccess! Found {len(df_result)} interactions.")
        print(df_result.head())
    else:
        print("No interactions found.")
        