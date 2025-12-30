from gprofiler import GProfiler
import pandas as pd
from scipy.stats import hypergeom

# https://biit.cs.ut.ee/gprofiler/page/apis

def run_gprofiler_with_raw_stats(gene_list):
    gp = GProfiler(return_dataframe=True)
    # Run enrichment (returns adjusted p-values in 'p_value' column)
    results = gp.profile(
        organism='hsapiens',
        query=gene_list,
        sources=[ "GO:BP", "GO:MF", "GO:CC",  "KEGG", "REAC", "WP", "TF", "MIRNA", "HPA", "CORUM", "HP" ],
        user_threshold=1.0,           # Allow non-significant
        significance_threshold_method='fdr', 
        all_results=True, # Keep everything
        no_evidences=False         
    )
    if results.empty:
        return results
    # --- CALCULATE RAW P-VALUE ---
    # We use the Hypergeometric Survival Function (1 - CDF)
    # M = Population size (effective_domain_size)
    # n = Number of successes in population (term_size)
    # N = Sample size (query_size)
    # k = Number of successes in sample (intersection_size)
    results['raw_p_value'] = results.apply(
        lambda row: hypergeom.sf(
            row['intersection_size'] - 1, 
            row['effective_domain_size'], 
            row['term_size'], 
            row['query_size']
        ), axis=1
    )
    return results



# --- Main Execution ---
if __name__ == "__main__":
    # 1. Input Genes
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "APOE", "IL6"]
    print(f"Input Genes: {my_genes}")
    # 2. Run Analysis
    print("Running g:Profiler enrichment...")
    df = run_gprofiler_with_raw_stats(my_genes)
    if not df.empty:
        # 3. Clean and Sort Data
        # Select relevant columns
        cols = ['source', 'native', 'name', 'p_value', 'term_size', 'query_size', 'intersection_size', 'intersections']
    
        # 'intersections' is a list of gene names involved in the term
        # Let's clean it up if it exists in the dataframe
        if 'intersections' in df.columns:
            # Flatten list to string
            df['intersections'] = df['intersections'].apply(lambda x: ", ".join(x) if isinstance(x, list) else str(x))

        # Filter columns that exist
        final_cols = [c for c in cols if c in df.columns]
        df_clean = df[final_cols].sort_values(by=['source', 'p_value'])

        # 4. Display Results
        print("\nTop 5 Results:")
        print(df_clean.head().to_string(index=False))

        # OPTIONAL: Save to CSV
        # df_clean.to_csv("gprofiler_results.csv", index=False)
        # print("\nResults saved to gprofiler_results.csv")
    else:
        print("No results returned.")