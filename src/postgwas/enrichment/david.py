import zeep
import pandas as pd
import time

def run_david_enrichment(gene_list, email, id_type='OFFICIAL_GENE_SYMBOL'):
    """
    Runs enrichment analysis using DAVID Web Service (SOAP).
    
    Args:
        gene_list (list): List of gene symbols.
        email (str): YOUR REGISTERED EMAIL (Must register at https://david.ncifcrf.gov/webservice/register.html).
        id_type (str): Identifier type (e.g., 'OFFICIAL_GENE_SYMBOL', 'ENTREZ_GENE_ID').
    """
    # 1. Initialize Client
    url = 'https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl'
    print(f"Connecting to DAVID Web Service...")
    
    try:
        client = zeep.Client(url)
        
        # 2. Authenticate
        # The API requires a registered email address
        result = client.service.authenticate(email)
        if result != 'true':
            print(f"Authentication failed for {email}. Please register at https://david.ncifcrf.gov/webservice/register.html")
            return pd.DataFrame()
        print("Authentication successful.")

        # 3. Add Gene List
        # 1st arg: list of genes (comma-separated string)
        # 2nd arg: identifier type
        # 3rd arg: list name
        # 4th arg: temporary list type (0)
        genes_str = ",".join(gene_list)
        list_name = "Python_List_Upload"
        
        val = client.service.addList(genes_str, id_type, list_name, 0)
        if val < 1:
            print("Failed to map any genes. Check your symbols or ID type.")
            return pd.DataFrame()
        print(f"Uploaded {val} mapped genes.")

        # 4. Get Chart Report (Enrichment)
        # Threshold: EASE score (P-value) usually 0.1
        # count: Min genes per term (usually 2)
        thd = 0.1
        ct = 2
        
        print("Fetching Chart Report (this may take a moment)...")
        # categories can be specified, strictly optional (uses default if skipped)
        # Common categories: 'GOTERM_BP_DIRECT', 'KEGG_PATHWAY', 'UP_KEYWORDS'
        
        # We fetch the default categories report
        report_records = client.service.getChartReport(thd, ct)
        
        if not report_records:
            print("No significant enrichment found.")
            return pd.DataFrame()

        # 5. Parse Results
        rows = []
        for row in report_records:
            rows.append({
                "Category": row['categoryName'],
                "Term": row['termName'],
                "Count": row['listHits'],
                "Percent": row['percent'],
                "P-Value": row['ease'],        # Raw P-Value (EASE Score)
                "Bonferroni": row['bonferroni'],
                "Benjamini": row['benjamini'],
                "FDR": row['fdr'],
                "Genes": row['geneIds']        # Returns IDs, sometimes Symbols depending on config
            })
            
        return pd.DataFrame(rows)

    except Exception as e:
        print(f"Error communicating with DAVID: {e}")
        return pd.DataFrame()

# --- Main Execution ---
if __name__ == "__main__":
    # 1. Your Gene List
    my_genes = ["TP53", "BRCA1", "EGFR", "TNF", "APOE", "IL6"]
    
    # 2. YOUR REGISTERED EMAIL
    # You MUST replace this with your actual registered email
    my_email = "YOUR_REGISTERED_EMAIL@example.com" 
    
    # 3. Run Analysis
    df = run_david_enrichment(my_genes, my_email)
    
    if not df.empty:
        # Sort by P-Value
        df = df.sort_values(by="P-Value")
        
        print("\nTop 5 Enriched Terms:")
        print(df[['Category', 'Term', 'P-Value', 'Benjamini']].head(5).to_string(index=False))
        
        # df.to_csv("david_results.csv", index=False)