from pathlib import Path
import pandas as pd
import subprocess
import sys
import os
from postgwas.enrichment.utils import load_gene_list
from postgwas.enrichment.cli import get_geneset_parser



def main():
    # ---------------------------------------------------------
    # 1. ENVIRONMENT CHECK & DISPATCHER
    # ---------------------------------------------------------
    # Get current environment name
    current_env = os.environ.get("CONDA_DEFAULT_ENV", "unknown")
    target_env = "enricher"

    # If we are NOT in the enricher environment, we must switch temporarily
    if current_env != target_env:
        print(f"🔄 Current Environment: '{current_env}'")
        print(f"🚀 Switching to '{target_env}' environment to run analysis...")
        
        # 1. Construct the command to run THIS same script inside the 'enricher' env
        #    We pass all arguments (sys.argv[1:]) exactly as received
        cmd = [
            "micromamba", "run", "-n", target_env,
            "python", sys.argv[0]
        ] + sys.argv[1:]

        try:
            # 2. Run the analysis in the child process and wait for it to finish
            subprocess.check_call(cmd)
            
            # 3. "Reactivate" / Return to PostGWAS
            print("-" * 60)
            print(f"✅ Enricher analysis complete.")
            print(f"🔙 Reactivating '{current_env}' environment...")
            print(f"🔹 You are now back in '{current_env}' and can run other postgwas modules.")
            return

        except subprocess.CalledProcessError as e:
            print(f"❌ Analysis failed with error code {e.returncode}")
            sys.exit(e.returncode)
        except KeyboardInterrupt:
            print("\n⚠️ Interrupted by user.")
            sys.exit(1)

    # =========================================================
    # 2. WORKER LOGIC (Runs ONLY inside 'enricher' env)
    # =========================================================
    # We use LAZY IMPORTS here. If we imported these at the top, 
    # the script would crash immediately in the 'postgwas' env.
    from postgwas.enrichment.workflows import run_multisource_enrichment_pipeline
    from postgwas.enrichment.utils import load_gene_list

    # Actually PARSE the arguments
    parser = get_geneset_parser()
    args = parser.parse_args()

    print(f"⚙️  Processing inside '{current_env}' environment...")

    # -------------------------
    # Load genes
    # -------------------------
    try:
        gene_list = load_gene_list(args.genes)
        print(f"🧬 Loaded {len(gene_list)} genes from input.")
    except Exception as e:
        print(f"❌ Error loading genes: {e}")
        sys.exit(1)

    # -------------------------
    # Run FULL pipeline
    # -------------------------
    run_multisource_enrichment_pipeline(
        gene_list=gene_list,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        biogrid_access_key=args.biogrid_key,
        david_email=args.david_email,
        dsigdb_gmt=args.dsigdb_gmt,
        score_threshold=args.string_score,
        fdr_thr=args.fdr_thr,
    )

if __name__ == "__main__":
    main()