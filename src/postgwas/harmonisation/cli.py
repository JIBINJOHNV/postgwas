import argparse
import multiprocessing
import sys
from rich_argparse import RichHelpFormatter 
from postgwas.utils.main import validate_path
import textwrap
from postgwas.clis.common_cli import get_defaultresourse_parser



def get_harmonisation_parser():
    """Returns the ArgumentParser for the PostGWAS harmonisation step, structured with a group."""
    
    # add_help=False is VITAL for inheritance
    parser = argparse.ArgumentParser(
        add_help=False, 
        formatter_class=RichHelpFormatter
    ) 

    # 💥 CRITICAL FIX: DEFINE ARGUMENT GROUP HERE 💥
    harmonisation_group = parser.add_argument_group('HARMONISATION Arguments')
    
    # --- UPDATED --CONFIG HELP MESSAGE ---
    config_help_text =textwrap.dedent(
        '''[bold bright_red]Required[/bold bright_red]: Configuration CSV file with full path. 
        It should contain the following columns:\n
            sumstat_file gwas_outputname chr_col pos_col snp_id_col ea_col oa_col eaf_col beta_or_col se_col\n
            imp_z_col pval_col ncontrol_col ncase_col ncontrol ncase imp_info_col delimiter\n
            infofile infocolumn eaffile eafcolumn liftover chr_pos_col resourse_folder output_folder
            '''
    )

    # All arguments now use the 'harmonisation_group' object
    harmonisation_group.add_argument(
        "--config",
        type=validate_path(must_exist=True, must_be_file=True),
        metavar='', 
        help=config_help_text 
    )

    harmonisation_group.add_argument(
        "--defaults",
        type=validate_path(must_exist=True, must_be_file=True),
        metavar='', 
        help="[bold bright_red]Required[/bold bright_red]: YAML file with full path"
    )
    
    return parser

# --- 2. EXECUTABLE LOGIC (The "Pipeline Step") ---
def run_harmonisation(args):
    """
    Executes the logic for step 1. 
    Accepts 'args' namespace, which must contain: config, defaults, nthreads.
    Returns the dataframe so Step 2 can use it.
    """
    # Move imports here to avoid circular dependency issues later in the pipeline
    from postgwas.harmonisation.io import read_sumstats, read_config
    from postgwas.harmonisation.main import run_harmonisation_pipeline

    # Load USER config
    try:
        cfg_list = read_config(args.config)
        if not cfg_list or len(cfg_list) == 0:
            raise ValueError("read_config returned empty output.")
        user_cfg = cfg_list[0]
    except Exception as e:
        # We raise error here instead of sys.exit so the pipeline can catch it if needed
        raise RuntimeError(f"❌ ERROR: Failed to load user config file '{args.config}'. Reason: {e}")

    # Run pipeline logic
    postgwas_qc_df = run_harmonisation_pipeline(
        sample_column_dict=user_cfg,
        default_cfg=args.defaults,
        nthreads=args.nthreads
    )
    
    return postgwas_qc_df

# --- 3. STANDALONE ENTRY POINT (The "CLI") ---
def main():
    """
    Entry point only used when running: python cli.py
    """
    parser = argparse.ArgumentParser(
        description="Run gwas sumstat harmonisation",
        parents=[get_defaultresourse_parser(), get_harmonisation_parser()],
        formatter_class=RichHelpFormatter
    )

    args = parser.parse_args()

    # -------- FIX: If no arguments → show help instead of running -----
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    try:
        df = run_harmonisation(args)
        print("--- Harmonisation Complete ---")
        print(df)

    except Exception as e:
        print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()