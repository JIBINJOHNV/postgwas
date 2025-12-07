import os
import polars as pl

def validate_gwas_config(config: dict, df: pl.DataFrame) -> bool:
    """
    Validate a GWAS configuration dictionary and check if required columns exist in the provided Polars DataFrame.
    
    Rules:
    - One of each OR pair must be non-NA/non-empty in config and present in DataFrame.
    - Both EA and OA must be non-NA/non-empty in config and present in DataFrame.
    - Single required keys (pval_col, gwas_outputname, sumstat_file) must be non-NA/non-empty in config.
    - Only pval_col is checked in DataFrame; gwas_outputname and sumstat_file are not columns.
    - Sumstat file must exist on disk.
    
    Args:
        config (dict): GWAS configuration dictionary.
        df (pl.DataFrame): Polars DataFrame to check for columns.
    
    Returns:
        bool: True if valid.
    
    Raises:
        ValueError: If any validation rule fails or columns are missing.
        FileNotFoundError: If sumstat_file doesn't exist.
    """
    def is_valid(value):
        return value not in [None, "", "NA", "na", "NaN"]
    
    # Define validation rules
    or_pairs = [
        ("chr_col", "chr_pos_col"),
        ("pos_col", "chr_pos_col"),
        ("beta_or_col", "imp_z_col"),
        ("eaf_col", "eaffile"),
        ("imp_info_col", "infofile"),
        ("ncontrol_col", "ncontrol"),
        ("ncase_col", "ncase"),
    ]
    and_pairs = [("ea_col", "oa_col")]
    required = ["pval_col", "gwas_outputname", "sumstat_file"]
    # Check OR pairs in config and DataFrame
    df_columns = set(df.columns)
    for a, b in or_pairs:
        a_val, b_val = config.get(a), config.get(b)
        if not (is_valid(a_val) or is_valid(b_val)):
            raise ValueError(f"Either '{a}' or '{b}' must be valid in config.")
        if is_valid(a_val) and a_val not in df_columns and is_valid(b_val) and b_val not in df_columns:
            raise ValueError(f"Either '{a_val}' or '{b_val}' must be present in DataFrame.")
    # Check AND pairs in config and DataFrame
    for a, b in and_pairs:
        a_val, b_val = config.get(a), config.get(b)
        if not (is_valid(a_val) and is_valid(b_val)):
            raise ValueError(f"Both '{a}' and '{b}' must be valid in config.")
        if a_val not in df_columns or b_val not in df_columns:
            raise ValueError(f"Both '{a_val}' and '{b_val}' must be present in DataFrame.")
    # Check required singles in config and pval_col in DataFrame
    for key in required:
        value = config.get(key)
        if not is_valid(value):
            raise ValueError(f"Missing or invalid '{key}' in config.")
        if key == "pval_col" and value not in df_columns:
            raise ValueError(f"Column '{value}' for '{key}' missing in DataFrame.")
    # Check file existence
    if not os.path.exists(config["sumstat_file"]):
        raise FileNotFoundError(f"File not found: {config['sumstat_file']}")
    #print("GWAS configuration and DataFrame columns validated.")
    return True