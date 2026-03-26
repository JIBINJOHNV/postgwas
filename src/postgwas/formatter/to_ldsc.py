import subprocess
from pathlib import Path
import polars as pl
import textwrap
import time


## reason to use effective sample size and sample prevalence 0.5 ; https://groups.google.com/g/ldsc_users/c/yJT-_qSh_44/m/MmKKJYsBAwAJ


## reason to use effective sample size and sample prevalence 0.5

# =========================================================
# LOGGER
# =========================================================
def log_message(log_file: Path, message: str):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    msg = f"[{ts}] {message}"
    print(msg)
    with open(log_file, "a") as lf:
        lf.write(msg + "\n")


# =========================================================
# STEP 1: DROP / FIX N_CAS / N_CON
# =========================================================
def drop_if_mostly_empty(input_file: str, log_file: str, threshold=0.8):

    input_path = Path(input_file)
    log_path = Path(log_file)

    log_message(log_path, f"[STEP] Checking N_CAS / N_CON in: {input_path}")

    df = pl.read_csv(input_path, separator="\t")
    total_rows = df.height

    cols = ["N_CAS", "N_CON"]
    empty_fraction_map = {}
    present_cols = []

    sample_prev = None  # ✅ NEW

    # ---------------------------
    # CHECK EMPTY FRACTION
    # ---------------------------
    for c in cols:
        if c not in df.columns:
            log_message(log_path, f"[SKIP] {c} not present")
            continue

        present_cols.append(c)

        empty_count = df.select(pl.col(c).is_null().sum()).item()
        empty_fraction = empty_count / total_rows

        empty_fraction_map[c] = empty_fraction

        log_message(
            log_path,
            f"[CHECK] {c}: {empty_count}/{total_rows} ({empty_fraction:.2%}) empty"
        )

    # ---------------------------
    # 🧬 SAMPLE PREVALENCE
    # ---------------------------
    if all(c in df.columns for c in ["N_CAS", "N_CON"]):
        try:
            n_cases = df.select(pl.col("N_CAS").median()).item()
            n_controls = df.select(pl.col("N_CON").median()).item()

            if n_cases is not None and n_controls is not None and (n_cases + n_controls) > 0:
                sample_prev = n_cases / (n_cases + n_controls)

                log_message(
                    log_path,
                    f"[SPREV] Sample prevalence = {sample_prev:.4f} "
                    f"(median N_CAS={n_cases}, N_CON={n_controls})"
                )
            else:
                log_message(log_path, "[SPREV] Unable to compute (invalid values)")

        except Exception as e:
            log_message(log_path, f"[SPREV] Error computing prevalence: {e}")

    else:
        log_message(log_path, "[SPREV] Not computed (missing N_CAS or N_CON)")

    # ---------------------------
    # DECISION LOGIC
    # ---------------------------
    cols_to_drop = []

    if len(present_cols) == 2:
        if all(empty_fraction_map[c] >= threshold for c in present_cols):
            log_message(log_path, "[DROP] Both N_CAS & N_CON mostly empty → dropping both")
            cols_to_drop = present_cols
        else:
            log_message(log_path, "[KEEP] Both N_CAS & N_CON retained")

    elif len(present_cols) == 1:
        c = present_cols[0]

        if c == "N_CON":
            log_message(log_path, "[FIX] Only N_CON → renaming to N")
            df = df.rename({"N_CON": "N"})

        elif c == "N_CAS":
            log_message(log_path, "[FIX] Only N_CAS → renaming to N")
            df = df.rename({"N_CAS": "N"})

    else:
        log_message(log_path, "[INFO] No N columns found")

    # ---------------------------
    # APPLY DROP
    # ---------------------------
    if cols_to_drop:
        df = df.drop(cols_to_drop)

    # ---------------------------
    # SAVE
    # ---------------------------
    df.write_csv(input_path, separator="\t")
    log_message(log_path, "[DONE] N column check complete")

    # ✅ RETURN BOTH
    return str(input_path), sample_prev


# =========================================================
# STEP 2: VCF → LDSC
# =========================================================
def vcf_to_ldssc_munge_input(sumstat_vcf: str, output_folder: str, sample_name: str):

    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{sample_name}_ldsc_input.tsv"
    log_file = output_dir / f"{sample_name}_vcf_to_ldsc.log"

    cmd = textwrap.dedent(f"""
        {{
            printf "SNP\\tA2\\tA1\\tZ\\tP\\tN_CAS\\tFRQ\\tINFO\\tN_CON\\n"
            bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
            bcftools query -f '%ID\\t%REF\\t%ALT\\t[%EZ]\\t[%LP]\\t[%NC]\\t[%AF]\\t[%SI]\\t[%NCO]\\n' | \
            sed 's|:|_|g' | \
            awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
                raw_p = exp(-$5 * log(10))
                $5 = sprintf("%.6g", raw_p)
                print
            }}'
        }} > "{output_file}"
    """)

    # write command
    with open(log_file, "w") as lf:
        lf.write("### COMMAND\n\n")
        lf.write(cmd + "\n\n")

    log_message(log_file, "[STEP] Running VCF → LDSC")

    try:
        with open(log_file, "a") as lf:
            subprocess.run(
                cmd,
                shell=True,
                executable="/bin/bash",
                check=True,
                stdout=lf,
                stderr=lf,
            )
    except subprocess.CalledProcessError:
        log_message(log_file, "❌ Conversion failed")
        raise RuntimeError(f"vcf_to_ldsc failed → {log_file}")

    log_message(log_file, "[DONE] VCF → LDSC complete")

    return {
        "ldsc_file": str(output_file),
        "log_file": str(log_file),
    }


# =========================================================
# PIPELINE WRAPPER (UNCHANGED STRUCTURE)
# =========================================================
def vcf_to_ldsc(sumstat_vcf, output_folder, sample_name):

    result = vcf_to_ldssc_munge_input(sumstat_vcf, output_folder, sample_name)

    result["ldsc_file"], sample_prev=drop_if_mostly_empty( 
        input_file=result["ldsc_file"],
        log_file=result["log_file"],
        threshold=0.8
    )

    return {
        "ldsc_file": str(result["ldsc_file"]),
        "log_file": str(result["log_file"]),
        "sample_prev": sample_prev,
    }