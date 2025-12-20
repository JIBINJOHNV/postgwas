def print_qc_summary(final_qc_df, columns=None, warn_threshold=0.05):
    """
    Print GWAS QC summary in a human-readable, biologically meaningful format.
    Warnings are shown ONLY when thresholds are exceeded.
    """

    PREFIX = "\t\t"

    if isinstance(columns, str):
        columns = [columns]

    if columns is None:
        columns = ["raw_variant_count"]

    # Convert long → index-based
    final_qc_df = final_qc_df.set_index("index")

    QC_LABELS = {
        "num_records": "Total variant records",
        "num_snps": "Total SNPs",
        "ts_tv_ratio": "Transition / Transversion ratio",
        "variants_with_missing_extrnal_af": "Variants missing external allele frequency",
        "variants_with_EAF_diff_gt_cutoff": "Variants with large EAF discrepancy",
    }

    available_cols = [c for c in columns if c in final_qc_df.columns]
    if not available_cols:
        print(f"{PREFIX}❌ No valid QC columns found.")
        return

    for col in available_cols:
        print(f"\n{PREFIX}📊 GWAS QC Summary — {col}")
        print(f"{PREFIX}" + "-" * (26 + len(col)))

        # Use num_records as denominator
        if "num_records" not in final_qc_df.index:
            print(f"{PREFIX}❌ Missing num_records — cannot compute QC ratios.")
            continue

        total_variants = final_qc_df.loc["num_records", col]

        def frac(val):
            return (val / total_variants) if total_variants and val is not None else 0.0

        # --------------------------------------------------
        # Core QC metrics
        # --------------------------------------------------
        for key, label in QC_LABELS.items():
            if key not in final_qc_df.index:
                continue

            val = final_qc_df.loc[key, col]

            if key == "ts_tv_ratio":
                print(f"{PREFIX}🔬 {label:45s}: {val:.2f}")
            else:
                print(f"{PREFIX}🧬 {label:45s}: {int(val):,}")

        # --------------------------------------------------
        # QC Warnings (only meaningful ones)
        # --------------------------------------------------
        warning_lines = []

        # Missing external AF
        if "variants_with_missing_extrnal_af" in final_qc_df.index:
            ratio = frac(final_qc_df.loc["variants_with_missing_extrnal_af", col])
            if ratio > warn_threshold:
                warning_lines.append(
                    f"{PREFIX}❗ {QC_LABELS['variants_with_missing_extrnal_af']}: "
                    f"{ratio*100:.2f}% "
                    f"({int(final_qc_df.loc['variants_with_missing_extrnal_af', col]):,} variants)"
                )

        # Large EAF discordance (allele flip warning)
        if "variants_with_EAF_diff_gt_cutoff" in final_qc_df.index:
            ratio = frac(final_qc_df.loc["variants_with_EAF_diff_gt_cutoff", col])
            if ratio > warn_threshold:
                warning_lines.append(
                    f"{PREFIX}❗ {QC_LABELS['variants_with_EAF_diff_gt_cutoff']}: "
                    f"{ratio*100:.2f}% "
                    f"({int(final_qc_df.loc['variants_with_EAF_diff_gt_cutoff', col]):,} variants)\n"
                    f"{PREFIX}   → Possible allele flips, strand mismatches, or population mismatch"
                )

        # --------------------------------------------------
        # Print warnings only if present
        # --------------------------------------------------
        if warning_lines:
            print(f"\n{PREFIX}⚠️ QC Warnings")
            print(f"{PREFIX}" + "-" * 13)
            for line in warning_lines:
                print(line)
    
