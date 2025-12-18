

from postgwas.ld_clump.ld_prune import ld_clump_by_regions,ld_clump_standard

            
def run_ld_clump_standard_direct():
    pass 


def run_ld_clump_direct(args, ctx=None):
    """
    Wrapper that dispatches the correct direct-mode LD clumping function.
    """
    if args.ld_mode == "by_regions":
        outputs=ld_clump_by_regions(
            sumstat_vcf=args.vcf,
            output_folder=args.outdir,
            sample_name=args.sample_id,
            population=args.population,
            nthreads=args.nthreads,
            bcftools=args.bcftools
        )
        if ctx is not None:
            ctx["ld_clump"] = outputs
        return outputs

    elif args.ld_mode == "standard":
        return run_ld_clump_standard_direct(
            vcf_path=args.vcf,
            output_dir=args.outdir,
            sample_id=args.sample_id,
            r2_cutoff=args.r2_cutoff,
            window_kb=args.window_kb,
            nthreads=args.nthreads
        )





def run_ld_clump_pipeline(args):
    """
    Wrapper that dispatches the correct PIPELINE mode function.
    """
    if args.ld_mode == "by_regions":
        pass
        #return run_ld_clump_by_regions_pipeline(args)

    elif args.ld_mode == "standard":
        pass
        #return run_ld_clump_standard_pipeline(args)




