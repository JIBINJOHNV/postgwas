import subprocess
import argparse
import sys
import os
import time
import yaml  # You may need to run: pip install pyyaml
from pathlib import Path

# ---------------------------------------------------------
# 1. ENHANCED DOCKER & RESOURCE CHECKS
# ---------------------------------------------------------

def check_docker_environment():
    """Checks if Docker is installed, and if it's running."""
    print("🔍 Checking Docker environment...")
    try:
        subprocess.run(["docker", "--version"], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ ERROR: Docker is not installed.")
        sys.exit(1)

    try:
        subprocess.run(["docker", "info"], check=True, capture_output=True)
        print("🐳 Docker is running.")
    except subprocess.CalledProcessError:
        print("⚠️ Docker is not running. Attempting to start...")
        if sys.platform == "darwin": subprocess.run(["open", "-a", "Docker"])
        elif sys.platform == "linux": subprocess.run(["sudo", "systemctl", "start", "docker"])
        
        for i in range(1, 7):
            print(f"⏳ Waiting for Docker... ({i}/6)")
            time.sleep(5)
            try:
                subprocess.run(["docker", "info"], check=True, capture_output=True)
                print("✨ Docker started!")
                return
            except subprocess.CalledProcessError: continue
        print("❌ ERROR: Could not start Docker.")
        sys.exit(1)

def update_yaml_resources(yaml_path, resource_folder):
    """Updates the YAML file with the user-provided resource path."""
    print(f"📝 Updating {os.path.basename(yaml_path)} with current resource paths...")
    
    # 1. Validation: Must end with gwas2vcf
    if not resource_folder.rstrip('/').endswith('gwas2vcf'):
        print(f"❌ ERROR: Resource folder must end with 'gwas2vcf'. Provided: {resource_folder}")
        sys.exit(1)

    # 2. Validation: Folder must not be empty
    if not os.listdir(resource_folder):
        print(f"❌ ERROR: Resource folder is empty: {resource_folder}")
        sys.exit(1)

    # 3. Read, Update, and Write YAML
    try:
        with open(yaml_path, 'r') as f:
            config = yaml.safe_load(f)

        # Standardizing trailing slash
        res_base = os.path.join(resource_folder, '') 
        
        config['resource_folder'] = res_base
        config['grch37_file'] = os.path.join(res_base, "GRCh37_38_check_files/GRCh37_check_file.tsv")
        config['grch38_file'] = os.path.join(res_base, "GRCh37_38_check_files/GRCh38_check_file.tsv")

        with open(yaml_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        print("✅ YAML updated successfully.")
    except Exception as e:
        print(f"❌ ERROR updating YAML: {e}")
        sys.exit(1)

# ---------------------------------------------------------
# 2. UTILITIES
# ---------------------------------------------------------

def validate_paths(paths_to_check, create_if_missing=False):
    for label, path in paths_to_check.items():
        p = Path(path)
        if not p.exists():
            if create_if_missing:
                p.mkdir(parents=True, exist_ok=True)
            else:
                print(f"❌ ERROR: {label} missing at {path}")
                sys.exit(1)
    print("✅ Local paths validated.")

def get_common_base(paths):
    valid_paths = [os.path.abspath(p) for p in paths if p and p != "NA"]
    valid_paths.append(os.getcwd())
    return os.path.commonpath(valid_paths)

def run_step(step_name, command):
    print(f"\n--- 🚀 [STEP: {step_name}] ---")
    try:
        subprocess.run(command, check=True)
        print(f"✅ {step_name} successful!")
    except subprocess.CalledProcessError as e:
        print(f"❌ FAILED: {step_name} (Code {e.returncode})")
        sys.exit(e.returncode)

# ---------------------------------------------------------
# 3. MAIN RUNNER
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="PostGWAS Pipeline Runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input-dir", help="REQUIRED: Input GWAS folder")
    parser.add_argument("-r", "--resource-folder", help="REQUIRED: Resource folder ending in gwas2vcf")
    parser.add_argument("-o", "--outdir", default=os.getcwd(), help="Output folder")
    parser.add_argument("-d", "--yaml-defaults", help="Path to harmonisation.yaml")
    parser.add_argument("--threads", default="10")
    parser.add_argument("--memory", default="50G")
    parser.add_argument("--image", default="jibinjv/postgwas:1.3")

    args = parser.parse_args()

    if not args.input_dir or not args.resource_folder:
        print("❌ ERROR: --input-dir and --resource-folder are required.")
        sys.exit(1)

    check_docker_environment()

    if args.yaml_defaults is None:
        args.yaml_defaults = os.path.join(args.input_dir, "harmonisation.yaml")

    # PRE-STEP: Update YAML paths before validation
    validate_paths({
        "Input Directory": args.input_dir,
        "Resource Folder": args.resource_folder,
        "YAML Config": args.yaml_defaults
    })
    
    update_yaml_resources(args.yaml_defaults, args.resource_folder)

    validate_paths({"Output Directory": args.outdir}, create_if_missing=True)

    config_path = os.path.join(args.outdir, "harmonisatio_example_input_file.csv")
    base_vol = get_common_base([args.input_dir, args.resource_folder, args.outdir, args.yaml_defaults])

    # STEP 1: Column Mapping
    run_step("Column Mapping", [
        "docker", "run", "--rm", "--platform", "linux/amd64", 
        "-v", f"{base_vol}:{base_vol}", args.image, 
        "python", "/opt/postgwas/src/postgwas/scripts/create_sumstat_map_pl.py", 
        "--input", args.input_dir, "--resource-folder", args.resource_folder, 
        "--output-path", config_path, "--harmonisation-output-path", args.outdir
    ])

    # STEP 2: Harmonisation
    run_step("Harmonisation", [
        "docker", "run", "--rm", "--platform", "linux/amd64", 
        "-v", f"{base_vol}:{base_vol}", args.image, 
        "postgwas", "harmonisation", "--nthreads", args.threads, 
        "--max-mem", args.memory, "--config", config_path, "--defaults", args.yaml_defaults
    ])

    print(f"\n🎉 COMPLETED. Results: {args.outdir}")

if __name__ == "__main__":
    main()