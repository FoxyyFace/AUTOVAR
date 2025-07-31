import os
import subprocess
import argparse
from pathlib import Path

def find_available_environments(workflow_dir):
    """Scan and return all available environment configurations"""
    env_configs = []
    workflow_path = Path(workflow_dir)
    
    if not workflow_path.exists():
        print(f"Error: Workflow directory does not exist {workflow_path}")
        return []
    
    # Scan all workflow subdirectories
    for sub_dir in workflow_path.iterdir():
        if not sub_dir.is_dir():
            continue
            
        envs_dir = sub_dir / "envs"
        if not envs_dir.exists():
            continue
            
        # Collect all environment files for this rule
        rule_envs = []
        for yaml_file in envs_dir.glob("*.yaml"):
            rule_envs.append({
                "yaml_path": yaml_file,
                "env_name": yaml_file.stem,
                "target_dir": None  # To be set later
            })
        
        if rule_envs:
            env_configs.append({
                "rule_name": sub_dir.name,
                "envs": rule_envs,
                "rule_path": sub_dir
            })
    
    return env_configs

def install_environment(env_config, envs_base = "conda_env"):
    """Install all environments for a single rule - using Mamba"""
    rule_name = env_config["rule_name"]
    for env_info in env_config["envs"]:
        env_name = env_info["env_name"]
        yaml_file = env_info["yaml_path"]
        
        # Determine target path
        target_dir = Path(envs_base) / rule_name / env_name
        target_dir.parent.mkdir(parents = True, exist_ok = True)
        env_info["target_dir"] = target_dir
        
        print(f"\nCreating environment for {rule_name}/{env_name}...")
        print(f"Configuration file: {yaml_file}")
        print(f"Environment path: {target_dir}")
        
        # Build the installation command - force use of Mamba
        cmd = [
            "mamba", "create", "--yes",
            "--prefix", str(target_dir),
            "--file", str(yaml_file)
        ]
        
        # Execute command
        try:
            print(f"Executing command: {' '.join(cmd)}")
            result = subprocess.run(cmd, check = True, capture_output = True, text = True)
            print(f"Successfully created environment: {target_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Error: Failed to create environment, exit code {e.returncode}")
            print(f"Error details: {e.stderr}")
            return False
        except FileNotFoundError:
            print(f"Error: Mamba command not found, please ensure Mamba is installed")
            print(f"Installation guide: conda install -n base -c conda-forge mamba")
            return False
    
    return True

def main():
    parser = argparse.ArgumentParser(description = "Interactive setup of Snakemake workflow Conda environments - using Mamba")
    parser.add_argument("--workflow", default = "workflow", 
                        help = "Workflow directory path (default: workflow)")
    parser.add_argument("--envs", default = "conda_env", 
                        help = "Conda environments storage path (default: conda_env)")
    parser.add_argument("-y", "--yes", action = "store_true",
                        help = "Automatically install all environments (skip user selection)")
    args = parser.parse_args()
    
    # Find all available environment configurations
    env_configs = find_available_environments(args.workflow)
    
    if not env_configs:
        print(f"No valid environment configurations found in {args.workflow}")
        return
    
    # Display available rules list
    print("\nDetected the following available rule environment configurations:")
    for i, config in enumerate(env_configs, 1):
        env_list = ", ".join([e["env_name"] for e in config["envs"]])
        print(f"{i}. {config['rule_name']} (contains environments: {env_list})")
    
    # Check if skipping interaction
    if args.yes:
        print("\nAutomatically installing all environments...")
        for config in env_configs:
            install_environment(config, args.envs)
        print("\nAll environments installed successfully!")
        return
    
    # User interaction selection
    print("\nPlease choose an option:")
    print("1. Install all rule environments")
    print("2. Select specific rules to install")
    print("3. Exit")
    
    choice = input("Enter option number: ").strip()
    
    if choice == "1":
        # Install all environments
        print("\nStarting installation of all rule environments...")
        for config in env_configs:
            install_environment(config, args.envs)
        print("\nAll environments installed successfully!")
        
    elif choice == "2":
        # Let user select specific rules
        print("\nEnter the rule numbers to install (or enter 0 to cancel/exit)")
        selected = input("Multiple numbers separated by spaces: ").split()
        
        # Check if user entered 0 (cancel)
        if "0" in selected:
            print("Installation canceled. Exiting program.")
            return
        
        # Convert to indices
        selected_indices = [int(idx) - 1 for idx in selected if idx.isdigit()]
        
        if not selected_indices:
            print("No valid rules selected, exiting.")
            return
        
        print("\nStarting installation of selected rule environments...")
        for idx in selected_indices:
            if 0 <= idx < len(env_configs):
                install_environment(env_configs[idx], args.envs)
            else:
                print(f"Warning: Skipping invalid index {idx + 1}")
        print("\nSelected environments installed successfully!")
        
    else:
        print("Exiting program.")

if __name__ == "__main__":
    # Check if Mamba is available
    try:
        subprocess.run(["mamba", "--version"], 
                       stdout = subprocess.DEVNULL, 
                       stderr = subprocess.DEVNULL)
    except FileNotFoundError:
        print("Error: Mamba command not found, please install Mamba first")
        print("Installation guide: conda install -n base -c conda-forge mamba")
        exit(1)
    
    main()