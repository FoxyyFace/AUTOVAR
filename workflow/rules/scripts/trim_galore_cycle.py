import subprocess
import shutil
from pathlib import Path
import uuid
import re

MAX_CYCLES = 15  # Maximum number of cycles

def has_adapter_contamination_count(report_path):
    # Check for adapter contamination (based on absolute count)
    with open(report_path) as f:
        for line in f:
            if "Reads with adapters" in line:
                count_str = line.split(":")[1].strip().split()[0].replace(",", "")
                return int(count_str) > 0  # Return True if any adapters are present
    return False

def has_adapter_contamination_percent(report_path):
    # Check for adapter contamination (based on percentage)
    with open(report_path) as f:
        for line in f:
            if "Reads with adapters" in line:
                match = re.search(r"\((\d+\.\d+)%\)", line)
                if match:
                    return float(match.group(1)) > 0.0  # Return True if percentage > 0
    return False

def check_contamination(report_path, strict_mode):
    # Dynamically select the detection strategy
    if strict_mode:
        return has_adapter_contamination_count(report_path)
    else:
        return has_adapter_contamination_percent(report_path)

def run_trim_galore(input_fq1, input_fq2, output_dir, cycle_num):
    # Execute a single Trim Galore processing cycle
    cmd = [
        "trim_galore",
        "--paired",
        "-q", str(snakemake.params.quality),
        "--max_n", str(snakemake.params.max_n),
        "-j", str(snakemake.params.threads),
        "-o", str(output_dir),
        str(input_fq1),
        str(input_fq2),
    ]
    
    result = subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    
    with open(snakemake.log[0], "a") as log_file:
        log_file.write(f"\n=== Cycle {cycle_num} Log ===\n")
        log_file.write(result.stdout)
    
    return {
        "fq1": next(output_dir.glob("*_1_val_1.fq.gz")),
        "fq2": next(output_dir.glob("*_2_val_2.fq.gz")),
        "report1": next(output_dir.glob("*_1.fq.gz_trimming_report.txt")),
        "report2": next(output_dir.glob("*_2.fq.gz_trimming_report.txt")),
    }

def main():
    strict_mode = snakemake.params.strict_qc
    base_temp_dir = Path(snakemake.output.output_dir) / f"tmp_{uuid.uuid4().hex[:6]}"
    base_temp_dir.mkdir(parents=True, exist_ok=True)
    
    current_files = {
        "fq1": Path(snakemake.input.fq1),
        "fq2": Path(snakemake.input.fq2)
    }
    
    final_output = {
        "fq1": Path(snakemake.output.fq1),
        "fq2": Path(snakemake.output.fq2),
        "report1": Path(snakemake.output[0]) / f"{Path(snakemake.output[0]).name}_1.trimming_report.txt",
        "report2": Path(snakemake.output[0]) / f"{Path(snakemake.output[0]).name}_2.trimming_report.txt",
    }

    try:
        for cycle in range(1, MAX_CYCLES + 1):
            cycle_temp_dir = base_temp_dir / f"cycle_{cycle}"
            cycle_temp_dir.mkdir(exist_ok=True)
            
            processed = run_trim_galore(
                current_files["fq1"], 
                current_files["fq2"], 
                cycle_temp_dir, 
                cycle
            )

            # Check contamination status
            r1_contaminated = check_contamination(processed["report1"], strict_mode)
            r2_contaminated = check_contamination(processed["report2"], strict_mode)
            
            # Record the current best results
            current_files = processed
            
            # Early termination condition
            if not r1_contaminated and not r2_contaminated:
                print(f"Cycle {cycle}: Adapters fully removed")
                break
            
            # Final cycle processing
            if cycle == MAX_CYCLES:
                print(f"Reached maximum cycles ({MAX_CYCLES}), using best available results")
                
        # Move final results (whether early termination or max cycles reached)
        for key in ["fq1", "fq2", "report1", "report2"]:
            shutil.move(current_files[key], final_output[key])
            
    finally:
        shutil.rmtree(base_temp_dir, ignore_errors=True)

if __name__ == "__main__":
    main()