import zipfile
import sys
import json
import os

def parse_fastqc_data(zip_path):

    # Parse the FastQC data from a ZIP file.
    # Extracts the 'Per base sequence content' section for further analysis.

    if not os.path.exists(zip_path):
        raise FileNotFoundError(f"File does not exist: {zip_path}")
    
    try:
        with zipfile.ZipFile(zip_path) as zf:
            print(f"Successfully opened ZIP file: {zip_path}")
            data_file = next((f for f in zf.namelist() if f.endswith('fastqc_data.txt')), None)
            if not data_file:
                raise ValueError(f"Invalid FastQC report: Missing fastqc_data.txt in {zip_path}")
            
            print(f"Found fastqc_data.txt file: {data_file}")
            with zf.open(data_file) as f:
                content = [line.decode('utf-8').strip() for line in f]
    except zipfile.BadZipFile:
        raise ValueError(f"Invalid ZIP file: {zip_path}")
    
    records = []
    in_target_section = False
    
    for i, line in enumerate(content):
        # Check for the start of the 'Per base sequence content' section
        if line.startswith('>>Per base sequence content'):
            in_target_section = True
            print(f"Found >>Per base sequence content section at line {i+1}")
            continue
        # Check for the end of the section
        if line.startswith('>>END_MODULE') and in_target_section:
            print(f"End of >>Per base sequence content section at line {i+1}")
            in_target_section = False
            break
        # Process lines within the section
        if in_target_section and line and not line.startswith('#'):
            # Handle comma-separated values
            parts = line.replace(',', '.').split()
            if len(parts) == 5:
                records.append(parts)
                print(f"Recorded data line: {parts}")
    
    if not in_target_section:
        print("Could not find or fully read the >>Per base sequence content section.")
    
    return records

def calculate_trim(records, direction='front'):

    #Calculate the number of bases to trim from the front or end of the reads.
    #Uses the 'Per base sequence content' data to identify abnormal regions.

    total_trim = 0
    iterator = records if direction == 'front' else reversed(records)
    counting = False
    previous_was_normal = False
    
    for parts in iterator:
        if len(parts) != 5:
            continue
            
        base_range = parts[0]
        if '-' in base_range:
            start, end = map(int, base_range.split('-'))
            span = end - start + 1
        else:
            span = 1
        
        try:
            g = float(parts[1])
            a = float(parts[2])
            t = float(parts[3])
            c = float(parts[4])
        except (ValueError, IndexError):
            continue  # Skip lines that cannot be parsed
            
        a_t_diff = abs(a - t)
        c_g_diff = abs(c - g)
        
        # Adjust threshold for abnormal conditions (keep at 1)
        is_abnormal = (a_t_diff > 1) or (c_g_diff > 1)
        
        if is_abnormal:
            counting = True
            total_trim += span
            print(f"Trim count increased by: {span}, Total trim: {total_trim}")
            previous_was_normal = False
        elif counting:
            if previous_was_normal:
                break  # Stop counting after encountering two consecutive normal regions
            else:
                previous_was_normal = True
    
    return total_trim

def main():
    # Main function to process FastQC reports and calculate trimming parameters.
    # Outputs the trimming parameters as a JSON file.
    try:
        r1_zip, r2_zip = snakemake.input
        output_file = snakemake.output[0]
        
        print(f"Processing input files: {r1_zip}, {r2_zip}")
        print(f"Output file: {output_file}")
        
        # Parse FastQC data for R1
        r1_data = parse_fastqc_data(r1_zip)
        f_value = calculate_trim(r1_data, 'front')
        t_value = calculate_trim(r1_data, 'end')
        
        # Parse FastQC data for R2
        r2_data = parse_fastqc_data(r2_zip)
        F_value = calculate_trim(r2_data, 'front')
        T_value = calculate_trim(r2_data, 'end')
        
        # Save trimming parameters to a JSON file
        params = {
            "f_value": f_value,
            "F_value": F_value,
            "t_value": t_value,
            "T_value": T_value
        }
        
        with open(output_file, 'w') as f:
            json.dump(params, f, indent=2)
            print(f"Results written to file: {output_file}")
            
    except Exception as e:
        error_msg = f"Error processing {snakemake.input}: {str(e)}"
        print(error_msg, file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()



