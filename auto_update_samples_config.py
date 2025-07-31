import os
import yaml

def get_samples_structure(samples_root="data/samples"):
    samples = {}
    species_list = []
    for species in sorted(os.listdir(samples_root)):
        species_path = os.path.join(samples_root, species)
        if not os.path.isdir(species_path):
            print(f"Skipping {species} as it is not a directory.")
            continue
        species_list.append(species)
        samples[species] = {}
        for sample in sorted(os.listdir(species_path)):
            sample_path = os.path.join(species_path, sample)
            if not os.path.isdir(sample_path):
                print(f"Skipping {sample} as it is not a directory.")
                continue
            fq_files = [
                os.path.join(sample_path, f)
                for f in os.listdir(sample_path)
                if f.endswith(".fq.gz") or f.endswith(".fastq.gz")
            ]
            fq_files.sort()
            if len(fq_files) == 0:
                print(f"Skipping {sample} as it does not have any FASTQ files.")
                continue
            # 判断SE还是PE
            sample_type = "PE" if len(fq_files) == 2 else "SE"
            samples[species][sample] = {
                "type": "RNA",  # 可根据需要修改
                "platform": "Illumina",  # 可根据需要修改
                "fq_paths": fq_files,
                "source": "unknown"  # 可根据需要修改
            }
            print(f"Added sample {sample} of type {sample_type} with files {fq_files}")
    return species_list, samples

def main():
    species_list, samples = get_samples_structure("data/samples")
    config = {
        "species_list": species_list,
        "samples": samples
    }
    with open("configure/samples_config.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False, allow_unicode=True)
    print("samples_config.yaml 已自动更新！")

if __name__ == "__main__":
    main()




