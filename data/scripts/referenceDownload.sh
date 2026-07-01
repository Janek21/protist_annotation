#!/bin/bash

#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --job-name=taxonDownloader

## qos determining
#SBATCH --qos=test

#SBATCH --time=240

## Mem+cpu
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


#record start
echo ">STARTING at $(date)"

sp_out="species"
taxlist="true_taxons.txt"

#download annoations+available assembly
annocli download --taxids-file "$taxlist" --add-asm --fix-alias --output "$sp_out" --ref-only
echo "annocli_done"
#remove aliasNaming
for dir in "$sp_out"/*/GC*/; do
    if [ -d "$dir" ]; then
        #Iterate through EVERY aliasMatch file found in this directory
        #Using .gff*.gz captures both .gff.gz and .gff3.gz safely
        for alias_file in "$dir"*aliasMatch.gff*.gz; do

            # Check if the file exists (handles the case where the glob finds nothing)
            if [ -f "$alias_file" ]; then

                # Derive the target filename (removing '.aliasMatch')
                target_file="${alias_file//.aliasMatch/}"

                # Remove the original "regular" GFF/GFF3 file
                rm -f "$target_file"

                # Rename the aliasMatch file to the standard name
                mv "$alias_file" "$target_file"

                echo "Converted: $(basename "$target_file")"
            fi
        done

        # Clean up the TSV files after the renames are done
        rm -f "$dir"*.tsv
    fi
done

mv species/Phaeodactylum_tricornutum_CCAP_1055/1_556484/GCF* species/Phaeodactylum_tricornutum_CCAP_1055/
rm -rf species/Phaeodactylum_tricornutum_CCAP_1055/1_556484
mv species/Phaeodactylum_tricornutum_CCAP_1055/ species/Phaeodactylum_tricornutum_CCAP_2850

echo "python time"
python3 scripts/seqGet.py "$taxlist" -o "$sp_out"

#remove GCA if GCF anotation is available
for folder in "$sp_out"/*/; do
    #ensure we are looking at a directory to avoid errors with loose files
    [[ -d "$folder" ]] || continue

    #identify GCA and GCF
    gca_dirs=("$folder"GCA_*)
    gcf_dirs=("$folder"GCF_*)

    #check if at least one directory exists for both prefixes
    #check if the first element is a real directory
    if [[ -d "${gca_dirs[0]}" && -d "${gcf_dirs[0]}" ]]; then
        echo "Match found in: $(basename "$folder")"
        echo "Removing: ${gca_dirs[@]}"
        
        #remove gca
        rm -rf "${gca_dirs[@]}"
    else
        echo "Skipping: $(basename "$folder") (GCA only or no match)"
    fi
done

#rm "$sp_out"/Durusdinium_trenchii* -rf
#rm "$sp_out"/Symbiodinium_natans* -rf
#rm "$sp_out"/Symbiodinium_necroappetens* -rf
#rm "$sp_out"/Brevolum_minutum* -rf
rm "$sp_out"/Cyanidioschyzon_merolae_strain_10d* -rf

#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
echo ">ENDING at $(date)"
