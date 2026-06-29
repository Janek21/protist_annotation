#!/bin/bash

species_name="$1"
#if no 2nd argument is given, it uses /no_backup...
busco_db="${2:-/no_backup/rg/references/busco_downloads}"
specie_tsv="${3:-../data/longread_protists.tsv}"
#use slurm cpus ionly if in slurm
cpus="${SLURM_CPUS_PER_TASK:-2}"

#species shortname
sp=$(echo "$species_name"|cut -f2 -d"_")
sp_extra=$(echo "$species_name"|cut -f3 -d"_")
echo "$sp + $sp_extra"

#cativate busco conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate buscomania

mkdir -p "species/$species_name"

#create storing folders and variables
#locate the genome FASTA, avoid .fai and prefer .fna(no gz)
raw_ref_fa=$(find ../data/species -type f -name 'GC*.fna*' ! -name '*.fai' ! -name '*.gzi' -path "*/$species_name*" 2>/dev/null | sort | head -n 1)
tmp_files="species/$species_name/files"
res_base="species/$species_name/busco_res"

rm -rf "$res_base/gen_L" "$res_base/gen_E"
mkdir -p "$res_base/gen_L" "$res_base/gen_E"
mkdir -p "$tmp_files"
mkdir -p summary/busco_lineage summary/busco_eukaryote summary/counts

##wath out for compression
#fa
if [[ "$raw_ref_fa" == *.gz ]]; then
    echo "FASTA is compressed. Decompressing to tmp folder"
    base_fa=$(basename "${raw_ref_fa%.gz}")
    gunzip -c "$raw_ref_fa" > "$tmp_files/$base_fa"
    ref_fa=$(realpath "$tmp_files/$base_fa")
else
    ref_fa=$(realpath "$raw_ref_fa")
fi


##resolve taxon id from the organism-name column (field 5), anchored so a loose
#substring can't steal another species' taxID (e.g. "gracilis" hitting
#"neogracilis"). progressive: full folder name -> binomial -> fail loud.
as_words=$(echo "$species_name" | tr '_' ' ')
binom=$(echo "$as_words" | awk '{print $1, $2}')
echo "Resolving taxon for '$as_words' in $specie_tsv"

#1) most specific: full folder name == organism field, case-insensitive
taxonID=$(awk -F'\t' -v q="$as_words" 'BEGIN{IGNORECASE=1} $5==q {print $3}' "$specie_tsv" | sort -u)
#2) fall back to species level (genus + epithet) if the strain isn't listed
if [ "$(printf '%s' "$taxonID" | grep -c .)" -ne 1 ]; then
    taxonID=$(awk -F'\t' -v q="$binom" 'BEGIN{IGNORECASE=1} $5==q {print $3}' "$specie_tsv" | sort -u)
fi

#fail loudly on 0 or >1 distinct taxIDs rather than guessing
n_tax=$(printf '%s' "$taxonID" | grep -c .)
if [ "$n_tax" -ne 1 ]; then
    echo ">ERROR: could not unambiguously resolve taxon for '$species_name' (matched $n_tax distinct taxIDs: $(echo $taxonID | tr '\n' ' '))." >&2
    exit 1
fi
echo "TAXON IS: $taxonID"

#genome size = total assembly length (exact; sum of contig lengths, incl. N gaps).
#emitted here (genome BUSCO runs for every species) so even species without a
#reference annotation - which skip protein_busco.sh and its gene/transcript counts
#- still get a genome size. ref_fa is already decompressed above.
fai="${ref_fa}.fai"
if [ -s "$fai" ]; then
    genome_size=$(cut -f2 "$fai" | awk '{s+=$1} END{print s+0}')
else
    genome_size=$(awk '/^>/{next} {s+=length($0)} END{print s+0}' "$ref_fa")
fi
echo "$genome_size" > "summary/counts/${species_name}_${taxonID}_gs.txt"
echo "Genome size: ${genome_size} bp"

#get lineage
busco_lineage=$(python3 scripts/get_busco_db.py -e "tgtvbraqmlgkpsgmxi@vtmpj.com" -t "$taxonID" -b "$busco_db/file_versions.tsv" -v odb12)
echo "BUSCO lineage for $taxonID is $busco_lineage"

#Run busco for genomes - by-lineage
busco -m genome -i "$ref_fa" --download_path "$busco_db" -l "$busco_lineage" -c "$cpus" -f --out_path "$res_base" -o gen_L --tar
#Run busco for genomes - eukaryote lineage
busco -m genome -i "$ref_fa" --download_path "$busco_db" -l eukaryota_odb12 -c "$cpus" -f --out_path "$res_base" -o gen_E --tar

#short summaries: rename and link to summary folders
mv "$res_base/gen_L"/*json "$res_base/gen_L/${species_name}_${taxonID}_gen_Lbusco.json"
ln -vf "$res_base/gen_L/${species_name}_${taxonID}_gen_Lbusco.json" summary/busco_lineage
mv "$res_base/gen_E"/*json "$res_base/gen_E/${species_name}_${taxonID}_gen_Ebusco.json"
ln -vf "$res_base/gen_E/${species_name}_${taxonID}_gen_Ebusco.json" summary/busco_eukaryote
#busco --plot summary/busco_lineage



