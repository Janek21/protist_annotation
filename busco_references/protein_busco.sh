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

#create storing folders and variablesi
#locate the genome FASTA, avoid .fai and prefer .fna(no gz)
raw_ref_gff=$(find ../data/species -type f -name "$species_name*.gff*" ! -name '*.tbi' ! -name '*.csi' ! -name '*.gzi' -path "*/$species_name*" 2>/dev/null | sort | head -n 1)
raw_ref_fa=$(find ../data/species -type f -name 'GC*.fna*' ! -name '*.fai' ! -name '*.gzi' -path "*/$species_name*" 2>/dev/null | sort | head -n 1)

#no reference annotation -> skip protein BUSCO gracefully (genome BUSCO still runs)
if [ -z "$raw_ref_gff" ]; then
    echo ">No reference annotation found for $species_name. Skipping protein BUSCO and gene/transcript counts."
    exit 0
fi

#accession prefix (GCA or GCF) from the matched assembly folder
acc_prefix=$(basename "$(dirname "$raw_ref_gff")" | cut -c1-3)
tmp_files="species/$species_name/files"
res_base="species/$species_name/busco_res"

rm -rf "$res_base/prot_L" "$res_base/prot_E"
mkdir -p "$res_base/prot_L" "$res_base/prot_E"
mkdir -p "$tmp_files"
mkdir -p summary/busco_lineage summary/busco_eukaryote summary/counts

##wath out for compression
#gff
if [[ "$raw_ref_gff" == *.gz ]]; then
    echo "GFF is compressed. Decompressing to tmp folder"
    base_gff=$(basename "${raw_ref_gff%.gz}")
    gunzip -c "$raw_ref_gff" > "$tmp_files/$base_gff"
    ref_gff=$(realpath "$tmp_files/$base_gff")
else
    ref_gff=$(realpath "$raw_ref_gff")
fi
#fa
if [[ "$raw_ref_fa" == *.gz ]]; then
    echo "FASTA is compressed. Decompressing to tmp folder"
    base_fa=$(basename "${raw_ref_fa%.gz}")
    gunzip -c "$raw_ref_fa" > "$tmp_files/$base_fa"
    ref_fa=$(realpath "$tmp_files/$base_fa")
else
    ref_fa=$(realpath "$raw_ref_fa")
fi

##run busco

#resolve taxon id from the organism-name column (field 5), anchored so a loose
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

#get lineage
busco_lineage=$(python3 scripts/get_busco_db.py -e "ibdyjsayzcllkyvjkc@nespf.com" -t "$taxonID" -b "$busco_db/file_versions.tsv" -v odb12)
echo "BUSCO lineage for $taxonID is $busco_lineage"

#fix gff if needed(replace ? in strand column for .)
questionPresence_ref=$(cut -f7 $ref_gff|grep -Fx "?"|wc -l)

if [ "$questionPresence_ref" -ne 0 ]; then
        echo "Replacing ? strand symbol in reference annotation"
        awk -F'\t' 'BEGIN{OFS="\t"} {$7=gensub(/\?/, ".", "g", $7); print}' "$ref_gff" > "$tmp_files/newRef_${acc_prefix}_${species_name}.gff"
        #replace reference variable
        ref_gff=$(realpath "$tmp_files/newRef_${acc_prefix}_${species_name}.gff")
fi

#count gene/transcript models with gffread. --keep-genes normalises the reference
#into gene + transcript records (real gene features preserved, one synthesised per
#id when absent); counting the normalised col3 types catches every biotype
#(mRNA, tRNA, rRNA, ncRNA, pseudogene, ...) so the reference baseline is counted
#the same way as the IsoQuant/LyRic predictions it is compared against.
read -r gene_count transcript_count < <(
    { gffread "$ref_gff" --keep-genes -o - 2>/dev/null || true; } | awk -F'\t' '
        /^#/ { next }
        $3 ~ /^([A-Za-z_]*gene)$/                { g++; next }
        $3 ~ /^(transcript|mRNA|[A-Za-z_]*RNA)$/ { t++ }
        END { print g + 0, t + 0 }'
)
echo "$gene_count" > "summary/counts/${species_name}_${taxonID}_gc.txt"
echo "$transcript_count" > "summary/counts/${species_name}_${taxonID}_tc.txt"
echo "Gene models: $gene_count | Transcript models: $transcript_count"

#resolve the correct NCBI translation table for this taxon.
#AGAT to use the right table
gcode=$(python3 scripts/get_genetic_code.py -e "ibdyjsayzcllkyvjkc@nespf.com" -t "$taxonID" 2>/dev/null)
if ! [[ "$gcode" =~ ^[0-9]+$ ]]; then
    echo ">Could not resolve genetic code for taxon $taxonID; defaulting to table 1."
    gcode=1
fi
echo "Translation table for $taxonID: $gcode"

#per-task AGAT config so parallel array jobs don't collide on agat_config.yaml
agat_cfg="$tmp_files/agat_${sp}_${SLURM_ARRAY_TASK_ID:-$$}.yaml"
agat config --expose --output "$agat_cfg" >/dev/null 2>&1

#extract proteins from annotation with AGAT, using the resolved table #if error go in the statement
if ! agat_sp_extract_sequences.pl -g "$ref_gff" -f "$ref_fa" -t cds -p --table "$gcode" --config "$agat_cfg" -o "$tmp_files/prot_$sp.fa"; then
    echo ">ERROR: AGAT protein extraction failed for $species_name." >&2
    rm -f "$agat_cfg"
    exit 1
fi

#clean up AGAT
rm -f "$agat_cfg" ./*.agat.log "${ref_gff}.agat.log" "${ref_fa}.index" "${ref_fa}.gdx"

#coding transcripts = unique reference transcripts that yielded a CDS protein.
#AGAT writes one record per CDS-bearing transcript; dedupe ids to be safe against
#duplicate names (the same reason seqkit rename runs below).
coding_tx=$( { grep '^>' "$tmp_files/prot_$sp.fa" 2>/dev/null || true; } \
    | awk '{sub(/^>/,"",$1); print $1}' | sort -u | wc -l)
echo "$coding_tx" > "summary/counts/${species_name}_${taxonID}_ct.txt"
echo "Coding transcripts: $coding_tx"

#fix duplciate sequence names
seqkit rename -n "$tmp_files/prot_$sp.fa" > "$tmp_files/ND_prot_${sp}.fa"

#Run busco - by-lineage
busco -m protein -i "$tmp_files/ND_prot_$sp.fa" --download_path "$busco_db" -l "$busco_lineage" -c "$cpus" -f --out_path "$res_base" -o prot_L --tar
#Run busco - eukaryote lineage
busco -m protein -i "$tmp_files/ND_prot_$sp.fa" --download_path "$busco_db" -l eukaryota_odb12 -c "$cpus" -f --out_path "$res_base" -o prot_E --tar

#short summaries: rename and link to summary folders
mv "$res_base/prot_L"/*json "$res_base/prot_L/${species_name}_${taxonID}_prot_Lbusco.json"
ln -vf "$res_base/prot_L/${species_name}_${taxonID}_prot_Lbusco.json" summary/busco_lineage
mv "$res_base/prot_E"/*json "$res_base/prot_E/${species_name}_${taxonID}_prot_Ebusco.json"
ln -vf "$res_base/prot_E/${species_name}_${taxonID}_prot_Ebusco.json" summary/busco_eukaryote
#busco --plot summary/busco_lineage



