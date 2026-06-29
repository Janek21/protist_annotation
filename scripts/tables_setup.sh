#!/usr/bin/env bash

#clear folder
rm protist_table/*html
rm protist_table/*tsv
rm protist_table/*txt

#Protist tables
python protist_table/protist_survey.py --email luejbpvzbzmvbqwnfb@gonrr.net --api-key $NCBI_API_KEY --html -o protist_table/strain_protists.tsv --no-assemblies-only
echo "Extracted all strains with RNA-seq"
python protist_table/protist_survey.py --email luejbpvzbzmvbqwnfb@gonrr.net --api-key $NCBI_API_KEY --html -o protist_table/protists.tsv
echo "Extracted all strains with RNA-seq and an assembly"

#Extraction
grep "long-read" protist_table/protists.tsv|cut -f4 | sort -u|grep -v taxid > data/protists_taxid.txt

#Longread info
cd data
python3 scripts/fetch_protist_rnaseq.py --reads long -o longread_protists.tsv --taxids-file protists_taxid.txt

cut -f3 longread_protists.tsv |sort -u > true_taxons.txt
