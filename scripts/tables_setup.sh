#!/usr/bin/env bash

#Protist tables
python protist_table/protist_survey.py --email luejbpvzbzmvbqwnfb@gonrr.net --api-key 25dd740e73e5a092d87db2f6e230e39eae09 --html -o protist_table/strain_protists.tsv --no-assemblies-only
python protist_table/protist_survey.py --email luejbpvzbzmvbqwnfb@gonrr.net --api-key 25dd740e73e5a092d87db2f6e230e39eae09 --html -o protist_table/protists.tsv

#Extraction
grep "long-read" protist_table/protists.tsv|cut -f4 | sort -u|grep -v taxid > data/protists_taxid.txt

#Longread info
cd data
python3 scripts/fetch_protist_rnaseq.py --reads long -o longread_protists.tsv --taxids-file protists_taxid.txt

cut -f3 longread_protists.tsv |sort -u > true_taxons.txt
