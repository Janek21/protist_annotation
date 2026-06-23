#!/usr/bin/env python

import os
import time
import random
import csv
import argparse
from urllib.error import HTTPError
from Bio import Entrez

def load_busco_lineages(busco_database):
    with open(busco_database, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        return {row[0].strip() for row in reader}

def get_taxonomy(taxon_id, query_email):
    cache_dir = os.path.join(os.getenv('TMPDIR', '/tmp'), 'biopython_cache')
    os.makedirs(cache_dir, exist_ok=True)
    Entrez.local_cache = cache_dir
    Entrez.email = query_email 
    
    #slight random delay to avoid bursts when many jobs run in parallel.
    #with an API key NCBI allows 10 req/s (vs 3 anonymous), so we can wait less.
    r_delay = random.uniform(0.1, 0.5) if Entrez.api_key else random.uniform(1, 10)
    time.sleep(r_delay)
    #retry loop
    for attempt in range(5):
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
            records = Entrez.read(handle, validate=False)
            handle.close()
            lineage = records[0]["Lineage"].split("; ")
            return lineage
        #catch error
        except HTTPError as e:
            if e.code==429: #overload ncbi error
                timeWait=(2*attempt)+random.uniform(1,3)#timeout slightly random
                print(f"NCBI rate limit hit. Task sleeping for {timeWait}s")
                time.sleep(timeWait)

def get_busco_lineage(taxon_id, query_email, busco_lineages, odb_version):
    lineage = get_taxonomy(taxon_id, query_email)
    for rank in reversed(lineage):
        database = f"{rank.lower()}_{odb_version}"
        if database in busco_lineages:
            return database
    return None

def get_taxon_id(species_name, query_email):
    """
    Takes a species name string and fetches the NCBI Taxon ID.
    """
    Entrez.email = query_email
    try:
        # Use esearch to find the taxon ID for the given species name
        handle = Entrez.esearch(db="taxonomy", term=species_name)
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            return record["IdList"][0]
        else:
            return None
    except Exception as e:
        print(f"Error fetching Taxon ID: {e}")
        return None

def species2lineage(species_name, query_email, busco_database, odb_version):
    #get taxon id
    taxon_id=get_taxon_id(species_name, query_email)
    if not taxon_id:
        return f"Error: No Taxon ID found for species {species_name}"
    
    #get busco lineages text
    busco_lineages=load_busco_lineages(busco_database)
    #scan lineage text for taxon classification
    busco_lineage_output=get_busco_lineage(taxon_id, query_email, busco_lineages, odb_version)
    if not busco_lineage_output:
        return f"No suitable BUSCO lineage found for {species_name}. (TaxonID found={taxon_id})"

    return busco_lineage_output

def taxon2lineage(taxon_id, query_email, busco_database, odb_version):
    
    #get busco lineages text
    busco_lineages=load_busco_lineages(busco_database)
    #scan lineage text for taxon classification
    busco_lineage_output=get_busco_lineage(taxon_id, query_email, busco_lineages, odb_version)
    if not busco_lineage_output:
        return f"No suitable BUSCO lineage found for {taxon_id}"

    return busco_lineage_output

def main():
    parser = argparse.ArgumentParser(description='Fetch protein sequences from NCBI for a given taxon ID.')
    parser.add_argument("-e", '--email', type=str, required=True, help='Email address to use for NCBI Entrez.')
    parser.add_argument("-k", "--api_key", type=str, default=os.getenv("NCBI_API_KEY"),
                        help="NCBI API key (defaults to the NCBI_API_KEY env var). "
                             "Raises the rate limit from 3 to 10 requests/second.")

    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", '--taxon_id', type=int, required=False, help='The taxon ID to start the search.')
    group.add_argument("-s", '--species_name', type=str, required=False, help='The species name to start the search.')

    parser.add_argument("-b", "--busco_lineage_database", type= str, required=True, help="Path to the database of current BUSCO lineages")
    parser.add_argument("-v", "--odb_version", type=str, required=True, help="odb database version 10 or 12")
    
    args=parser.parse_args()

    #set credentials once (module-global, so it applies to every Entrez call)
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    #if by species
    if args.species_name:
        busco_lineage_output=species2lineage(args.species_name, args.email, args.busco_lineage_database, args.odb_version)

    elif args.taxon_id:
        busco_lineage_output=taxon2lineage(args.taxon_id, args.email, args.busco_lineage_database, args.odb_version)

    print(busco_lineage_output)

if __name__ == "__main__":
    main()
