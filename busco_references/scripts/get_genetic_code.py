#!/usr/bin/env python
"""Resolve the NCBI nuclear genetic code (translation table) for a taxon.

NCBI Taxonomy assigns every taxon a nuclear genetic code id (GCId), e.g. 1 for
the standard code and 6 for the ciliate nuclear code used by Paramecium and
Tetrahymena. That id is exactly the translation table to feed to a protein
extractor. Prints the integer code to stdout on success.
"""

import os
import time
import random
import argparse
from urllib.error import HTTPError
from Bio import Entrez


def get_genetic_code(taxon_id):
    cache_dir = os.path.join(os.getenv('TMPDIR', '/tmp'), 'biopython_cache')
    os.makedirs(cache_dir, exist_ok=True)
    Entrez.local_cache = cache_dir

    #slight random delay to avoid bursts when many jobs run in parallel
    r_delay = random.uniform(0.1, 0.5) if Entrez.api_key else random.uniform(1, 5)
    time.sleep(r_delay)

    #retry loop (mirrors get_busco_db.py for NCBI rate-limit handling)
    for attempt in range(5):
        try:
            handle = Entrez.efetch(db="taxonomy", id=str(taxon_id), retmode="xml")
            records = Entrez.read(handle, validate=False)
            handle.close()
            #nuclear genetic code: GeneticCode/GCId (MitoGeneticCode is separate)
            return int(records[0]["GeneticCode"]["GCId"])
        except HTTPError as e:
            if e.code == 429:  #ncbi overload
                time_wait = (2 * attempt) + random.uniform(1, 3)
                print(f"NCBI rate limit hit. Sleeping for {time_wait}s")
                time.sleep(time_wait)
            else:
                raise
    return None


def species_to_taxon(species_name):
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0] if record["IdList"] else None


def main():
    parser = argparse.ArgumentParser(
        description='Fetch the NCBI nuclear genetic code (translation table) for a taxon.')
    parser.add_argument("-e", "--email", type=str, required=True,
                        help="Email address to use for NCBI Entrez.")
    parser.add_argument("-k", "--api_key", type=str, default=os.getenv("NCBI_API_KEY"),
                        help="NCBI API key (defaults to the NCBI_API_KEY env var).")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--taxon_id", type=int, help="The taxon ID to look up.")
    group.add_argument("-s", "--species_name", type=str, help="The species name to look up.")

    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    taxon_id = args.taxon_id
    if args.species_name:
        taxon_id = species_to_taxon(args.species_name)
        if not taxon_id:
            raise SystemExit(f"Error: No Taxon ID found for species {args.species_name}")

    code = get_genetic_code(taxon_id)
    if code is None:
        raise SystemExit(f"Error: Could not resolve genetic code for taxon {taxon_id}")
    print(code)


if __name__ == "__main__":
    main()
