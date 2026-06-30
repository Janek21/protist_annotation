import subprocess
import zipfile
import shutil
import os
import sys
import re
import json
import argparse
import gzip

def download_genome_for_taxid(taxid, output_root):
    taxid_str = str(taxid).strip()
    
    if not os.path.exists(output_root):
        os.makedirs(output_root)
    
    # 1. FAST CHECK: Check if the folder already exists
    for d in os.listdir(output_root):
        if d == taxid_str or d.endswith(f"_{taxid_str}"):
            print(f"[{taxid}] Folder '{d}' already exists. Skipping...")
            return

    print(f"[{taxid}] Fetching metadata...")
    
    # 2. Attempt 1: Look for Reference Genome
    cmd_ref = ["datasets", "summary", "genome", "taxon", taxid_str, "--reference"]
    report = None
    is_fallback = False

    try:
        result = subprocess.run(cmd_ref, capture_output=True, text=True, check=True)
        data = json.loads(result.stdout)
        
        if 'reports' in data and data['reports']:
            report = data['reports'][0]
        else:
            # 2b. Attempt 2: FALLBACK - Get best available (Representative or otherwise)
            print(f"[{taxid}] No Reference found. Searching for best alternative...")
            cmd_all = ["datasets", "summary", "genome", "taxon", taxid_str]
            result_all = subprocess.run(cmd_all, capture_output=True, text=True, check=True)
            data_all = json.loads(result_all.stdout)
            
            if 'reports' in data_all and data_all['reports']:
                # NCBI usually lists the 'Representative' or highest quality first
                report = data_all['reports'][0]
                is_fallback = True
            else:
                print(f"[{taxid}] Error: No assemblies found at all for this TaxID.")
                return

    except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
        print(f"[{taxid}] Metadata error: {e}")
        return

    # 3. Parse Metadata
    accession = report['accession']
    org_name = report.get('organism', {}).get('organism_name', 'UnknownSpecies')
    org_name_safe = org_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
    
    category = "REFERENCE" if not is_fallback else "REPRESENTATIVE/OTHER"
    print(f"[{taxid}] Selected {category} assembly: {accession} ({org_name})")

    # 4. Folder Logic
    base_folder_name = None
    for d in os.listdir(output_root):
        if os.path.isdir(os.path.join(output_root, d)):
            d_name_no_taxid = re.sub(r'_\d+$', '', d)
            if d_name_no_taxid.replace(' ', '_').lower() == org_name_safe.lower():
                base_folder_name = d
                break
    
    if not base_folder_name:
        base_folder_name = f"{org_name_safe}_{taxid}"
        
    full_base_path = os.path.join(output_root, base_folder_name)
    assembly_folder = os.path.join(full_base_path, accession)
    
    # Check for existing assembly, including gzip-compressed FASTA files
    if os.path.exists(assembly_folder) and any(
        f.endswith('.fna') or f.endswith('.fna.gz') for f in os.listdir(assembly_folder)
    ):
        print(f"[{taxid}] Assembly already exists. Skipping.")
        return

    os.makedirs(assembly_folder, exist_ok=True)
    
    # 5. Download
    zip_filepath = os.path.join(assembly_folder, f"{accession}.zip")
    cmd_download = [
        "datasets", "download", "genome", "accession", accession,
        "--include", "genome", "--filename", zip_filepath
    ]

    try:
        print(f"[{taxid}] Downloading...")
        subprocess.run(cmd_download, capture_output=True, text=True, check=True)
        
        # 6. Extraction (write FASTA files directly as gzip-compressed .fna.gz)
        with zipfile.ZipFile(zip_filepath, 'r') as z:
            fasta_files = [f for f in z.namelist() if f.endswith('.fna')]
            for fasta_file in fasta_files:
                filename = os.path.basename(fasta_file) + ".gz"
                target_path = os.path.join(assembly_folder, filename)
                with z.open(fasta_file) as source, gzip.open(target_path, "wb") as target:
                    shutil.copyfileobj(source, target)
        
        os.remove(zip_filepath)
        print(f"[{taxid}] Success! Saved to {assembly_folder}\n")
        
    except Exception as e:
        print(f"[{taxid}] Download/Extraction failed: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download genomes with fallback support.")
    parser.add_argument("input", help="TaxID(s) or .txt file", nargs="+")
    parser.add_argument("-o", "--out", help="Output folder", default=".")
    args = parser.parse_args()

    taxids = []
    if len(args.input) == 1 and os.path.isfile(args.input[0]):
        with open(args.input[0], 'r') as f:
            taxids = [line.strip() for line in f if line.strip()]
    else:
        taxids = args.input

    for tid in taxids:
        download_genome_for_taxid(tid, args.out)
