#!/usr/bin/env python3

import os 
import json
import random
import requests
import argparse
from tqdm import tqdm

# ------------------------------------------------------------------------------
# some known taxon IDs (you can expand this dictionary as needed with the add function)
# and the direction to save the id's
# ------------------------------------------------------------------------------

taxon_map = "taxon_map.json" 
ids_dir = "./ids"

#------------------------------------------------------------------------------
# functions :p
#------------------------------------------------------------------------------

def load_taxon_map(filename=taxon_map) -> dict:
    if not os.path.exists(filename):
        return {}
    with open(filename, "r", encoding="utf-8") as f:
        return json.load(f)

# add new taxon id's 

def save_taxon_map(taxon_map: dict, filename=taxon_map):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(taxon_map, f, indent=2)


def download_uniprot_ids(species: str, taxon_id: str, cache_file: str) -> list:

    if species == "random":
        exit("Species cannot be 'random' for download_uniprot_ids")

    url = "https://rest.uniprot.org/uniprotkb/stream"
    # what the end it's suppoused to look like
    # "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_alphafolddb%2Cxref_pdb%2Cxref_ensembl&format=tsv&query=%28organism_id%3A9606%29"

    # add more fields here if you want. don't forget about the xref_ if they come from other databases
    #fields = "accession,alphafolddb,pdb,ensembl"
    
    params = {
        "query": f"organism_id:{taxon_id}",
        "fields": "accession,xref_alphafolddb,xref_pdb,xref_ensembl",
        "format": "tsv"
    }
    
    print(f"Fetching IDs for {species} (organism_id={taxon_id})")
    r = requests.get(url, params=params,  stream=True)
    r.raise_for_status()


    total_size = int(r.headers.get("content-length", 0))  # in bytes
    block_size = 1024  # 1KB chunks

    # the bar is not showing in my machine, dk why, but it shows the size :P
    with tqdm(total=total_size, unit="B", unit_scale=True, desc=f"Downloading {species}") as pbar:
        with open(cache_file, "wb") as f:
            for chunk in r.iter_content(block_size):
                pbar.update(len(chunk))
                f.write(chunk)
                continue

    # if total_size != 0 and pbar.n != total_size:
    #     raise RuntimeError("Could not download file (incomplete)")

    # text_data = r.text.strip()
    # lines = text_data.split("\n")

    # if len(lines) < 2:
    #     print("Warning: got fewer than 2 lines from UniProt (header only).")
    
    # save to cache_file exactly as received (including header line)
    with open(cache_file, "rw") as f:
            text_data = r.text.strip()
            lines = text_data.split("\n")
            f.write(text_data + "\n")
            if len(lines) < 2:
                print("Warning: got fewer than 2 lines from UniProt (header only).")




def parse_tsv_for_db(tsv_file: str, db: str):

    # in the TSV header columns are in this order:
    # Entry,AlphaFoldDB,PDB,Ensembl
    # if you want to add more stuff, just add it to the request and update this logic :)
    # direct mapping for clarity:

    col_index_map = {
        "uniprot":  0,  # the "Entry" column
        "af":       1,  # the "AlphaFoldDB" column
        "pdb":      2,  # the "PDB" column
        "enst":     3,  # the "Ensembl" column
    }
    
    idx = col_index_map.get(db)
    if idx is None:
        raise ValueError(f"Unsupported db type '{db}'. Must be one of {list(col_index_map.keys())}")

    with open(tsv_file, "r", encoding="utf-8") as f:
        lines = f.read().strip().split("\n")
    
    if len(lines) < 2:
        # either empty file or only header
        return []
    
    # header
    data_lines = lines[1:]
    
    ids = []
    for line in data_lines:
        columns = line.split("\t")
        if len(columns) < 4:
            continue
        
        col_data = columns[idx].strip()  # this is the relevant cell
        if not col_data:
            # no data in this column for this row
            continue
        
        # cells might have multiple IDs separated by semicolons, e.g. "ENST00000380596.10;"
        # split by ";" and filter out empty strings
        splitted = [x.strip() for x in col_data.split(";") if x.strip()]
        
        # extend our global list
        ids.extend(splitted)
    
    return ids


def get_random_uniprot_id_for_species(species: str, db: str,taxon_map: dict) -> str:

    if species == "random":
        if not taxon_map:
            raise ValueError("No species in taxon_map.json, cannot pick a random one!")
        species = random.choice(list(taxon_map.keys()))
        print(f"Selected random species: {species}")

    if species not in taxon_map:
        raise ValueError(f"Species '{species}' not found in local taxon map. "
                         f"Use 'radid add <species> <taxon_id>' to add it.")

    taxon_id = taxon_map[species]

    #  cache file path
    os.makedirs(ids_dir, exist_ok=True)
    cache_file = os.path.join(ids_dir, f"{species}.tsv")
    
    # load from cache if present
    if os.path.isfile(cache_file):
        print(f"Loading IDs from cache file: {cache_file}")
        with open(cache_file, "r", encoding="utf-8") as f:
            lines = f.read().strip().split("\n")
        # if empty or only header, we might need to re-download
        if len(lines) < 2:
            print("Cache file is empty or invalid; re-downloading...")
            accessions = download_uniprot_ids(species, taxon_id, cache_file)
        else:
            accessions = lines[1:]  # skip header
    else:
        print(f"No cache found for {species}; downloading from UniProt...")
        accessions = download_uniprot_ids(species, taxon_id, cache_file)

    relevant_ids = parse_tsv_for_db(cache_file, db)

    # no need for this now!
    # if not accessions:
    #     raise ValueError(f"No valid IDs found for species '{species}' (taxon_id={taxon_id}).")

    return random.choice(relevant_ids)

#------------------------------------------------------------------------------
# main and help
#------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Get a random ID (UniProt, AlphaFold, PDB, Ensembl ENST) or manage species->taxon map.\n"
    )
    parser.add_argument("args", nargs="*", help="""
Examples:
  radid add <species> <taxon_id>     # Add new species/taxon mapping
  radid list                         # List local species->taxon mappings
  radid homo_sapiens uniprot         # Random UniProt from Homo sapiens
  radid mus_musculus af              # Random AlphaFold from mouse
  radid mus_musculus pdb             # Random PDB ID from mouse
  radid homo_sapiens enst            # Random Ensembl transcript
  radid random enst                  # Random species, then random transcript
""".strip())

    parsed = parser.parse_args()
    if not parsed.args:
        parser.print_help()
        return

    # species->taxon map
    taxon_map = load_taxon_map()



    if parsed.args[0].lower() == "add":
        if len(parsed.args) != 3:
            print("Usage: radid add <species_name> <taxon_id>")
            return
        species_name, taxon_id = parsed.args[1], parsed.args[2]
        taxon_map[species_name] = taxon_id
        save_taxon_map(taxon_map)
        print(f"Added/updated: {species_name} -> {taxon_id}")
        return

    if parsed.args[0].lower() == "list":
        print("Local files:\n")
        for i in os.listdir(ids_dir):

            print(i[:-4])
        return

    if len(parsed.args) > 2:
        print("ERROR: invalid usage. Example: radid homo_sapiens uniprot")
        return

    species, db = parsed.args
    db = db.lower()
    
    #  { "uniprot", "af", "pdb", "enst" }
    valid_dbs = ["uniprot", "af", "pdb", "enst"]
    if db not in valid_dbs:
        print(f"ERROR: unknown db '{db}'. Choose from {valid_dbs}")
        return

    try:
        random_id = get_random_uniprot_id_for_species(species, db, taxon_map)
        print(random_id)
    except Exception as e:
        print(f"ERROR: {e}")

if __name__ == "__main__":
    main()