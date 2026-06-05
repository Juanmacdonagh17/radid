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
loaded_ids = {}

#------------------------------------------------------------------------------
# functions :p
#------------------------------------------------------------------------------

### function for taxon map and writting 

def load_taxon_map(filename=taxon_map) -> dict:
    if not os.path.exists(filename):
        return {}
    with open(filename, "r", encoding="utf-8") as f:
        return json.load(f)

def save_taxon_map(taxon_map: dict, filename=taxon_map):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(taxon_map, f, indent=2)


### functions for downloading and parsing proteomes 

def download_uniprot_ids(species: str, taxon_id: str, cache_file: str) -> list:

    if species == "random":
        exit("Species cannot be 'random' for download_uniprot_ids")

    url = "https://rest.uniprot.org/uniprotkb/stream"


    # add more fields here if you want. don't forget about the xref_ if they come from other databases

    params = {
        "query": f"organism_id:{taxon_id}",
        "fields": "accession,xref_alphafolddb,xref_pdb,xref_ensembl,xref_ensemblbacteria,xref_ensemblfungi,xref_ensemblprotists", # all of the ensembl fields
        "format": "tsv"
    }
    
    print(f"Fetching IDs for {species} (organism_id={taxon_id})")
    r = requests.get(url, params=params,  stream=True, timeout=30)
    r.raise_for_status()


    total_size = int(r.headers.get("content-length", 0))  # in bytes
    block_size = 1024  

    # the bar is not showing in my machine, dk why, but it shows the size :P
    with tqdm(total=total_size, unit="B", unit_scale=True, desc=f"Downloading {species}") as pbar:
        with open(cache_file, "wb") as f:
            for chunk in r.iter_content(block_size):
                pbar.update(len(chunk))
                f.write(chunk)
                continue

### function for size checking, this are not large files at all but just in case

def file_size_kb(path: str) -> float:
    return os.path.getsize(path) / 1024.0

def dir_size_kb(path: str) -> float:
    total = 0
    if not os.path.isdir(path):
        return 0.0
    for root, _, files in os.walk(path):
        for fn in files:
            total += os.path.getsize(os.path.join(root, fn))
    return total / 1024.0


### function for getting a random ID 

def parse_tsv_for_db(tsv_file: str, db: str):

    # in the TSV header columns are in this order:
    # Entry,AlphaFoldDB,PDB,Ensembl, EnsemblBacteria, etc
    # gotta try this
    # if you want to add more stuff, just add it to the request and update this logic :)
    # direct mapping for clarity:

    col_index_map = {
        "uniprot":  0,  # the "Entry" column
        "af":       1,  # the "AlphaFoldDB" column (this is the same as the uniprot, but I keep it because there are UP w/o AF)
        "pdb":      2,  # the "PDB" column
        "enst":     3#,  # the "Ensembl" column
        #"enst_bact": 4
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
  
    # number of columns
    for i in lines[0:1]:
        len_cols = len(i.split("\t"))
    # print(len_cols)
   

    ids = []
    for line in data_lines:
        columns = line.split("\t")
        #print(columns)
        # if len(columns) < 4: # this should match the length of the request fields!
        #     continue

        if db == "enst":
            #print("AA")
            #print(len_cols)
            
            possible_ensembl_cols = range(3,len_cols)
            found_data = None
            for c in possible_ensembl_cols:
                if c >= len(columns):      
                    break
                #print("BB")
                #print(c)
                val = columns[c].strip()
                #print(val)
                if val:
                    #print("CC")
                    found_data = val
                    break
            if not found_data:
                #print("DD")
                # no ensembl data in this row
                continue
            splitted = [x.strip() for x in found_data.split(";") if x.strip()]
            ids.extend(splitted)
        else:
            #print("BB")
            # for uniprot, af, pdb
            # if db not in col_index_map:
            #     raise ValueError(f"Unknown db '{db}'.")
            col_idx = col_index_map[db]
            if col_idx >= len(columns):     # <-- add: same guard for uniprot/af/pdb
                continue
            cell = columns[col_idx].strip()
            if not cell:
                continue
            splitted = [x.strip() for x in cell.split(";") if x.strip()]
            ids.extend(splitted)

    return ids

### function for getting a random ID 

def load_ids_for_species(species: str, db: str, taxon_map: dict) -> list:
    if species == "random":
        # fix from before: pool map + cached-only species, not just the map
        cached = [f[:-4] for f in os.listdir(ids_dir)] if os.path.isdir(ids_dir) else []
        pool = sorted(set(taxon_map) | set(cached))
        if not pool:
            raise ValueError("No species available to pick from.")
        species = random.choice(pool)
        print(f"Selected random species: {species}")

    os.makedirs(ids_dir, exist_ok=True)
    cache_file = os.path.join(ids_dir, f"{species}.tsv")
    need_download = not os.path.isfile(cache_file)

    if need_download:
        if species not in taxon_map:
            raise ValueError(
                f"Species '{species}' not found in taxon map and no cache at {cache_file}."
                " Use 'radid add <species> <taxon_id>' or drop a TSV at ids/<species>.tsv."
            )
        taxon_id = taxon_map[species]
        print(f"No cache found for {species}; downloading from UniProt...")
        download_uniprot_ids(species, taxon_id, cache_file)
        print(f"Saved {file_size_kb(cache_file):.1f} KB. Total cache: {dir_size_kb(ids_dir):.1f} KB\n")
    elif (species, db) not in loaded_ids:
        print(f"Loading IDs from cache file: {cache_file}")

    if (species, db) not in loaded_ids:
        relevant_ids = parse_tsv_for_db(cache_file, db)
        if not relevant_ids:
            raise ValueError(f"No IDs found in {cache_file} for db='{db}'.")
        loaded_ids[(species, db)] = relevant_ids

    return loaded_ids[(species, db)]


def get_random_ids(species: str, db: str, taxon_map: dict, n: int) -> list:
    ids = load_ids_for_species(species, db, taxon_map)
    if n > len(ids):
        print(f"Only {len(ids)} unique IDs available; returning all.")
        n = len(ids)
    return random.sample(ids, n)


#------------------------------------------------------------------------------
# main and help
#------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Get  random IDs (UniProt, AlphaFold, PDB, Ensembl ENST) or manage species->taxon map.\n",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("args", nargs="*", help="""
Examples:
  radid add <species> <taxon_id>     # Add new species/taxon mapping                \n
  radid list                         # List local species->taxon mappings           \n
  
  radid example uniprot      10      # 10 random UniProt from example species (already loaded)      \n  
  radid uniprot example      10      # same as above, order does not matter :)     \n                      
  radid homo_sapiens uniprot 1       # 1 random UniProt from Homo sapiens           \n
  radid mus_musculus af      2       # 2 random AlphaFold IDs from mouse            \n
  radid mus_musculus pdb     3       # 3 random PDB IDs from mouse                  \n
  radid homo_sapiens enst    4       # 4 random Homo Sapiens Ensembl transcript IDs \n
  radid random enst          5       # 5 random species, then 5 random transcript   \n
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
        #print("\n")
        print("Local files available at: ", ids_dir)
        print("Downloaded species:")
        print("-----------------------------------------------------------------------------")
        files = os.listdir(ids_dir) if os.path.isdir(ids_dir) else []
        for i in files:
            print(i[:-4]) # so it does not print the .tsv each time, just the file name
        print("\n")
        print("-----------------------------------------------------------------------------")
        print("\nSpecies in taxon_map.json (added):") # they can differ (you first add and when you request the species it gets downloaded)
        for species, taxon_id in taxon_map.items():
            print(f"  {species} -> {taxon_id}")
        return
    

    if len(parsed.args) != 3:
        print("ERROR: invalid usage. Example: ./radid.py example uniprot 10")
        return

    # species, db, number = parsed.args
    # db = db.lower()
    # number = int(number)



    a1, a2, a3 = parsed.args

    try:
        number = int(a3)
    except ValueError:
        print(f"ERROR: count must be an integer, got '{a3}'")
        return
    if number < 1:
        print("ERROR: count must be >= 1")
        return
    

    valid_dbs = {"uniprot", "af", "pdb", "enst"}

    # arguments work both ways. less confusing
    x1 = a1.lower()
    x2 = a2.lower()

    if x1 in valid_dbs and x2 not in valid_dbs:
        db, species, number = x1, a2, int(a3)
    elif x2 in valid_dbs and x1 not in valid_dbs:
        species, db, number = a1, x2, int(a3)
    elif x1 in valid_dbs and x2 in valid_dbs:
        print(f"ERROR: both '{a1}' and '{a2}' look like dbs. Choose one of {sorted(valid_dbs)}.")
        return
    else:
        print(f"ERROR: neither '{a1}' nor '{a2}' looks like a db. Choose from {sorted(valid_dbs)}.")
        return
    
    valid_dbs = ["uniprot", "af", "pdb", "enst"]
    if db not in valid_dbs:
        print(f"Unknown db: '{db}'. Choose from {valid_dbs}")
        return

    try:
        for rid in get_random_ids(species, db, taxon_map, number):
            print(rid)
    except Exception as e:
        print(f"Error: {e}")
    

if __name__ == "__main__":
    main()