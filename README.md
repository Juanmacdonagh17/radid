# radid

`radid` is a small command-line tool to fetch random IDs from UniProt proteomes.

It can return random:

- **UniProt accessions**
- **AlphaFoldDB IDs**
- **PDB IDs**
- **Ensembl transcript IDs (ENST)**

It also lets you manage a local `species -> taxon_id` mapping and caches downloaded UniProt TSV files locally for reuse.

---

## Features

- Get random IDs from a given species
- Supports `uniprot`, `af`, `pdb`, and `enst`
- Accepts arguments in either order:
  - `radid homo_sapiens uniprot 5`
  - `radid uniprot homo_sapiens 5`
- Local cache of downloaded TSV files in `./ids`
- Local editable taxon mapping via `taxon_map.json`
- Random species mode

---

## Installation

Clone the repo, install the required packages and chmodit:

```bash
pip install requests tqdm
chmod +x radid.py
```


Examples:

```bash
radid add arabidopsis_thaliana 3702
radid arabidopsis_thaliana uniprot 5
radid pdb homo_sapiens 3
radid random af 1
```
