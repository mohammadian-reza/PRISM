import os
import gzip
import shutil
import urllib.request

TARGET_DIR = "processed/pdbs"
os.makedirs(TARGET_DIR, exist_ok=True)

def download_pdb_file(pdb_name, pdb_dir):
    final_pdb = f"{pdb_dir}/{pdb_name}.pdb"
    gz_file = f"{pdb_dir}/{pdb_name}.ent.gz"

    try:
        url = f"https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{pdb_name[:4].lower()}.ent.gz"
        response = urllib.request.urlopen(url)

        with open(gz_file, "wb") as fh:
            fh.write(response.read())

        with gzip.open(gz_file, "rb") as f_in, open(final_pdb, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_file)

        return True
    except Exception as e:
        print(f"PDB download failed for {pdb_name}: {e}")
        return False

def pdb_downloader():
    targets = []

    # Read pair list
    with open("input_pair_list.txt", "r") as filehnd:
        for line in filehnd:
            line = line.strip().split()
            if len(line) == 2:
                if len(line[0]) == 5 and len(line[1]) == 5:
                    targets.append(line[0])
                    targets.append(line[1])
                else:
                    print(f"Skipping target {line[0]} or {line[1]} because it is not a 5-letter PDB ID."
                          f"We only accept single chain PDB IDs.")
                    continue
    
    # Download and process PDBs
    for target in list(set(targets)):
        if not os.path.exists(f"{TARGET_DIR}/{target[:4].lower()}.pdb"):
            if not download_pdb_file(target[:4].lower(), TARGET_DIR):
                print(f"Failed to download PDB {target}")
                continue
        else:
            print(f"PDB {target} already exists")

    return list(set(targets))
