import argparse
import json
import os
import shutil
import subprocess
import urllib.request

def main():
    parser = argparse.ArgumentParser(description="Download fragment files from ENCODE")
    parser.add_argument("sample_id", type=str, help="ENCODE ENCSR sample ID")
    parser.add_argument("tmp_base_dir", type=str, help="base temp dir for staging")
    parser.add_argument("out_path", type=str, help="Final output path for the fragments.tsv.gz file")
    args = parser.parse_args()
    
    # Look up URL for the fragment file
    metadata = encode_fetch_json(args.sample_id)
    url = [f["href"] for f in metadata["files"] if f["output_type"] == "fragments"]
    assert len(url) == 1
    url = url[0]
    
    # Download the fragment file
    tmp_dir = os.path.join(args.tmp_base_dir, args.sample_id)
    tar_path = os.path.join(tmp_dir, "fragments.tar.gz")
    os.makedirs(tmp_dir, exist_ok=True)
    subprocess.run(["curl", "-L", f"https://www.encodeproject.org/{url}", "-o", tar_path])

    # Untar the files and move them to the destination location
    subprocess.run(["tar", "-xzf", tmp_dir + "/fragments.tar.gz", "--directory", tmp_dir], check=True)
    fragment_dir = tmp_dir + "/encode_scatac_dcc_2/results/" + args.sample_id + "-1/fragments/"
    shutil.move(fragment_dir + "/fragments.tsv.gz", args.out_path)
    shutil.move(fragment_dir + "/fragments.tsv.gz.tbi", args.out_path + ".tbi")
    shutil.rmtree(os.path.join(tmp_dir, "encode_scatac_dcc_2"))
    os.remove(tar_path)

def encode_fetch_json(id):
    req = urllib.request.Request(f"https://www.encodeproject.org/{id}", headers={'accept': 'application/json'})
    resp = urllib.request.urlopen(req)
    return json.load(resp)

if __name__ == "__main__":
   main()