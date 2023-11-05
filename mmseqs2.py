from argparse import ArgumentParser
from pathlib import Path
import http.client as httplib
import sys

from colabfold.batch import (get_msa_and_templates,
                             get_queries,
                             msa_to_str)
from colabfold.utils import (DEFAULT_API_SERVER,
                             safe_filename)

def main(args):

    # First check for internet connection
    conn = httplib.HTTPSConnection("8.8.8.8", timeout=5)
    try: 
        conn.request("HEAD", "/")
    except:
        sys.exit(f"Unable to connect to internet, terminating.")
    conn.close()

    result_dir = "mmseqs2_a3m"
    result_dir = Path(result_dir)
    result_dir.mkdir(exist_ok=True)

    queries, is_complex = get_queries(args.query, None)

    print(f"Generating MSAs for {len(queries)} queries")
    for job_number, (raw_jobname, query_sequence, a3m_lines) in enumerate(queries):
        print(f"{raw_jobname}: {job_number}/{len(queries)}")
        jobname = safe_filename(raw_jobname)
        filepath = result_dir.joinpath(f"{jobname}.a3m")

        if filepath.is_file():
            print(f"Skipping '{jobname}', file already exists: {filepath}")
            continue

        msa_mode="mmseqs2_uniref_env"
        use_templates=False
        custom_template_path=None
        pair_mode="unpaired"
        host_url=DEFAULT_API_SERVER

        (unpaired_msa,
        paired_msa, 
        query_seqs_unique, 
        query_seqs_cardinality, 
        template_features) \
                        = get_msa_and_templates(jobname, 
                                                query_sequence, 
                                                result_dir, 
                                                msa_mode, 
                                                use_templates, 
                                                custom_template_path, 
                                                pair_mode, 
                                                host_url)

        msa = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
        result_dir.joinpath(f"{jobname}.a3m").write_text(msa)
    print("Done!")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--query",
        type=Path,
        help="Reads a directory of FASTA files, a single FASTA file or a csv file",
    )
    args = parser.parse_args()
    main(args)