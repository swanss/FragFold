// Simple script to test submitting jobs in a SLURM env

// Declare syntax version
nextflow.enable.dsl=2

// Define parameters
params.job_name = "test"
params.query_seq = "/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/example/ftsZ.fasta"
params.fragment_ntermres_start = 160
params.fragment_ntermres_final = 170
params.fragment_length = 30
params.protein_nterm_res=10
params.protein_cterm_res=316
params.protein_copies=1

// Define each process
process build_msa {
    label 'cpu_network'

    input:
        path query_seq
    
    output:
        path 'mmseqs2_a3m/*.a3m'

    shell:
    '''
    export PATH="!{colabfold_dir}/bin:$PATH"
    python !{repo_dir}/scripts/mmseqs2.py --query !{query_seq}
    '''
}

process process_msa {
    label 'cpu'

    input:
        path a3m
        val fragment_ntermres_start
        val fragment_ntermres_final
        val fragment_length
        val protein_a3m_input
        val protein_ntermres
        val protein_ctermres
        val protein_copies

    output:
        path '*/*.a3m'

    shell:
    '''
    python !{repo_dir}/fragfold/create_fragment_msa.py \
        --fragment_a3m_input !{a3m} \
        --fragment_ntermres_start !{fragment_ntermres_start} \
        --fragment_ntermres_final !{fragment_ntermres_final}  \
        --fragment_length !{fragment_length} \
        --protein_a3m_input !{protein_a3m_input} \
        --protein_ntermres !{protein_ntermres} \
        --protein_ctermres !{protein_ctermres} \
        --protein_copies !{protein_copies}
    '''
}

process colabfold {
    label 'gpu'

    input:
        path a3m_concat

    output:
        path 'data/log.txt', emit: log
        path 'data/*_unrelaxed_rank_00?_*.pdb', emit: pdb
        
    shell:
    '''
    export PATH="!{colabfold_dir}/bin:$PATH"
    colabfold_batch !{a3m_concat} data \
        --data !{alphafold_params_dir} \
        --model-type 'alphafold2_ptm' \
        --pair-mode 'unpaired'
    '''
}

process create_summary_csv {
    label 'cpu'
    publishDir '.', saveAs: { csv -> "$csv" } 

    input:
        path 'log_file_*.txt'
        path pdb_file

    output:
        path '*.csv'

    shell:
    '''
    python !{repo_dir}/fragfold/colabfold_process_output.py --predicted_pdbs !{pdb_file} --confidence_logs log_file_*.txt --full_protein eh --fragment_protein uh
    '''
}

// Define workflow
workflow {
    query_seq = file(params.query_seq)
    a3m = build_msa(query_seq)
    process_msa(a3m,
               params.fragment_ntermres_start,
               params.fragment_ntermres_final,
               params.fragment_length,
               a3m,
               params.protein_nterm_res,
               params.protein_cterm_res,
               params.protein_copies)
        | flatten
        | colabfold
    create_summary_csv(colabfold.out.log | collect , colabfold.out.pdb | collect )
}