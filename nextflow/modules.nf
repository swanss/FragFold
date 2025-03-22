def getName(file_path,file_name) {
    if ( file_name == "" ) {
        return file(file_path).baseName
    } else {
        return file_name
    }
}


// Define each process
process build_msa {
    label 'cpu_network'
    cache 'lenient'

    input:
        path query_seq
        val ready
    
    output:
        path 'mmseqs2_a3m/*.a3m', emit: a3m
        val true, emit: done

    shell:
    '''
    export PATH="!{colabfold_dir}/bin:$PATH"
    python !{repo_dir}/scripts/mmseqs2.py --query !{query_seq}
    '''
}

process process_msa {
    label 'cpu'
    cache 'lenient'

    input:
        path a3m
        val fragment_ntermres_start
        val fragment_ntermres_final
        val fragment_length
        val protein_a3m_input
        val protein_ntermres
        val protein_ctermres
        val protein_copies
        val fragment_parent_name
        val protein_name

    output:
        path '*/*.a3m'

    shell:
    '''
    # set arg vars
    if [[ !{fragment_single_sequence} == true ]]; then arg1="--fragment_single_sequence"; else arg1=""; fi
    if [[ !{fragment_shuffle_sequence} == true ]]; then arg2="--fragment_shuffle_sequence"; else arg2=""; fi

    python !{repo_dir}/fragfold/create_fragment_msa.py \
        --fragment_a3m_input !{a3m} \
        --fragment_ntermres_start !{fragment_ntermres_start} \
        --fragment_ntermres_final !{fragment_ntermres_final}  \
        --fragment_length !{fragment_length} \
        --protein_a3m_input !{protein_a3m_input} \
        --protein_ntermres !{protein_ntermres} \
        --protein_ctermres !{protein_ctermres} \
        --protein_copies !{protein_copies} \
        --fragment_parent_name !{fragment_parent_name} \
        --protein_name !{protein_name} \
        $arg1 \
        $arg2 \
    '''
}

process colabfold {
    label 'gpu'
    cache 'lenient'
    publishDir "$colabfold_outdir/${a3m_concat.baseName}", overwrite: true

    input:
        path a3m_concat

    output:
        path 'log.txt', emit: log
        path '*_unrelaxed_rank_00?_*.pdb', emit: pdb
        path '*.png', emit: png

    shell:
    '''
    export PATH="!{colabfold_dir}/bin:$PATH"
    colabfold_batch !{a3m_concat} . \
        --data !{alphafold_params_dir} \
        --model-type !{model_type} \
        --pair-mode !{pair_mode} \
        --num-models !{num_models}
    '''
}

process create_summary_csv {
    label 'cpu'
    cache 'lenient'
    publishDir '.', saveAs: { csv -> "${output_name}__${csv}" }, overwrite: true
    publishDir path: "$peakprediction_outdir", pattern: '*.png', overwrite: true

    input:
        path 'log_file_*.txt'
        path pdb_file
        val protein_name
        val fragment_parent_name
        path experimental_data
        val output_name
        val contact_distance_cutoff

    output:
        path 'colabfold_predictions.csv', emit: csv
        path '*.png', emit: png

    script:
    shell:
    '''
    exp_data=!{experimental_data}
    filename=$(basename -- "$exp_data")
    EXP_DATA_ARG=""
    if [[ $filename != "NO_FILE" ]]; then EXP_DATA_ARG="--experimental_data "$exp_data; fi
    python !{repo_dir}/fragfold/colabfold_process_output.py \
        --predicted_pdbs !{pdb_file} \
        --confidence_logs log_file_*.txt \
        --protein_name !{protein_name} \
        --fragment_parent_name !{fragment_parent_name} \
        --contact_distance_cutoff !{contact_distance_cutoff} \
        --generate_plots \
        $EXP_DATA_ARG
    '''
}

process create_summary_csv_fromjson {
    label 'cpu'
    cache 'lenient'
    publishDir '.', saveAs: { csv -> "${output_name}__${csv}" }, overwrite: true
    publishDir path: "$peakprediction_outdir", pattern: '*.png', overwrite: true

    input:
        path json_file
        path experimental_data
        val output_name
        val contact_distance_cutoff

    output:
        path 'colabfold_predictions.csv', emit: csv
        path '*.png', emit: png

    shell:
    '''
    python !{repo_dir}/fragfold/colabfold_process_output.py \
        --import_json !{json_file} \
        --experimental_data !{experimental_data} \
        --contact_distance_cutoff !{contact_distance_cutoff} \
        --generate_plots
    '''
}

process predict_peaks {
    label 'cpu_small'
    cache 'lenient'
    publishDir path: "$peakprediction_outdir", pattern: '*.csv', saveAs: { x -> "${output_name}__${x}" }, overwrite: true
    publishDir path: "$peakprediction_outdir", pattern: 'cluster_info/*/*_mergedpeaks.png', saveAs: { x -> "${output_name}__${file(x).name}" }, overwrite: true

    input:
        path csv
        val output_name
        val n_contacts
        val n_weighted_contacts
        val iptm
        val contact_distance_cluster_cutoff
        val cluster_peaks_frac_overlap

    output:
        path '*.csv', optional: true, emit: csv
        path 'cluster_info/*/*_mergedpeaks.png', emit: png

    shell:
    '''
    python !{repo_dir}/fragfold/predict_alphafold_peaks.py \
        --colabfold_data_csv !{csv} \
        --n_contacts !{n_contacts} \
        --n_weighted_contacts !{n_weighted_contacts} \
        --iptm !{iptm} \
        --contact_distance_cluster_cutoff !{contact_distance_cluster_cutoff} \
        --cluster_peaks_frac_overlap !{cluster_peaks_frac_overlap} \
        --verbose
    '''
}
