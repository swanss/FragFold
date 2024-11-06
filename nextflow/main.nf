// Import processes from module
include {getName; build_msa; build_msa as build_msa_fragparent; process_msa; colabfold; create_summary_csv; predict_peaks} from './modules'

// Declare syntax version
nextflow.enable.dsl=2

// Define workflow
workflow {
    build_msa(file(params.protein_query_seq),true)
    
    protein_msa = build_msa.out.a3m
    protein_name = getName(params.protein_query_seq,params.protein_name)
    if (params.heteromeric_mode) {
        build_msa_fragparent(file(params.fragment_query_seq),build_msa.out.done)
        fragment_msa = build_msa_fragparent.out.a3m
        fragment_parent_name = getName(params.fragment_query_seq,params.fragment_parent_name)
    } else {
        fragment_msa = protein_msa
        fragment_parent_name = protein_name
    }
    process_msa(fragment_msa,
                params.fragment_ntermres_start,
                params.fragment_ntermres_final,
                params.fragment_length,
                protein_msa,
                params.protein_nterm_res,
                params.protein_cterm_res,
                params.protein_copies)
                | flatten
                | colabfold
    create_summary_csv(colabfold.out.log | collect , 
                       colabfold.out.pdb | collect ,
                       protein_name,
                       fragment_parent_name,
                       params.experimental_data,
                       params.job_name,
                       params.contact_distance_cutoff)
    predict_peaks(create_summary_csv.out.csv,
                      params.job_name,
                      params.n_contacts,
                      params.n_weighted_contacts,
                      params.iptm,
                      params.contact_distance_cluster_cutoff,
                      params.cluster_peaks_frac_overlap)
}