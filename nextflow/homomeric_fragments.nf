// Import processes from module
include {build_msa, process_msa, colabfold, create_summary_csv} from './modules'

// Declare syntax version
nextflow.enable.dsl=2

// Define parameters
params.job_name = "ftsZ_test"
params.query_seq = "/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/example/ftsZ.fasta"
params.fragment_ntermres_start = 160
params.fragment_ntermres_final = 170
params.fragment_length = 30
params.protein_nterm_res = 150
params.protein_cterm_res = 200
params.protein_copies = 1

// Define workflow
workflow {
    build_msa(file(params.protein_query_seq),true)
    process_msa(build_msa.out.a3m,
               params.fragment_ntermres_start,
               params.fragment_ntermres_final,
               params.fragment_length,
               build_msa.out.a3m,
               params.protein_nterm_res,
               params.protein_cterm_res,
               params.protein_copies)
        | flatten
        | colabfold
    protein_name = file(params.protein_query_seq).baseName
    create_summary_csv(colabfold.out.log | collect ,
                       colabfold.out.pdb | collect , 
                       protein_name,
                       protein_name)
}