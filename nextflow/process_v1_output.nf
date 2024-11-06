// Import processes from module
include {create_summary_csv_fromjson; predict_peaks} from './modules'

// Declare syntax version
nextflow.enable.dsl=2

params.json_file = "/data1/groups/keatinglab/swans/savinovCollaboration/FragFold/input/json/process_colabfold_output/colabfold115_output_combined.json"
params.experimental_data = "/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/input/inhibitory_data/Savinov_2022_inhib_peptide_mapping.csv"
params.job_name = "colabfold115_v1"
params.contact_distance_cutoff = 3.5
params.n_contacts = 3
params.n_weighted_contacts = 3
params.iptm = 0.3
params.contact_distance_cluster_cutoff = 0.6
params.cluster_peaks_frac_overlap = 0.7

// Define workflow
workflow {
    create_summary_csv_fromjson(
        params.json_file,
        params.experimental_data,
        params.job_name,
        params.contact_distance_cutoff)
    predict_peaks(create_summary_csv_fromjson.out.csv,
            params.job_name,
            params.n_contacts,
            params.n_weighted_contacts,
            params.iptm,
            params.contact_distance_cluster_cutoff,
            params.cluster_peaks_frac_overlap)
}