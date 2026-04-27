#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// core-microbiome pipeline — generates heatmap prevalence/detection data
// ============================================================================

// ---- Helper: exclude filter argument ----------------------------------------
def exclude_arg(String col, String vals) {
    def parts = []
    if (col  != null && col  != "") parts << "--exclude_column '${col}'"
    if (vals != null && vals != "") parts << "--exclude_values '${vals}'"
    parts.join(" ")
}

// ============================================================================
// Process: CORE_MICROBIOME
// ============================================================================
process CORE_MICROBIOME {
    tag "${params.label}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path feature_table
    path meta_table
    path taxonomy_table

    output:
    path "Core_heatmap_*.csv", emit: heatmap_csv
    path "Core_members_*.csv", emit: members_csv

    script:
    def tax_arg  = (taxonomy_table.name != 'NO_TAXONOMY') ? "--taxonomy_table ${taxonomy_table}" : ""
    def excl_arg = exclude_arg(params.exclude_column, params.exclude_values)
    def det_min  = (params.detection_range_min != null) ? params.detection_range_min : -1
    def det_max  = (params.detection_range_max != null) ? params.detection_range_max : -1

    """
    Rscript ${projectDir}/src/R/core_microbiome.R \
        --feature_table   "${feature_table}" \
        --input_format    "${params.input_format}" \
        ${tax_arg} \
        --meta_table      "${meta_table}" \
        --output_dir      "." \
        --taxon_rank     "${params.taxon_rank}" \
        --label           "${params.label}" \
        --min_library_size ${params.min_library_size} \
        --prevalence_min   ${params.prevalence_min} \
        --what_detection  "${params.what_detection}" \
        --detection_steps  ${params.detection_steps} \
        --detection_range_min ${det_min} \
        --detection_range_max ${det_max} \
        ${excl_arg} \
        ${params.group ? "--group '${params.group}'" : ""} \
        ${params.type  ? "--type  '${params.type}'"  : ""}
    """
}

// ============================================================================
// Process: Merge all CSVs into a single Parquet file
// ============================================================================
process MERGE_PARQUET {
    tag "${params.label}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path csvs

    output:
    path "core_microbiome_${params.label}.parquet"

    script:
    """
    Rscript ${projectDir}/src/R/merge_parquet.R \\
        --label      '${params.label}' \\
        --output_dir '.'
    """
}

// ============================================================================
// Workflow
// ============================================================================
workflow {

    // ---- Validate required params ------------------------------------------
    if (!params.feature_table) {
        error "params.feature_table must be set"
    }
    if (!params.meta_table) {
        error "params.meta_table must be set"
    }

    // ---- Channels ----------------------------------------------------------
    ch_feature_table = Channel.fromPath(params.feature_table, checkIfExists: true)
    ch_meta_table    = Channel.fromPath(params.meta_table,    checkIfExists: true)
    ch_taxonomy = params.taxonomy_table
        ? Channel.fromPath(params.taxonomy_table, checkIfExists: true).first()
        : Channel.value(file('NO_TAXONOMY'))

    cm_out = CORE_MICROBIOME(ch_feature_table, ch_meta_table, ch_taxonomy)

    if (params.merge_parquet) {
        all_csvs = cm_out.heatmap_csv
            .mix(cm_out.members_csv)
            .collect()
        MERGE_PARQUET(all_csvs)
    }
}
