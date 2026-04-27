#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// core-microbiome pipeline — generates heatmap prevalence/detection data
// ============================================================================

// ---- Helper: resolve optional taxonomy_table argument ----------------------
def taxonomy_arg(String path) {
    (path != null && path != "") ? "--taxonomy_table ${path}" : ""
}

def exclude_arg(String col, String vals) {
    def parts = []
    if (col  != null && col  != "") parts << "--exclude_column '${col}'"
    if (vals != null && vals != "") parts << "--exclude_values '${vals}'"
    parts.join(" ")
}

def groups_arg(String col, String paste_cols) {
    def parts = []
    if (col        != null && col        != "") parts << "--groups_column '${col}'"
    if (paste_cols != null && paste_cols != "") parts << "--groups_paste_columns '${paste_cols}'"
    parts.join(" ")
}

def type_arg(String col, String col2, String levels2) {
    def parts = []
    if (col    != null && col    != "") parts << "--type_column '${col}'"
    if (col2   != null && col2   != "") parts << "--type2_column '${col2}'"
    if (levels2 != null && levels2 != "") parts << "--type2_levels '${levels2}'"
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
    val  taxonomy_table_path

    output:
    path "Core_heatmap_*.csv", emit: heatmap_csv
    path "Core_members_*.csv", emit: members_csv

    script:
    def tax_arg   = taxonomy_arg(taxonomy_table_path)
    def excl_arg  = exclude_arg(params.exclude_column, params.exclude_values)
    def grp_arg   = groups_arg(params.groups_column, params.groups_paste_columns)
    def tp_arg    = type_arg(params.type_column, params.type2_column, params.type2_levels)
    def det_min   = (params.detection_range_min != null) ? params.detection_range_min : -1
    def det_max   = (params.detection_range_max != null) ? params.detection_range_max : -1

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
        ${grp_arg} \
        ${tp_arg}
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
    ch_taxonomy      = Channel.value(
        (params.taxonomy_table != null && params.taxonomy_table != "")
        ? params.taxonomy_table
        : ""
    )

    cm_out = CORE_MICROBIOME(ch_feature_table, ch_meta_table, ch_taxonomy)

    if (params.merge_parquet) {
        all_csvs = cm_out.heatmap_csv
            .mix(cm_out.members_csv)
            .collect()
        MERGE_PARQUET(all_csvs)
    }
}
