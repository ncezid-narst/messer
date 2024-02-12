// TDOD nanoq: implement options for outputting (e.g. options like -j for json format or -s for stats only and no reads


process NANOQ {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::nanoq=0.10.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoq:0.10.0--h031d066_2'}"

    input:
    tuple val(meta), path(reads)
    // val min_len // TODO nanoq: make these optional? as in include them in task.ext.args?
    // val max_len

    output:
    tuple val(meta), path("*fastq.gz")   , optional: true, emit: fastq
    tuple val(meta), path('nanoq.log')   , emit: log
    // tuple val(meta), path('report.txt')  , optional: true, emit: report_txt
    // tuple val(meta), path('report.json') , optional: true, emit: report_json
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    nanoq \\
        $args \\
        -i $reads \\
        -o ${prefix}._nanoq.fastq.gz \\
        >> nanoq.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed 's/nanoq //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    touch ${prefix}._nanoq.fastq.gz
    touch nanoq.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed 's/nanoq //g')
    END_VERSIONS
    """
}
