process NANOQ {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::nanoq=0.10.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoq:0.10.0--h031d066_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*fastq.gz')   , emit: fastq
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    nanoq \\
        $args \\
        -i $reads \\
        -o ${prefix}._nanoq.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed 's/nanoq //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}._nanoq.fastq.gz
    touch nanoq.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed 's/nanoq //g')
    END_VERSIONS
    """
}
