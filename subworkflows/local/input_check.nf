//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastx_channel(it) }
        .set { files }

    emit:
    files                                     // channel: [ val(meta), [ files ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq, fasta ] ]
def create_fastx_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.contig     = row.contig
    meta.size       = row.target_size
    meta.replicon   = row.replicon

    // add paths of the fastx files to the meta map
    def fastx_meta = []
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> FastQ file does not exist!\n${row.fastq}"
    }
    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> FastA file does not exist!\n${row.fasta}"
    }
    fastx_meta = [ meta, [ file(row.fastq), file(row.fasta) ] ]
    return fastx_meta
}
