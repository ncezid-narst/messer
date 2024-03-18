//
// extract contig as fasta format from fasta file
//

include { SEQTK_SEQ as SEQTK_SEQ1     } from '../../modules/nf-core/seqtk/seq/main'
include { SEQTK_SEQ as SEQTK_SEQ2     } from '../../modules/nf-core/seqtk/seq/main'
include { SEQTK_SUBSEQ                } from '../../modules/nf-core/seqtk/subseq/main'

workflow EXTRACT_CONTIG {
    take:
        fasta // channel: tuple val (meta), path (fasta)

    main:

        // fasta.view()

        ch_versions = Channel.empty()

        // convert contig name to file as required by seqtk subseq
        fasta
            .map { meta, fasta -> meta.contig }
            .collectFile(name: 'filter_list_${meta.contig}.txt', newLine: true)
        filter_list = Channel.fromPath('filter_list_${meta.contig}.txt')
        meta = fasta.map { meta, fasta -> meta }

        SEQTK_SEQ1(fasta)
        ch_versions = ch_versions.mix(SEQTK_SEQ1.out.versions)
        ch_extracted_fastq = SEQTK_SEQ1.out.fastx.collect{it[1]}
        ch_extracted_fastq.view()

        SEQTK_SUBSEQ(ch_extracted_fastq, filter_list)
        ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions )
        // SEQTK_SUBSEQ.out.view()
        // ch_fastq_contig = meta.map { meta -> [meta] }
        //     .concat(SEQTK_SUBSEQ.out.sequences)

        // SEQTK_SEQ2(
        //     meta.map { meta -> [meta] }
        //     .concat(SEQTK_SUBSEQ.out.sequences)
        // )
        // ch_versions = ch_versions.mix(SEQTK_SEQ2.out.versions )

    emit:
        contig_fasta = SEQTK_SEQ1.out.fastx
        versions     = ch_versions          
    }
