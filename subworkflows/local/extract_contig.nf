//
// extract contig as fasta format from fasta file
//

include { SEQTK_SEQ1                  } from '../modules/nf-core/seqtk/seq/main'
include { SEQTK_SEQ2                  } from '../modules/nf-core/seqtk/seq/main'
include { SEQTK_SUBSEQ                } from '../modules/nf-core/seqtk/subseq/main'

workflow EXTRACT_CONTIG {
    take:
        tuple val (meta), path (files)

    main:

        // convert contig name to file as required by seqtk subseq
        filter_list = file('filter_list.txt')
        filter_list.text = $meta.contig

        SEQTK_SEQ1(meta, files[1])
        ch_versions = ch_versions.mix(SEQTK_SEQ1.out.versions)
        ch_extracted_fastq = SEQTK_SEQ1.out.fastx[1]

        SEQTK_SUBSEQ(ch_extracted_fastq, filter_list)
        ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions )

        SEQTK_SEQ2(meta, SEQTK_SUBSEQ.out.sequences)
        ch_versions = ch_versions.mix(SEQTK_SEQ2.out.versions )

    emit:
        contig_fasta = SEQTK_SEQ2.out.fastx
        versions     = ch_versions          
    }
