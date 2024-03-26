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

        // extract meta channel
        meta = fasta.map { meta, fasta -> meta } // [id, contig, size, replicon, single_end, filter_list]
        // create filter_list for subseq to use
        filter_list = fasta
            .map{it[0]}
            .collectFile() {
                item -> [ item.id + ".txt", item.contig + "\n" ]
            }

        SEQTK_SEQ1(fasta)
        ch_versions = ch_versions.mix(SEQTK_SEQ1.out.versions)


        SEQTK_SUBSEQ(SEQTK_SEQ1.out.fastx, filter_list)
        ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions )
        

        SEQTK_SEQ2(SEQTK_SUBSEQ.out.sequences)
        ch_versions = ch_versions.mix(SEQTK_SEQ2.out.versions )

    emit:
        contig_fasta = SEQTK_SEQ2.out.fastx
        versions     = ch_versions          
    }
