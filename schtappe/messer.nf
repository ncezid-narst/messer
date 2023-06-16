#!/usr/bin/env nextflow

/*
Nextflow pipeline to recover and reassemble small plasmids with ONT data
*/

/*
PROCESSES
*/
process EXTRACT {
	tag "Extracting contig for ${replicon}"
	publishDir (
		path: "${params.outdir}/${sample_id}/${replicon}/",
		mode: 'copy',
		saveAs: { filename -> 
			if (filename.matches "${replicon}.fasta") "$filename"
			else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(assembly), path(reads)
	
	output:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(assembly), path(reads), path("${replicon}.fasta")
	
	script:
	"""
	touch name.lst
	printf ${contig_name} >> name.lst
	seqtk seq -F '#' ${assembly} > "${sample_id}.fastq"
	seqtk subseq "${sample_id}.fastq" name.lst > "${replicon}.fastq"
	seqtk seq -a "${replicon}.fastq" > "${replicon}.fasta"
	"""
}

process ALIGN {
	tag "Aligning reads to ${replicon} contig"
	publishDir (
		path: "${params.outdir}/${sample_id}/${replicon}/",
		mode: 'copy',
		saveAs: { filename -> 
			if (filename.matches 'aln.sam') "$filename"
			else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(assembly), path(reads), path(replicontig)
	
	output:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(assembly), path(reads), path(replicontig), path("aln.sam")
	
	script:
	"""
	minimap2 -ax map-ont --sam-hit-only ${replicontig} ${reads} > "aln.sam"
	"""
}

process READSGET {
	tag "Formatting reads for ${replicon}"
	publishDir (
		path: "${params.outdir}/${sample_id}/${replicon}/reads/",
		mode: 'copy',
		saveAs: { filename -> 
			if (filename.matches "mapped.fastq") "$filename"
			else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(assembly), path(reads), path(replicontig), path(samfile)
	
	output:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(assembly), path(reads), path(replicontig), path(samfile), path("mapped.fastq")
	
	script:
	"""
	samtools fastq ${samfile} > "mapped.fastq"
	"""
}

process READFILTERING {
	tag "Nanoq on ${replicon}"
	publishDir(
		path: "${params.outdir}/${sample_id}/${replicon}/reads/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith('.fastq.gz')) "$filename"
					else if (filename.endsWith('.log')) "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(assembly), path(reads), path(replicontig), path(samfile), path(mapped_reads)
	
	output:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path("${replicon}_nanoq.fastq.gz"), path("nanoq.log")
	
	script:
	"""
	nanoq -i $mapped_reads -l ${params.nanoq_min_length} -m ${params.nanoq_max_length} -o ${replicon}_nanoq.fastq.gz -vvv -H >> nanoq.log 2>&1
	"""
}

process DOWNSAMPLE {
	tag "Rasusa on $replicon"
	publishDir(
		path: "${params.outdir}/${sample_id}/${replicon}/reads/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith('.fastq.gz')) "$filename"
					else if (filename.endsWith('.log')) "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(reads)
	
	output:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path("${replicon}_nanoq_rasusa.fastq.gz"), path("rasusa.log")
	
	script:
	"""
	rasusa -g ${contig_size} -c ${params.rasusa_coverage} -i ${reads} -o ${replicon}_nanoq_rasusa.fastq.gz >> rasusa.log 2>&1
	"""
}

process ASSEMBLE {
	tag "Flye on ${replicon}"
	publishDir (
		path: "${params.outdir}/${sample_id}/${replicon}/",
		mode: 'copy',
		saveAs: { filename -> 
			if (filename.matches 'flye') "$filename"
			else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path(reads)
	
	output:
	tuple val(sample_id), val(contig_name), val(contig_size), val(replicon), path("flye")
	
	script:
	"""
	flye ${params.flye_read_type} ${reads} -g ${contig_size} -o flye --threads ${params.flye_threads} --min-overlap ${params.flye_minoverlap}
	"""
}

/*
WORKFLOW
*/
workflow {
	Channel //set channel for intermediate/template assembly
		.fromPath(params.assembly, checkIfExists: true)
		.set{assembly_ch}
	Channel //set channel for reads
		.fromPath(params.reads, checkIfExists: true)
		//.view()
		.set{reads_ch}
	Channel //set channel for tsv to collect sample ID, contig name, correct contig size
		.fromPath(params.contiginfo, checkIfExists: true)
		.splitCsv(sep: '\t', header:true, strip:true)
		.map{ row -> tuple(row.WGSID, row.CONTIG, row.SIZE, row.REPLICON) }
		.map { tuple -> 
		def id = tuple[0]
		def contig = tuple[1]
		def size = tuple[2]
		def replicon = tuple[3]
		
		[id, contig, size, replicon]
		}
		.combine(assembly_ch)
		.combine(reads_ch)
		//.view()
		.set{combined_ch}
	
	replicon_ch = EXTRACT(combined_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6] }
		//.view()
	alignment_ch = ALIGN(replicon_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7] }
		//.view()
	readsget_ch = READSGET(alignment_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7,8] }
		//.view()
	nanoq_ch = READFILTERING(readsget_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4] }
		//.view()
	rasusa_ch = DOWNSAMPLE(nanoq_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4] }
		//.view()
	plasmidassembly_ch = ASSEMBLE(rasusa_ch)
}
