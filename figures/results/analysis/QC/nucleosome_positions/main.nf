!#/usr/bin/env nextflow

nextflow.enable.dsl=2








process splitReads {

	publishDir "/proj/sens2021596/nobackup/results/results/analysis/QC/nucleosome_positions/", mode: 'copy', overwrite: true
	cpus 4
	time '2h'
	clusterOptions '-A sens2021596'
	module 'bioinfo-tools:R_packages/4.0.0'

	input:
	tuple val()

	output:


	script:
	"""
	!#/usr/bin/env Rscript

	library(ATACseqQC)
	library(BSgenome.Hsapiens.UCSC.hg38)
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)

	txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
	genome <- Hsapiens

	objs <- splitGAlignmentsByCut()
	"""
}
