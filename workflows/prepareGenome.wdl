version 1.0

import "../tasks/bwa.wdl" as bwa
import "../tasks/bash.wdl" as bash
import "../tasks/samtools.wdl" as samtools

workflow prepareGenome {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
	}

	input {
		File fasta
		String outputPath
		Int threads = 1
	}

	call bash.makeLink as GenomeLn {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = threads
	}

	call bwa.index as BwaIndexGenome {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = threads
	}

	call samtools.faidx as SamtoolsIndexGenome {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = threads
	}

	call samtools.dict as SamtoolsDictGenome {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = threads
	}

	output {
		File refFasta = GenomeLn.outputFile
		File refFai = SamtoolsIndexGenome.outputFile
		File refDict = SamtoolsDictGenome.outputFile
		File refAmb = BwaIndexGenome.refAmb
		File refAnn = BwaIndexGenome.refAnn
		File refBwt = BwaIndexGenome.refBwt
		File refPac = BwaIndexGenome.refPac
		File refSa = BwaIndexGenome.refSa
	}

	parameter_meta {
		fasta: {
			description : 'Path to the reference file (format: fasta)',
			category : 'Required'
		}
		outputPath: {
			description : 'Output path where bam file was generated.',
			category : 'Required'
		}
		threads: {
			description : 'Sets the number of threads [default: 1]',
			category : 'System'
		}
	}
}
