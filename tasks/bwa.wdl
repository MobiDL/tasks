version 1.0

# Options optimizations :
#	- It seems that option "-V" didn't do anything
#	- Since minimap2 replace BWA for long-reads it seems deprecated to use "-x"
#	- It seems unnecessary to change all options in the human context

task mem {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-29"
	}

	input {
		String path_exe_bwa = "bwa"
		String path_exe_samtools = "samtools"

		String? outputPath
		String sample = sub(basename(fastqR1),"(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?","")

		File fastqR1
		File? fastqR2

		File refFasta
		File refFai
		File refAmb
		File refAnn
		File refBwt
		File refPac
		File refSa

		String platformReads = "ILLUMINA"

		Boolean markShorter = true
		Int minScore = 30

		Int threads = 1
	}

	command <<<

		if [[ ! -f ~{outputPath + "/"}~{sample}.bam ]]; then
			mkdir -p $(dirname ~{outputPath + "/"}~{sample}.bam)
		fi

		~{path_exe_bwa} mem \
			-R "@RG\tID:~{sample}\tSM:~{sample}\tPL:~{platformReads}" \
			-T ~{minScore} \
			~{true="-M" false="" markShorter} \
			-t ~{threads} \
			~{refFasta} \
			~{fastqR1} ~{default="" fastqR2} \
			| ~{path_exe_samtools} sort -@ ~{threads-1} -m 768M -o ~{outputPath + "/"}~{sample}.bam

	>>>

	output {
		File outputFile = '~{outputPath + "/"}~{sample}.bam'
	}

	parameter_meta {
		path_exe_bwa: {
			description: 'Path used as executable [default: "bwa"]',
			category: 'optional'
		}
		path_exe_samtools: {
			description: 'Path used as executable [default: "samtools"]',
			category: 'optional'
		}
		outputPath: {
			description: 'Output path where bam file was generated. [default: pwd()]',
			category: 'optional'
		}
		sample: {
			description: 'Sample name to use for output file name [default: sub(basename(fastqR1),"(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?","")]',
			category: 'optional'
		}
		fastqR1: {
			description: 'Input file with reads 1 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		fastqR2: {
			description: 'Input file with reads 2 (fastq, fastq.gz, fq, fq.gz).',
			category: 'optional'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'required'
		}
		platformReads: {
			description: 'Type of plateform that produce reads [default: ILLUMINA]',
			category: 'optional'
		}
		markShorter: {
			description: 'Mark shorter split hits as secondary (for Picard compatibility). [default: true]',
			category: 'optional'
		}
		minScore: {
			description: 'Don’t output alignment with score lower than this value. [default: 30]',
			category: 'optional'
		}
		threads: {
			description: "Sets the number of threads [default: 1]",
			category: "optional"
		}
	}
}
