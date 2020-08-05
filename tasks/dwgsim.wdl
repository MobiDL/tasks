version 1.0

task simulateReadsIllumina {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-05"
	}

	input {
		String path_exe = "dwgsim"

		File in
		String? outputPath
		String name = "dwgsim-illumina"

		Int sizeR1 = 150
		Int sizeR2 = sizeR1

		Float rateMut = 0.0010
		Float freqMut = 0.5000
		Float rateIndels = 0.1000
		Float probExt = 0.3000

		File? target
		String? readPrefix

		Int threads = 1
	}

	String outputPrefix = if defined(outputPath) then "~{outputPath}/~{name}" else "~{name}"

	command <<<

		if [[ ! -d $(dirname ~{outputPrefix}) ]]; then
			mkdir -p $(dirname ~{outputPrefix})
		fi


		~{path_exe} \
			-1 ~{sizeR1} \
			-2 ~{sizeR2} \
			-r ~{rateMut} \
			-F ~{freqMut} \
			-R ~{rateIndels} \
			-X ~{probExt} \
			~{default="" "-x " + target} \
			~{default="" "-P " + readPrefix} \
			~{in} \
			~{outputPrefix}
	>>>

	output {
		File bfastFQ = outputPrefix + ".bfast.fastq"
		File fastqR1 = outputPrefix + ".bwa.read1.fastq"
		File fastqR2 = outputPrefix + ".bwa.read2.fastq"
		File mutationsTXT = outputPrefix + ".mutations.txt"
		File mutationsVCF = outputPrefix + ".mutations.vcf"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "dwgsim"]',
			category: 'optional'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: dwgsim-illumina]',
			category: 'optional'
		}
		in: {
			description: 'Fasta file used to create fastq.',
			category: 'Required'
		}
		sizeR1: {
			description: 'Length of the first read [default: 150]',
			category: 'optional'
		}
		sizeR2: {
			description: 'Length of the second read [default: sizeR1]',
			category: 'optional'
		}
		rateMut: {
			description: 'Rate of mutations [default: 0.0010]',
			category: 'optional'
		}
		freqMut: {
			description: 'Frequency of given mutation to simulate low fequency somatic mutations [default: 0.5000]',
			category: 'optional'
		}
		rateIndels: {
			description: 'Fraction of mutations that are indels [default: 0.1000]',
			category: 'optional'
		}
		probExt: {
			description: 'Probability an indel is extended [default: 0.3000]',
			category: 'optional'
		}
		target: {
			description: 'Bed of regions to cover [default: WGS]',
			category: 'optional'
		}
		readPrefix: {
			description: 'Read prefix to prepend to each read name',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}
