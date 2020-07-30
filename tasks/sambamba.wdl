version 1.0

task markdup {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-29"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample
		String suffix = ".markdup"

		Boolean removeDuplicates = false
		Int compressionLevel = 1

		Int threads = 1
		String? tempDir

		Int? hashTableSize
		Int? overflowListSize
		Int? sortBufferSize
		Int? bufferSize
	}

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.sam|\.cram)","")
	String outputBam = if defined(outputPath) then "~{outputPath}/~{sampleName}~{suffix}.bam" else "~{sampleName}~{suffix}.bam"
	String outputBai = if defined(outputPath) then "~{outputPath}/~{sampleName}~{suffix}.bai" else "~{sampleName}~{suffix}.bai"

	command <<<

		if [[ ! -d $(dirname ~{outputBam}) ]]; then
			mkdir -p $(dirname ~{outputBam})
		fi

		~{path_exe} markdup \
			~{true="--remove-duplicates" false="" removeDuplicates} \
			--nthreads ~{threads} \
			--compression-level ~{compressionLevel} \
			~{default="" "--tmpdir " + tempDir} \
			~{default="" "--hash-table-size " + hashTableSize} \
			~{default="" "--overflow-list-size " + overflowListSize} \
			~{default="" "--sort-buffer-size " + sortBufferSize} \
			~{default="" "--io-buffer-size " + bufferSize} \
			~{in} ~{outputBam}

	>>>

	output {
		File outputBam = outputBam
		File outputBai = outputBai
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "sambamba"]',
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
		in: {
			description: 'Bam file to mark or remove duplicates.',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. sample.suffix.bam) [default: ".markdup"]',
			category: 'optional'
		}
		removeDuplicates: {
			description: 'Remove duplicates instead of just marking them [default: false]',
			category: 'optional'
		}
		compressionLevel: {
			description: 'Specify compression level of the resulting file (from 0 to 9) [default: 1]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
		hashTableSize: {
			description: 'Size of hash table for finding read pairs',
			category: 'optional'
		}
		overflowListSize: {
			description: 'Size of the overflow list where reads, thrown from the hash table, get a second chance to meet their pairs',
			category: 'optional'
		}
		sortBufferSize: {
			description: 'Total amount of memory (in *megabytes*) used for sorting purposes',
			category: 'optional'
		}
		bufferSize: {
			description: 'Two buffers of BUFFER_SIZE *megabytes* each are used for reading and writing BAM during the second pass',
			category: 'optional'
		}
	}
}
