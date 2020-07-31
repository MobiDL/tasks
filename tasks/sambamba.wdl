version 1.0

task markdup {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
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
	String outputBai = if defined(outputPath) then "~{outputPath}/~{sampleName}~{suffix}.bam.bai" else "~{sampleName}~{suffix}.bam.bai"

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
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
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

task index {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample

		Boolean checkBins = false

		Int threads = 1
	}

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.cram)","")
	String ext = sub(basename(in),"^.*(bam|cram)","$1")
	Boolean cram = if ext=="cram" then true else false
	String extOut = if cram then ".crai" else "" # fix the fact that sambamba write automatically
	String extFile = if cram then "" else ".bai" # extension for output crai but not for bai
	String outputIdx = if defined(outputPath) then "~{outputPath}/~{sampleName}~{extFile}" else "~{sampleName}~{extFile}"

	command <<<

		if [[ ! -d $(dirname ~{outputIdx}) ]]; then
			mkdir -p $(dirname ~{outputIdx})
		fi

		~{path_exe} index \
			--nthreads ~{threads} \
			~{true="--check-bins" false="" checkBins} \
			~{true="--cram-input" false="" cram} \
			~{in} \
			~{outputIdx}

	>>>

	output {
		File outputIdx = outputIdx + extOut
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
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam or cram file to index.',
			category: 'Required'
		}
		checkBins: {
			description: 'check that bins are set correctly',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task sort {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-31"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample
		String suffix = ".sort"

		String? filter
		Boolean? sortByReadName

		Int compressionLevel = 1
		Boolean uncompressedChuncks = false

		String? memory
		Int threads = 1
		String? tempDir
	}

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.sam|\.cram)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{sampleName}~{suffix}.bam" else "~{sampleName}~{suffix}.bam"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} sort \
			--nthreads ~{threads} \
			~{default="" "--tmpdir " + tempDir} \
			~{default="" "--memory-limit " + memory} \
			~{default="" "--filter " + filter} \
			~{default="" true="--sort-by-name" false="--natural-sort" sortByReadName} \
			--compression-level ~{compressionLevel} \
			~{true="--uncompressed-chunks" false="" uncompressedChuncks} \
			--out ~{outputFile} \
			~{in}

	>>>

	output {
		File outputFile = outputFile
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
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam file to sort.',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. sample.suffix.bam) [default: ".sort"]',
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
		memory: {
			description: 'Sets the max memory (e.g. "2G")',
			category: 'optional'
		}
		tempDir: {
			description: 'Directory for storing intermediate files; default is system directory for temporary files',
			category: 'optional'
		}
		filter: {
			description: 'Keep only reads that satisfy FILTER',
			category: 'optional'
		}
		sortByReadName: {
			description: 'Sort by read name instead of coordinate (true: lexicographical; false: natural) (default: null)',
			category: 'optional'
		}
		uncompressedChuncks: {
			description: 'Write sorted chunks as uncompressed BAM (default is writing with compression level 1), that might be faster in some cases but uses more disk space',
			category: 'optional'
		}
	}
}
