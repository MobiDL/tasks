version 1.0

task sort {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name
		String suffix = ".sort"
		String format = "bam"

		Int compressionLevel = 1
		String? memory
		Boolean sortByReadName = false
		String? tag
		File? refFasta
		Int threads = 1
	}

	Int memTemp = 768*threads
	String totalMem = if defined(memory) then memory else "~{memTemp}M"
	Int mem = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Boolean giga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryByThreads = if giga then floor((mem*1024)/threads) else floor((mem)/threads)

	String baseName = if defined(name) then name else sub(basename(in),"(\.sam|\.bam|\.cram)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}.~{format}" else "~{baseName}~{suffix}.~{format}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} sort \
			-l ~{compressionLevel} \
			-m ~{memoryByThreads}M \
			~{true="-n" false="" sortByReadName} \
			~{default="" "-t " + tag} \
			--output-fmt ~{format} \
			~{default="" "--reference " + refFasta} \
			--threads ~{threads - 1} \
			-o ~{outputFile} \
			~{in}

	>>>

	output {
		File outputFile = outputFile
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "samtools"]',
			category: 'optional'
		}
		outputPath: {
			description: 'Output path where bam file was generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
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
		format: {
			description: 'Specify a single output file format option [default: "bam"]',
			category: 'optional'
		}
		compressionLevel: {
			description: 'Specify compression level of the resulting file (from 0 to 9) [default: 1]',
			category: 'optional'
		}
		sortByReadName: {
			description: 'Sort by read name [default: false]',
			category: 'optional'
		}
		tag: {
			description: 'Sort by value of TAG. Uses position as secondary index (or read name if -n is set)',
			category: 'optional'
		}
		refFasta: {
			description: 'Reference sequence FASTA FILE',
			category: 'optional'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (768M*threads)]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task dict {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name
		String suffix = ".dict"

		String? assembly
		Boolean header = true
		String? species
		String? uri
	}

	String baseName = if defined(name) then name else basename(in)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}" else "~{baseName}~{suffix}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} dict \
			~{default="" "--assembly " + assembly} \
			~{true="" false="--no-header"} \
			~{default="" "--species " + species} \
			~{default="" "--uri " + uri} \
			-o ~{outputFile} \
			~{in}

	>>>

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "samtools"]',
			category: 'optional'
		}
		outputPath: {
			description: 'Output path where dict file was generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: basename(in)]',
			category: 'optional'
		}
		in: {
			description: 'Fasta file.',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. mygenome.fasta.dict) [default: ".dict"]',
			category: 'optional'
		}
		assembly: {
			description: 'Assembly',
			category: 'optional'
		}
		header: {
			description: 'Print the header (@HD line) [default: true]',
			category: 'optional'
		}
		uri: {
			description: 'URI (e.g. file:///abs/path/to/file.fa)',
			category: 'optional'
		}
	}
}
