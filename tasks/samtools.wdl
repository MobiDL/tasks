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
		Boolean sortByReadName = false
		String? tag
		File? refFasta

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	Int memTemp = memoryByThreads*threads
	String totalMem = if defined(memory) then memory else "~{memTemp}M"
	Int mem = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Boolean giga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memByThreads = if giga then floor((mem*1024)/threads) else floor((mem)/threads)

	String baseName = if defined(name) then name else sub(basename(in),"(\.sam|\.bam|\.cram)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}.~{format}" else "~{baseName}~{suffix}.~{format}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} sort \
			-l ~{compressionLevel} \
			~{true="-n" false="" sortByReadName} \
			~{default="" "-t " + tag} \
			--output-fmt ~{format} \
			~{default="" "--reference " + refFasta} \
			--threads ~{threads - 1} \
			-m ~{memByThreads}M \
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
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'optional'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
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
		String ext = ".dict"
		String subString = "\.fa(sta)?(\.gz)?"

		String? assembly
		Boolean header = true
		String? species
		String? uri
	}

	String baseName = if defined(name) then name else sub(basename(in),subString,"")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{ext}" else "~{baseName}~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} dict \
			~{default="" "--assembly " + assembly} \
			~{true="" false="--no-header" header} \
			~{default="" "--species " + species} \
			~{default="" "--uri " + uri} \
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
			description: 'Output path where dict file was generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),subString,"")]',
			category: 'optional'
		}
		in: {
			description: 'Fasta file.',
			category: 'Required'
		}
		ext: {
			description: 'Extension of the output file (e.g. mygenome.dict) [default: ".dict"]',
			category: 'optional'
		}
		subString: {
			description: 'Substring to remove to get file basename [default: "\.fa(sta)?(\.gz)?"]',
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

task index {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-31"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name

		Int? minIntervalSize
		Boolean csi = false
		Int threads = 1
	}

	String extFile = sub(basename(in),"(.*)\.(bam|cram)$","$2")
	Boolean csiOpt = if (defined(minIntervalSize) || csi) then true else false
	String ext = if extFile=="cram" then ".crai" else if csiOpt then ".csi" else ".bai"

	String baseName = if defined(name) then name else basename(in)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{ext}" else "~{baseName}~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} index \
			~{true="-c" false="-b" csiOpt} \
			~{default="" "-m " + minIntervalSize} \
			-@ ~{threads} \
			~{in} \
			~{outputFile}

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
			description: 'Output path where index was generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: basename(in)]',
			category: 'optional'
		}
		in: {
			description: 'Input bam or cram to index.',
			category: 'Required'
		}
		minIntervalSize: {
			description: 'Set minimum interval size for CSI indices to 2^INT.',
			category: 'optional'
		}
		csi: {
			description: 'Generate CSI-format index for BAM files (not functionnal with cram files) [default: false]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task view {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-31"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name
		Boolean cram = false

		File? faidx
		File? bed
		String? readGroup
		File? readGroupFile
		Int minQuality = 0
		String? library
		Int minCigarOp = 0

		Boolean multiRegioOperator = false
		Boolean collapseCigarOp =false

		File? refFasta

		Int threads = 1
	}

	String ext = if cram then "cram" else "bam"
	String baseName = if defined(name) then name else sub(basename(in),"(\.bam|\.cram|\.sam)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}.~{ext}" else "~{baseName}.~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} view \
			~{true="-C" false="-b" cram} \
			-o ~{outputFile} \
			~{default="" "-t " + faidx} \
			~{default="" "-L " + bed} \
			~{default="" "-r " + readGroup} \
			~{default="" "-R " + readGroupFile} \
			-q ~{minQuality} \
			~{default="" "-l " + library} \
			-m ~{minCigarOp} \
			~{true="-M" false="" multiRegioOperator} \
			~{true="-B" false="" collapseCigarOp} \
			--output-fmt ~{ext} \
			~{default="" "--reference " + refFasta} \
			--threads ~{threads - 1} \
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
			description: 'Output path where index was generated. [default: pwd()]',
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
		cram: {
			description: 'Output to cram format. (bam otherwise) [default: false]',
			category: 'optional'
		}
		faidx: {
			description: 'Listing reference names and lengths (faidx)',
			category: 'optional'
		}
		bed: {
			description: 'Only include reads overlapping this BED FILE',
			category: 'optional'
		}
		readGroup: {
			description: 'Only include reads in read group.',
			category: 'optional'
		}
		readGroupFile: {
			description: 'Only include reads with read group listed in a file',
			category: 'optional'
		}
		minQuality: {
			description: 'Only include reads with mapping quality >= value [default: 0]',
			category: 'optional'
		}
		library: {
			description: 'Only include reads in library.',
			category: 'optional'
		}
		minCigarOp: {
			description: 'Only include reads with number of CIGAR operations consuming query sequence >= value [default: 0]',
			category: 'optional'
		}
		multiRegioOperator: {
			description: 'Use the multi-region iterator (increases the speed, removes duplicates and outputs the reads as they are ordered in the file) [default: false]',
			category: 'optional'
		}
		collapseCigarOp: {
			description: 'Collapse the backward CIGAR operation [default: false]',
			category: 'optional'
		}
		refFasta: {
			description: 'Reference sequence FASTA FILE [null]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task faidx {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-05"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name

		Int threads = 1
	}

	String baseName = if defined(name) then name else basename(in)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}.fai" else "~{baseName}.fai"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} faidx ~{in}
		mv ~{in + ".fai"} ~{outputFile}

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
			description: 'Output path where index was generated. [default: pwd()]',
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
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task fqidx {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-07"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name

		Int threads = 1
	}

	String baseName = if defined(name) then name else basename(in)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}.fai" else "~{baseName}.fai"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} fqidx ~{in}
		mv ~{in + ".fai"} ~{outputFile}

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
			description: 'Output path where index was generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: basename(in)]',
			category: 'optional'
		}
		in: {
			description: 'Fastq file.',
			category: 'Required'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task flagstat {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-31"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? sample
		String ext = ".flagstats"

		Int threads = 1
	}

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.sam|\.cram)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{sampleName}~{ext}" else "~{sampleName}~{ext}"


	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} flagstat \
			--threads ~{threads - 1} \
			~{in} > ~{outputFile}

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
			description: 'Output path where flagstat file was generated. [default: pwd()]',
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
		ext: {
			description: 'Extension of the output file [default: ".flagstats"]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task bedcov {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		String path_exe = "samtools"

		File inBam
		File? inBamIdx
		File inBed
		String? outputPath
		String? name
		String suffix = ".bedcov"

		Int qualityThreshold = 0
		Boolean includeDel = true

		Int threads = 1
	}

	String baseName = if defined(name) then name else basename(inBam)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}.bed" else "~{baseName}~{suffix}.bed"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} bedcov \
			-Q ~{qualityThreshold} \
			~{true="" false="-j" includeDel} \
			~{inBed} \
			~{inBam} > ~{outputFile}

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
			description: 'Output path where bedcov file will be generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: basename(in)]',
			category: 'optional'
		}
		inBed: {
			description: 'Bed file.',
			category: 'Required'
		}
		inBam: {
			description: 'Bam file.',
			category: 'Required'
		}
		inBamIdx: {
			description: 'Index of the input bam file.',
			category: 'optional'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. sample.bedcov.bed) [default: ".bedcov"]',
			category: 'optional'
		}
		qualityThreshold: {
			description: 'Mapping quality threshold [default: 0]',
			category: 'optional'
		}
		includeDel: {
			description: 'Include deletions (D) and ref skips (N) in bedcov computation [default: true]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task markdup {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-11"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name
		String? outExt
		String? outputOpt
		String suffix = ".markdup"

		Boolean remove = false
		Int maxReadsLength = 300
		Boolean markSupAl = false
		Boolean reportStat = false
		Boolean markPrimDup = false

		Int threads = 1
	}

	String ext = if defined(outExt) then outExt else sub(basename(in),"(.*)\.(sam|bam|cram)$","$2")
	String baseName = if defined(name) then name else sub(basename(in),"(.*)\.(sam|bam|cram)$","$1")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}.~{ext}" else "~{baseName}~{suffix}.~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} markdup \
			~{true="-r" false="" remove} \
			-l ~{maxReadsLength} \
			~{true="-S" false="" markSupAl} \
			~{true="-s" false="" reportStat} \
			~{true="-t" false="" markPrimDup} \
			--output-fmt ~{ext}~{default="" "," + outputOpt} \
			--threads ~{threads - 1} \
			~{in} \
			~{outputFile}

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
			description: 'Output path where bedcov file will be generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: basename(in)]',
			category: 'optional'
		}
		outExt: {
			description: 'Specify output format (SAM, BAM, CRAM) [default: same as input]',
			category: 'optional'
		}
		in: {
			description: 'Bam file.',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. sample.markdup.bam) [default: ".markdup"]',
			category: 'optional'
		}
		outputOpt: {
			description: 'Specify output file format option in the form',
			category: 'optional'
		}
		remove: {
			description: 'Remove duplicates instead of just marking them [default: false]',
			category: 'optional'
		}
		maxReadsLength: {
			description: 'Max read length [default: 300]',
			category: 'optional'
		}
		markSupAl: {
			description: 'Mark supplemenary alignments of duplicates as duplicates (slower) [default: false]',
			category: 'optional'
		}
		reportStat: {
			description: 'Report stats. [default: false]',
			category: 'optional'
		}
		markPrimDup: {
			description: 'Mark primary duplicates with the name of the original in a "do" tag. [default: false]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}

task fixmate {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-11"
	}

	input {
		String path_exe = "samtools"

		File in
		String? outputPath
		String? name
		String? outExt
		String? outputOpt
		String suffix = ".fixmate"

		Boolean remove = false
		Boolean disableFR = false
		Boolean addTemplateCigar = false
		Boolean addMateScoreTag = false

		Int threads = 1
	}

	String ext = if defined(outExt) then outExt else sub(basename(in),"(.*)\.(sam|bam|cram)$","$2")
	String baseName = if defined(name) then name else sub(basename(in),"(.*)\.(sam|bam|cram)$","$1")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}.~{ext}" else "~{baseName}~{suffix}.~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} fixmate \
			~{true="-r" false="" remove} \
			~{true="-p" false="" disableFR} \
			~{true="-c" false="" addTemplateCigar} \
			~{true="-m" false="" addMateScoreTag} \
			--output-fmt ~{ext}~{default="" "," + outputOpt} \
			--threads ~{threads - 1} \
			~{in} \
			~{outputFile}

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
			description: 'Output path where bedcov file will be generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output file name [default: basename(in)]',
			category: 'optional'
		}
		outExt: {
			description: 'Specify output format (SAM, BAM, CRAM) [default: same as input]',
			category: 'optional'
		}
		in: {
			description: 'Bam file.',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. sample.markdup.bam) [default: ".fixmate"]',
			category: 'optional'
		}
		outputOpt: {
			description: 'Specify output file format option in the form',
			category: 'optional'
		}
		remove: {
			description: 'Remove unmapped reads and secondary alignments [default: false]',
			category: 'optional'
		}
		disableFR: {
			description: 'Disable FR proper pair check [default: false]',
			category: 'optional'
		}
		addTemplateCigar: {
			description: 'Add template cigar ct tag [default: false]',
			category: 'optional'
		}
		addMateScoreTag: {
			description: 'Add mate score tag [default: false]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}
