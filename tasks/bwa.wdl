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
		String? sample
		String subString = "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"

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
		Int memoryByThreads = 768
		String? memory
	}

	Int memTemp = memoryByThreads*threads
	String totalMem = if defined(memory) then memory else "~{memTemp}M"
	Int mem = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Boolean giga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryByThreadsCalc = if giga then floor((mem*1024)/threads) else floor((mem)/threads)

	String baseName = if defined(sample) then sample else sub(basename(fastqR1),subString,"")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}.bam" else "~{baseName}.bam"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe_bwa} mem \
			-R "@RG\tID:~{baseName}\tSM:~{baseName}\tPL:~{platformReads}" \
			-T ~{minScore} \
			~{true="-M" false="" markShorter} \
			-t ~{threads} \
			~{refFasta} \
			~{fastqR1} ~{default="" fastqR2} \
			| ~{path_exe_samtools} sort -@ ~{threads-1} -m ~{memoryByThreadsCalc}M -o ~{outputFile}

	>>>

	output {
		File outputFile = outputFile
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
			description: 'Sample name to use for output file name [default: sub(basename(fastqR1),subString,"")]',
			category: 'optional'
		}
		subString: {
			description: 'Substring to remove to get sample name [default: "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"]',
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

task index {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		String path_exe = "bwa"

		File in
		String? outputPath

		String? algo
		String? name
		Int? blockSize
		Boolean sixtyFour = false

		Int threads = 1
	}

	String baseName = if defined(name) then name else basename(in)
	String outputBaseName = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d $(dirname ~{outputBaseName}) ]]; then
			mkdir -p $(dirname ~{outputBaseName})
		fi

		~{path_exe} index \
			~{default="" "-a " + algo} \
			~{default="" "-b " + blockSize} \
			~{true="-6" false="" sixtyFour} \
			-p ~{outputBaseName} \
			~{in}

	>>>

	output {
		File refAmb = outputBaseName + ".amb"
		File refAnn = outputBaseName + ".ann"
		File refPac = outputBaseName + ".pac"
		File refBwt = outputBaseName + ".bwt"
		File refSa = outputBaseName + ".sa"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bwa"]',
			category: 'optional'
		}
		outputPath: {
			description: 'Output path where bam file was generated. [default: pwd()]',
			category: 'optional'
		}
		in:	{
			description: 'Fasta file to index by bwa.',
			category: 'required'
		}
		algo: {
			description: 'BWT construction algorithm: bwtsw, is or rb2 [auto]',
			category: 'optional'
		}
		name: {
			description: 'Prefix of the index files names [fasta name]',
			category: 'optional'
		}
		blockSize: {
			description: 'Block size for the bwtsw algorithm (effective with algo: bwtsw)',
			category: 'optional'
		}
		sixtyFour: {
			description: 'Index files named as <in.fasta>.64.* instead of <in.fasta>.*',
			category: 'optional'
		}
	}
}
