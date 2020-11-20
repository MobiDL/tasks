version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2020 MoBiDiC
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

task get_version {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-11-20"
	}

	input {
		String path_exe = "bcftools"

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	command <<<
		~{path_exe} --version
	>>>

	output {
		String output = read_string(stdout())
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bcftools"]',
			category: 'System'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}

task index {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-24"
	}

	input {
		String path_exe = "bcftools"

		File in
		String? outputPath

		Boolean tabix = false
		Boolean force = false
		Int minShift = 14

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String ext = if tabix then ".tbi" else ".csi"
	Boolean outputDef = defined(outputPath)
	String outputRaw = if outputDef then outputPath + "/" + basename(in) else in
	String outputIndex = outputRaw + ext

	String minShiftOpt = if ! tabix then "--min-shift " + minShift else ""

	command <<<

		if [[ ! -f ~{outputRaw} ]]; then
			mkdir -p $(dirname ~{outputRaw})
			ln ~{in} ~{outputRaw}
		fi

		~{path_exe} index \
			~{true="--tbi" false="--csi" tabix} \
			~{true="--force" false="" force} \
			~{minShiftOpt} \
			--threads ~{threads - 1} \
			-o ~{outputIndex} \
			~{outputRaw}

	>>>

	output {
		File vcf = "~{outputRaw}"
		File index = "~{outputIndex}"
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bcftools"]',
			category: 'System'
		}
		in: {
			description: "VCF/BCF file to index (extension: '.vcf.gz|.bcf')",
			category: 'Required'
		}
		outputPath: {
			description: "Path where was generated output [default: VCF path]",
			category: 'Output path/name option'
		}
		tabix: {
			description: "Generate TBI-format index for VCF files [default: CSI-format]",
			category: 'Tool option'
		}
		force: {
			description: "Overwrite index if it already exists [default: false]",
			category: 'Tool option'
		}
		minShift: {
			description: "Set minimal interval size for CSI indices to 2^INT (skipped if tabix used) [default: 14]",
			category: 'Tool option'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}

# TODO:
# 	- simplify options :
#		- remove opt : file-list (in worklow read file to convert into file list or reverse)
#		- correction for options -i default value
#		- need to test gvcf option more accuratly
task merge {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-24"
	}

	input {
		String path_exe = "bcftools"

		Array[File] in
		Array[File] index
		String? outputPath
		String prefix = "merge"

		Boolean forceSamples = false
		Boolean printHeader = false
		File? useHeader
		Boolean missingToRef = false
		Array[String]? applyFilter
		String filterLogic = "+"
		File? referenceFasta
		Boolean gvcf = false
		String infoRules = "-"
		File? fileList
		String merge = "both"
		Boolean noVersion = false
		String outputType = "v"
		String? regions
		File? regionsFile

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String useHeaderOpt = if defined(useHeader) then "--use-header ~{useHeader} " else ""
	String applyFilterOpt = if defined(applyFilter) then "--apply-filter '~{sep=',' applyFilter}' " else ""
	String gvcfOpt = if gvcf then if defined(referenceFasta) then "--gvcf ~{referenceFasta} " else "--gvcf - " else ""
	String fileListOpt = if defined(fileList) then "--file-list ~{fileList} " else ""
	String regionsOpt = if defined(regions) then "--regions ~{regions} " else ""
	String regionsFileOpt = if defined(regionsFile) then "--regions-file ~{regionsFile} " else ""

	String ext = if outputType == "v" then ".vcf" else if outputType == "b" then ".bcf.gz" else if outputType == "u" then ".bcf" else ".vcf.gz"

	String outputFile = if defined(outputPath) then "~{outputPath}/~{prefix}~{ext}" else "~{prefix}~{ext}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} merge \
			~{true="--force-samples" false="" forceSamples} \
			~{true="--print-header" false="" printHeader} \
			~{useHeaderOpt} \
			~{true="--missing-to-ref" false="" missingToRef} \
			~{applyFilterOpt} \
			~{gvcfOpt} \
			~{fileListOpt} \
			~{true="--no-version" false="" noVersion} \
			~{regionsOpt} \
			~{regionsFileOpt} \
			--threads ~{threads - 1} \
			--info-rules '~{infoRules}' \
			--merge '~{merge}' \
			--filter-logic '~{filterLogic}' \
			-O '~{outputType}' \
			-o ~{outputFile} \
			~{sep=" \\\n\t" in}

	>>>

	output {
		File outputMerge = "~{outputFile}"
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bcftools"]',
			category: 'System'
		}
		in: {
			description: 'Array of files to merge together (extension: ".vcf.gz|.bcf")',
			category: 'Required'
		}
		outputPath: {
			description: 'Path where was generated output [default: VCF path]',
			category: 'Output path/name option'
		}
		forceSamples: {
			description: 'Resolve duplicate sample names [default: false]',
			category: 'Tool option'
		}
		printHeader: {
			description: 'Print only the merged header and exit [default: false]',
			category: 'Tool option'
		}
		useHeader: {
			description: 'Use the provided header [default: -]',
			category: 'Tool option'
		}
		missingToRef: {
			description: 'Assume genotypes at missing sites are 0/0 [default: false]',
			category: 'Tool option'
		}
		applyFilter: {
			description: 'Require at least one of the listed FILTER strings (e.g. "PASS,.")',
			category: 'Tool option'
		}
		filterLogic: {
			description: 'Remove filters if some input is PASS ("x"), or apply all filters ("+") [default: "+"]',
			category: 'Tool option'
		}
		gvcf: {
			description: 'Merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,I16:sum,IDV:max,IMF:max (#doubt) [default: false]',
			category: 'Tool option'
		}
		referenceFasta: {
			description: 'Reference used to merge in gvcf mode',
			category: 'Tool option'
		}
		infoRules: {
			description: 'Rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [-]',
			category: 'Tool option'
		}
		fileList: {
			description: 'Read file names from the file',
			category: 'Tool option'
		}
		merge: {
			description: 'Allow multiallelic records for <snps|indels|both|all|none|id> [default: "both"]',
			category: 'Tool option'
		}
		noVersion: {
			description: 'Do not append version and command line to the header [default: false]',
			category: 'Tool option'
		}
		outputType: {
			description: '"b" compressed BCF; "u" uncompressed BCF; "z" compressed VCF; "v" uncompressed VCF [default: "v"]',
			category: 'Tool option'
		}
		regions: {
			description: "Restrict to comma-separated list of regions.",
			category: 'Tool option'
		}
		regionsFile: {
			description: "Restrict to regions listed in a file.",
			category: 'Tool option'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}

task norm {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.2"
		date: "2020-08-27"
	}

	input {
		String path_exe = "bcftools"

		File in
		String outputPath
		String? name
		String subString = "\.(vcf|bcf)(\.gz)?$"
		String subStringReplace = ""
		String suffix = ".norm"

		File refFasta

		String? checkRef

		Boolean removeDuplicates = false
		String? rmDupType

		Boolean splitMA = false
		String multiallelicType = "both"

		Boolean version = true

		Boolean normalize = true

		String? regions
		File? regionsFile

		String? targets
		File? targetsFile

		Boolean strictFilter = false
		String outputType = "v"

		Int siteWin = 1000

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	Map[String,String] extType = {"v" : ".vcf","u" : ".bcf","z" : ".vcf.gz","b" : ".bcf.gz"}

	String ext = extType[outputType]
	String baseName = if defined(name) then name else sub(basename(in),subString,subStringReplace)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}~{ext}" else "~{baseName}~{suffix}~{ext}"

	String multiallelics = if splitMA then "-~{multiallelicType}" else "+~{multiallelicType}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} norm \
			~{default="" "--check-ref " + checkRef} \
			~{true="--remove-duplicates" false="" removeDuplicates} \
			~{default="" "--rm-dup " + rmDupType} \
			--fasta-ref ~{refFasta} \
			--multiallelics ~{multiallelics} \
			~{true="" false="--no-version" version} \
			~{true="" false="--do-not-normalize" normalize} \
			~{default="" "--regions " + regions} \
			~{default="" "--regions-file " + regionsFile} \
			~{default="" "--targets " + targets} \
			~{default="" "--targets-file " + targetsFile} \
			~{true="--strict-filter" false="" strictFilter} \
			--output-type ~{outputType} \
			--output ~{outputFile} \
			--threads ~{threads - 1} \
			--site-win ~{siteWin} \
			~{in}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bcftools"]',
			category: 'System'
		}
		in: {
			description: "VCF/BCF file to left-align and normalize indels (extension: '.vcf.gz|.bcf')",
			category: 'Required'
		}
		outputPath: {
			description: 'Path where was generated output',
			category: 'Output path/name option'
		}
		name: {
			description: 'Output file base name [default: sub(basename(in),subString,"")].',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Extension to remove from the input file [default: "\.(vcf|bcf)(\.gz)?$"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		suffix: {
			description: 'Suffix to add for the output file (e.g sample.suffix.bam)[default: ".merge"]',
			category: 'Output path/name option'
		}
		refFasta: {
			description: 'Reference used to merge in gvcf mode',
			category: 'Required'
		}
		checkRef: {
			description: 'Check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [default: e]',
			category: 'Tool option'
		}
		removeDuplicates: {
			description: 'Remove duplicate lines of the same type.',
			category: 'Tool option'
		}
		rmDupType: {
			description: 'Remove duplicate snps|indels|both|all|none (implies removeDuplicates)',
			category: 'Tool option'
		}
		splitMA: {
			description: "Split (true) or join (false) multiallelics sites [default= false]",
			category: 'Tool option'
		}
		multiallelicType: {
			description: "Type of Multiallelics to treat for split/join (type: snps|indels|both|any) [default: both]",
			category: 'Tool option'
		}
		version: {
			description: 'Append version and command line to the header [default: true]',
			category: 'Tool option'
		}
		normalize: {
			description: 'Normalize indels (with -m or -c s) [default: true]',
			category: 'Tool option'
		}
		outputType: {
			description: '"b" compressed BCF; "u" uncompressed BCF; "z" compressed VCF; "v" uncompressed VCF [default: "v"]',
			category: 'Tool option'
		}
		regions: {
			description: "Restrict to comma-separated list of regions",
			category: 'Tool option'
		}
		regionsFile: {
			description: "Restrict to regions listed in a file",
			category: 'Tool option'
		}
		targets: {
			description: "Similar to 'regions' but streams rather than index-jumps",
			category: 'Tool option'
		}
		targetsFile: {
			description: "Similar to 'regionsFile' but streams rather than index-jumps",
			category: 'Tool option'
		}
		strictFilter: {
			description: "When merging (-m+), merged site is PASS only if all sites being merged PASS [default: false]",
			category: 'Tool option'
		}
		siteWin: {
			description: "Buffer for sorting lines which changed position during realignment [default: 1000]",
			category: 'Tool option'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}

task stats {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.2"
		date: "2020-08-27"
	}

	input {
		String path_exe = "bcftools"

		File in
		String? outputPath
		String? name
		String subString = "\.(vcf|bcf)(\.gz)?$"
		String subStringReplace = ""
		String ext = ".stats"

		File refFasta
		File refFai

		Array[Float]? afBins
		String? afTag

		Boolean firstAllele = false

		String collapse = "none"

		Int depthMin = 0
		Int depthMax = 500
		Int depthBin = 1

		File? exons

		Array[String]? applyFilters

		String? exclude
		String? include

		Boolean splitByID = false

		Array[String] samples = ["\"-\""]

		Array[String]? userTsTv

		String? regions
		File? regionsFile

		String? targets
		File? targetsFile

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(in),subString,subStringReplace)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{ext}" else "~{baseName}~{ext}"

	Boolean BoolAFBins = if defined(afBins) then true else false
	Boolean BoolApplyFilters = if defined(applyFilters) then true else false
	Boolean BoolUserTsTv = if defined(userTsTv) then true else false

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} stats \
			~{true="--af-bins " false="" BoolAFBins}~{default="" sep="," afBins} \
			~{default="" "--af-tag " + afTag} \
			~{true="--1st-allele-only" false ="" firstAllele} \
			--collapse ~{collapse} \
			--depth ~{depthMin},~{depthMax},~{depthBin} \
			~{default="" "--exons " + exons} \
			~{true="--apply-filters " false="" BoolApplyFilters}~{default="" sep="," applyFilters} \
			~{default="" "--exclude " + exclude} \
			~{default="" "--include " + include} \
			~{true="--split-by-ID" false = "" splitByID} \
			--samples ~{sep="," samples} \
			~{true="--user-tstv " false="" BoolUserTsTv}~{default="" sep="," userTsTv} \
			--fasta-ref ~{refFasta} \
			~{default="" "--regions " + regions} \
			~{default="" "--regions-file " + regionsFile} \
			~{default="" "--targets " + targets} \
			~{default="" "--targets-file " + targetsFile} \
			--threads ~{threads - 1} \
			~{in} > ~{outputFile}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bcftools"]',
			category: 'System'
		}
		in: {
			description: "VCF/BCF file to left-align and normalize indels (extension: '.vcf.gz|.bcf')",
			category: 'Required'
		}
		outputPath: {
			description: 'Path where was generated output',
			category: 'Output path/name option'
		}
		name: {
			description: 'Output file base name [default: sub(basename(in),subString,"")].',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Extension to remove from the input file [default: "\.(vcf|bcf)(\.gz)?$"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		ext: {
			description: 'Extension use for the output file [default: ".stats"]',
			category: 'Output path/name option'
		}
		refFasta: {
			description: 'Reference in fasta format',
			category: 'Required'
		}
		refFai: {
			description: 'Path to the reference file index (format: fai)',
			category: 'Required'
		}
		afBins: {
			description: 'Allele frequency bins',
			category: 'Tool option'
		}
		afTag: {
			description: 'Allele frequency tag to use',
			category: 'Tool option'
		}
		firstAllele: {
			description: 'Include only 1st allele at multiallelic sites',
			category: 'Tool option'
		}
		collapse: {
			description: 'Treat as identical records with <snps|indels|both|all|some|none> [default: none]',
			category: 'Tool option'
		}
		depthMin: {
			description: 'Minimum depth [default: 0]',
			category: 'Tool option'
		}
		depthMax: {
			description: 'Maximum depth [default: 500]',
			category: 'Tool option'
		}
		depthBin: {
			description: 'Bin size depth [default: 1]',
			category: 'Tool option'
		}
		exons: {
			description: 'Tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based, inclusive, bgzip compressed)',
			category: 'Tool option'
		}
		applyFilters: {
			description: 'Require at least one of the listed FILTER strings',
			category: 'Tool option'
		}
		exclude: {
			description: 'Exclude sites for which the expression is true',
			category: 'Tool option'
		}
		include: {
			description: 'Select sites for which the expression is true',
			category: 'Tool option'
		}
		splitByID: {
			description: 'Collect stats for sites with ID separately (known vs novel)',
			category: 'Tool option'
		}
		samples: {
			description: 'List of samples for sample stats ("-" -> all samples)[default: ["-"]]',
			category: 'Tool option'
		}
		userTsTv: {
			description: 'Collect Ts/Tv stats for any tag using the given binning',
			category: 'Tool option'
		}
		regions: {
			description: "Restrict to comma-separated list of regions",
			category: 'Tool option'
		}
		regionsFile: {
			description: "Restrict to regions listed in a file",
			category: 'Tool option'
		}
		targets: {
			description: "Similar to 'regions' but streams rather than index-jumps",
			category: 'Tool option'
		}
		targetsFile: {
			description: "Similar to 'regionsFile' but streams rather than index-jumps",
			category: 'Tool option'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}
