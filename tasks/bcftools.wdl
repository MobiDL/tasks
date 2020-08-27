version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2020  MoBiDiC
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
		Int minShift  = 14

		Int threads   = 1
	}

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

    parameter_meta {
        path_exe: {
			description: "Path used as executable [default: 'bcftools']",
			category: "optional"
		}
		in: {
			description: "VCF/BCF file to index (extension: '.vcf.gz|.bcf')",
			category: "required"
		}
		outputPath: {
			description: "Path where was generated output [default: VCF path]",
			category: "optional"
		}
		tabix: {
			description: "Generate TBI-format index for VCF files [default: CSI-format]",
			category: "optional"
		}
		force: {
			description: "Overwrite index if it already exists [default: false]",
			category: "optional"
		}
		minShift: {
			description: "Set minimal interval size for CSI indices to 2^INT (skipped if tabix used) [default: 14]",
			category: "optional"
		}
		threads: {
			description: "Sets the number of threads [default: 1]",
			category: "optional"
		}
    }
}

# TODO:
# 	- simplify options :
#		- remove opt : file-list (in worklow read file to convert into file list or reverse)
#		- correction for options -i default value
#       - need to test gvcf option more accuratly
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
	}

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

    parameter_meta {
        path_exe: {
			description: 'Path used as executable [default: "bcftools"]',
			category: "optional"
		}
		in: {
			description: 'Array of files to merge together (extension: ".vcf.gz|.bcf")',
			category: "required"
		}
		outputPath: {
			description: 'Path where was generated output [default: VCF path]',
			category: "optional"
		}
		forceSamples: {
			description: 'Resolve duplicate sample names [default: false]',
			category: "optional"
		}
		printHeader: {
			description: 'Print only the merged header and exit [default: false]',
			category: "optional"
		}
		useHeader: {
			description: 'Use the provided header [default: -]',
			category: "optional"
		}
		missingToRef: {
			description: 'Assume genotypes at missing sites are 0/0 [default: false]',
			category: "optional"
		}
		applyFilter: {
			description: 'Require at least one of the listed FILTER strings (e.g. "PASS,.")',
			category: "optional"
		}
		filterLogic: {
			description: 'Remove filters if some input is PASS ("x"), or apply all filters ("+")  [default: "+"]',
			category: "optional"
		}
		gvcf: {
			description: 'Merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,I16:sum,IDV:max,IMF:max (#doubt) [default: false]',
			category: "optional"
		}
		referenceFasta: {
			description: 'Reference used to merge in gvcf mode',
			category: "optional"
		}
		infoRules: {
			description: 'Rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [-]',
			category: "optional"
		}
		fileList: {
			description: 'Read file names from the file',
			category: "optional"
		}
		merge: {
			description: 'Allow multiallelic records for <snps|indels|both|all|none|id> [default: "both"]',
			category: "optional"
		}
		noVersion: {
			description: 'Do not append version and command line to the header [default: false]',
			category: "optional"
		}
		outputType: {
			description: '"b" compressed BCF; "u" uncompressed BCF; "z" compressed VCF; "v" uncompressed VCF [default: "v"]',
			category: "optional"
		}
		regions: {
			description: "Restrict to comma-separated list of regions.",
			category: "optional"
		}
		regionsFile: {
			description: "Restrict to regions listed in a file.",
			category: "optional"
		}
		threads: {
			description: "Sets the number of threads [default: 1]",
			category: "optional"
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
	}

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
		File vcf = "~{outputFile}"
  	}

    parameter_meta {
        path_exe: {
			description: "Path used as executable [default: 'bcftools']",
			category: "optional"
		}
		in: {
			description: "VCF/BCF file to left-align and normalize indels (extension: '.vcf.gz|.bcf')",
			category: "required"
		}
		outputPath: {
			description: 'Path where was generated output',
			category: "required"
		}
		name: {
			description: 'Output file base name [default: sub(basename(in),subString,"")].',
			category: 'optional'
		}
		subString: {
			description: 'Extension to remove from the input file [default: "\.(vcf|bcf)(\.gz)?$"]',
			category: 'optional'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'optional'
		}
		suffix: {
			description: 'Suffix to add for the output file (e.g sample.suffix.bam)[default: ".merge"]',
			category: 'optional'
		}
		refFasta: {
			description: 'Reference used to merge in gvcf mode',
			category: "required"
		}
		checkRef: {
			description: 'Check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [default: e]',
			category: "optional"
		}
		removeDuplicates: {
			description: 'Remove duplicate lines of the same type.',
			category: "optional"
		}
		rmDupType: {
			description: 'Remove duplicate snps|indels|both|all|none (implies removeDuplicates)',
			category: "optional"
		}
		splitMA: {
			description: "Split (true) or join (false) multiallelics sites [default= false]",
			category: "optional"
		}
		multiallelicType: {
			description: "Type of Multiallelics to treat for split/join (type: snps|indels|both|any) [default: both]",
			category: "optional"
		}
		version: {
			description: 'Append version and command line to the header [default: true]',
			category: "optional"
		}
		normalize: {
			description: 'Normalize indels (with -m or -c s) [default: true]',
			category: "optional"
		}
		outputType: {
			description: '"b" compressed BCF; "u" uncompressed BCF; "z" compressed VCF; "v" uncompressed VCF [default: "v"]',
			category: "optional"
		}
		regions: {
			description: "Restrict to comma-separated list of regions",
			category: "optional"
		}
		regionsFile: {
			description: "Restrict to regions listed in a file",
			category: "optional"
		}
		targets: {
			description: "Similar to 'regions' but streams rather than index-jumps",
			category: "optional"
		}
		targetsFile: {
			description: "Similar to 'regionsFile' but streams rather than index-jumps",
			category: "optional"
		}
		strictFilter: {
			description: "When merging (-m+), merged site is PASS only if all sites being merged PASS [default: false]",
			category: "optional"
		}
		siteWin: {
			description: "Buffer for sorting lines which changed position during realignment [default: 1000]",
			category: "optional"
		}
		threads: {
			description: "Sets the number of threads [default: 0]",
			category: "optional"
		}
	}
}
