version 1.0

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
		Int threads   = 0
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
		~{path_exe} index ~{true="--tbi" false="--csi" tabix} ~{true="--force" false="" force} ~{minShiftOpt} --threads ~{threads} -o ~{outputIndex} ~{outputRaw}

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
			description: "Sets the number of extra-threads [default: 0]",
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
		String? region
		File? regionFile
		Int threads   = 0
	}

	String useHeaderOpt = if defined(useHeader) then "--use-header ~{useHeader} " else ""
	String applyFilterOpt = if defined(applyFilter) then "--apply-filter '~{sep=',' applyFilter}' " else ""
	String gvcfOpt = if gvcf then if defined(referenceFasta) then "--gvcf ~{referenceFasta} " else "--gvcf - " else ""
	String fileListOpt = if defined(fileList) then "--file-list ~{fileList} " else ""
	String regionOpt = if defined(region) then "--region ~{region} " else ""
	String regionFileOpt = if defined(regionFile) then "--region-file ~{regionFile} " else ""

	String ext = if outputType == "v" then ".vcf" else if outputType == "b" then ".bcf.gz" else if outputType == "u" then ".bcf" else ".vcf.gz"

	String outputFile = if defined(outputPath) then "~{outputPath}/~{prefix}~{ext}" else "~{prefix}~{ext}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi
		~{path_exe} merge ~{true="--force-samples " false="" forceSamples}~{true="--print-header " false="" printHeader}~{useHeaderOpt}~{true="--missing-to-ref " false="" missingToRef}~{applyFilterOpt}~{gvcfOpt}~{fileListOpt}~{true="--no-version " false="" noVersion}~{regionOpt}~{regionFileOpt}\
			--info-rules '~{infoRules}' \
			--merge '~{merge}' \
			--filter-logic '~{filterLogic}' \
			-O '~{outputType}' \
			-o ~{outputFile} \
			~{sep=" \\\n\t" in}

	>>>

	output {
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
		region: {
			description: "Restrict to comma-separated list of regions",
			category: "optional"
		}
		regionFile: {
			description: "Restrict to regions listed in a file",
			category: "optional"
		}
		threads: {
			description: "Sets the number of extra-threads [default: 0]",
			category: "optional"
		}
    }
}
