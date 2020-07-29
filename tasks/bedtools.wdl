version 1.0

# Options optimizations
#	- loj, wao: prefer option -v to have complement
#	- wo: prefer script to count
#	- u : prefer 'uniq'
#	- s and S can be merged together (not defined -> null ; true s; false S)
#	- r and e can be merged together (not defined -> null ; true r; false e)
#	- c and C : prefer uniq -c
#	- sortout and g : prefer bedtools sort
task intersect {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-28"
	}

	input {
		String path_exe = "bedtools"

		File bedA
		Array[File]+ bedB

		String? outputPath
		String? name = "intersect.bed"

		Boolean wa = false
		Boolean wb = false

		Boolean v = false

		Float? f
		Float? F

		Boolean? reciprocal
		Boolean? strandness
	}

	Boolean filenames = wb

	String outputFile = if defined(outputPath) then "~{outputPath}/~{name}" else "~{name}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} intersect ~{true="-wa " false="" wa}~{true="-wb " false="" wb}\
			~{default="" true="-r " false="-e " reciprocal}~{default="" true="-s " false="-S " strandness}\
			~{true="-filenames " false="" filenames}\
			~{default="" "-f " + f} \
			~{default="" "-F " + F} \
			-a ~{bedA} \
			-b ~{sep=" " bedB} > ~{outputFile}

	>>>

	output {
		File out = "~{outputFile}"
	}

	parameter_meta {
        path_exe: {
			description: "Path used as executable [default: 'bedtools']",
			category: "Optional"
		}
		outputPath: {
			description: 'Path where was generated output. [default: pwd(script)]',
			category: "Optional"
		}
		bedA: {
			description: 'BAM/BED/GFF/VCF file "A".',
			category: "Required"
		}
		bedB: {
			description: 'BAM/BED/GFF/VCF files "B" (One or more).',
			category: "Required"
		}
        wa: {
			description: 'Write the original entry in A for each overlap.',
			category: "Optional"
		}
		wb: {
			description: 'Write the original entry in B for each overlap.',
			category: "Optional"
		}
        v: {
			description: 'Only report those entries in A that have no overlap in B.',
			category: "Optional"
		}
		f: {
			description: 'Minimum overlap required as a fraction of A (e.g 0.1). (default: null = 1bp)',
			category: "Optional"
		}
        F: {
			description: 'Minimum overlap required as a fraction of B(e.g 0.9). (default: null = 1bp)',
			category: "Optional"
		}
		strandness: {
			description: 'Force "strandedness" (true) or "different strandness" (false). [default: null]',
			category: "Optional"
		}
        reciprocal: {
			description: 'true : F = f ; false : fileter OR (f OR F is OK) ; default: need to respect f AND F ',
			category: "Optional"
		}
	}
}

task sort {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-28"
	}

	input {
		String path_exe = "bedtools"

		File in
		String ext = ".bed"

		String? outputPath
		String name = basename(in, ext)
		String suffix = "sort"

		Boolean? sortBySizeAsc
		Boolean? sortByChrSizeAsc
		Boolean? sortByChrScoreAsc

		File? idx
	}

	String outputFile = if defined(outputPath) then "~{outputPath}/~{name}.sort~{ext}" else "~{name}.~{suffix}~{ext}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} sort ~{default="" true="-sizeA " false="-sizeD " sortBySizeAsc}\
			~{default="" true="-chrThenSizeA " false="-chrThenSizeD " sortByChrSizeAsc}\
			~{default="" true="-chrThenScoreA " false="-chrThenScoreD " sortByChrScoreAsc}\
			~{"-faidx " + idx} \
			-i ~{in}  > ~{outputFile}

	>>>

	output {
		File out = "~{outputFile}"
	}

	parameter_meta {
        path_exe: {
			description: "Path used as executable [default: 'bedtools']",
			category: "Optional"
		}
		outputPath: {
			description: 'Path where was generated output. [default: pwd(script)]',
			category: "Optional"
		}
		in: {
			description: 'BED/GFF/VCF file.',
			category: "Required"
		}
		ext: {
			description: 'Extension of the input file (BED/GFF/VCF) [default: ".bed"]',
			category: 'optional'
		}
		name: {
			description: 'Prefix for the output file [default: basename(in, ext)]',
			category: 'optional'
		}
		suffix: {
			description: 'Suffix for the output file (e.g. name.suffix.ext) [default: "sort"]',
			category: 'optional'
		}
		sortBySizeAsc: {
			description: 'true: sort by feature size asc; false: sort by feature size desc; default: null',
			category: 'optional'
		}
		sortByChrSizeAsc: {
			description: 'true: sort by chrom then by feature size asc; false: sort by chrom then by feature size desc; default: null',
			category: 'optional'
		}
		sortByChrScoreAsc: {
			description: 'true: sort by chrom then by score asc; false: sort by chrom then by score desc; default: null',
			category: 'optional'
		}
		idx: {
			description: 'sort according to chromosome in file; default: null',
			category: 'optional'
		}
	}
}
