version 1.0

task index {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-12"
	}

	input {
		String path_exe = "tabix"

		File in
		String? outputPath
		String? name

		Boolean csi = false
		String? preset
		Int? minShift
		Boolean tes =true
		Boolean tes2 =true

		Int threads = 1
	}

	String extIn = sub(basename(in),".*\.(gff|bed|sam|vcf)\.gz$","$1")
	Boolean csiUsed = (csi || defined(minShift))
	String extOut = if csiUsed then ".csi" else ".tbi"

	String baseName = if defined(name) then name else basename(in)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{extOut}" else "~{baseName}~{extOut}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} \
			~{default="" "--preset " + preset} \
			~{default="" "--min-shift " + minShift} \
			~{true="--csi" false="" csiUsed} \
			~{in}

		mv ~{in}~{extOut} ~{outputFile}

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
		name: {
			description: 'Name to use for output file name [default: basename(in)]',
			category: 'optional'
		}
		in: {
			description: 'File to index with tabix.',
			category: 'Required'
		}
		csi: {
			description: 'Generate CSI index for VCF (otherwise is TBI) [default: false]',
			category: 'optional'
		}
		preset: {
			description: 'Input format for indexing. Valid values are: gff, bed, sam, vcf.',
			category: 'optional'
		}
		minShift: {
			description: 'Set minimal interval size for CSI indices to 2^INT (if  defined force to use csi)',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}
