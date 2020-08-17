version 1.0

task bgzip {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-24"
	}

	input {
		File in
		String outputPath

		Boolean decompress = false
		Boolean force = false
		Boolean index = false
		Boolean keepFile = false
		Int threads = 1
	}

	String outputFile = if ! decompress then "~{outputPath}/" + basename(in) + ".gz" else "~{outputPath}/" + basename(in,".gz")
	String outputIndex = if (index && ! decompress) then "~{outputFile}.gzi" else "~{outputFile}"
	String indexOpt = if (index && ! decompress) then "--index --index-name ~{outputIndex}" else ""

	command <<<

		bgzip \
			~{true="--stdout" false="" keepFile} \
			~{true="--decompress" false="" decompress} \
			~{true="--force" false="" force}\
			~{indexOpt} \
			~{in} \
			--threads ~{threads} > ~{outputFile}

	>>>

	output {
		File out = "~{outputFile}"
		File index = "~{outputIndex}"
	}

    parameter_meta {
        in: {
			description: "File to compres/decompress.",
			category: "required"
		}
        outputPath: {
			description: "Path where was generated output.",
			category: "required"
		}
        decompress: {
			description: "Decompress file (incompatible with index) [default: false]",
			category: "optional"
		}
        force: {
			description: "Overwrite files without asking [default: false]",
			category: "optional"
		}
        index: {
			description: "Compress and create BGZF index [default: false]",
			category: "optional"
		}
		threads: {
			description: "Sets the number of threads [default: 1]",
			category: "optional"
		}
	}
}
