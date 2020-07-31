version 1.0

# Create a task crumbleAdvanced for all options
task crumble {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-31"
	}

	input {
		String path_exe = "crumble"

		File in
		String? outputPath
		String? name
		String suffix = ".crumble"

		String? outputFormat
		Boolean addPGHeader = true

		Int compressionLevel = 8

		Int threads = 1
	}

	String ext = if defined(outputFormat) then outputFormat else sub(basename(in),"(.*)\.(bam|cram|sam)","$2")
	String baseName = if defined(name) then name else sub(basename(in),"(\.bam|\.cram|\.sam)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}.~{ext}" else "~{baseName}~{suffix}.~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		echo '''
		~{path_exe} \
			-~{compressionLevel} \
			-O ~{ext},nthreads=~{threads} \
			~{true="" false="-z" addPGHeader} \
			~{in} \
			~{outputFile}
		'''

	>>>

}
