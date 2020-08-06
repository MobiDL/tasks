version 1.0

task bamqc {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-05"
	}

	input {
		String path_exe = "qualimap"

		File in
		String? outputPath
		String? name

		Boolean chromLimit = true
		File? featureFile
		Int minHomopolymerSize = 3
		Boolean collectOverlapPairs = false

		Boolean outsideStats = true
		Boolean pdf = true


		Int nr = 1000
		Int nWindows = 400
		Int threads = 1
		String? memory
	}

	Int memTemp = 768*threads
	String totalMem = if defined(memory) then memory else "~{memTemp}M"

	String baseName = if defined(name) then name else sub(basename(in),"(\.sam|\.bam|\.cram)","")
	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	String pdfReport = outputRep + "/report.pdf"
	String htmlReport = outputRep + "/qualimapReport.html"
	String htmlFullReport = outputRep + "/qualimapReportOutsideRegions.html"

	command <<<

		if [[ ! -d ~{outputRep} ]]; then
			mkdir -p ~{outputRep}
		fi

		~{path_exe} bamqc \
			~{true="--paint-chromosome-limits" false="" chromLimit} \
			~{default="" "--feature-file " + featureFile} \
			-hm ~{minHomopolymerSize} \
			~{true="--collect-overlap-pairs" false="" collectOverlapPairs} \
			~{true="--outside-stats" false="" outsideStats} \
			-outdir ~{outputRep} \
			-outformat ~{true="PDF:HTML" false="HTML" pdf} \
			-nr ~{nr} \
			-nt ~{threads} \
			-nw ~{nWindows} \
			--java-mem-size=~{totalMem} \
			-bam ~{in}

	>>>

	output {
		File htmlReport = htmlReport
		File? pdfReport = pdfReport
		File? htmlFullReport = htmlFullReport
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "qualimap"]',
			category: 'optional'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'optional'
		}
		name: {
			description: 'Name to use for output repertory name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam file to analyze.',
			category: 'Required'
		}
		chromLimit: {
			description: 'Paint chromosome limits inside charts [default: true]',
			category: 'optional'
		}
		featureFile: {
			description: 'Feature file with regions of interest in GFF/GTF or BED format',
			category: 'optional'
		}
		minHomopolymerSize: {
			description: 'Minimum size for a homopolymer to be considered in indel analysis[default: 3]',
			category: 'optional'
		}
		collectOverlapPairs: {
			description: 'Activate this option to collect statistics of overlapping paired-end reads [default: false]',
			category: 'optional'
		}
		outsideStats: {
			description: 'Report information for the regions outside those defined by feature-file (ignored if no feature file defined) [default: true]',
			category: 'optional'
		}
		pdf: {
			description: 'Specify if a pdf report will be generated [default: true]',
			category: 'optional'
		}
		nr: {
			description: 'Sets number of reads analyzed in a chunk [default: 1000]',
			category: 'optional'
		}
		nWindows: {
			description: 'Sets number of reads analyzed in a windows [default: 400]',
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