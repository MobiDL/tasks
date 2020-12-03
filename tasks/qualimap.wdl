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
		String path_exe = "qualimap"

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
		~{path_exe} --help | grep "QualiMap\|Built"
	>>>

	output {
		String version = read_string(stdout())
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "qualimap"]',
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
		Boolean skipDuplicates = false

		Int nr = 1000
		Int nWindows = 400
		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(in),"(\.sam|\.bam|\.cram)","")
	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

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
			~{true="-sd" false="" skipDuplicates} \
			-nr ~{nr} \
			-nt ~{threads} \
			-nw ~{nWindows} \
			--java-mem-size=~{totalMem} \
			-bam ~{in}

	>>>

	output {
		File? pdfReport = outputRep + "/report.pdf"
		File? htmlReport = outputRep + "/qualimapReport.html"
		File? htmlFullReport = outputRep + "/qualimapReportOutsideRegions.html"
		Array[File] reports = select_all([
			htmlReport,
			pdfReport,
			htmlFullReport
		])
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "qualimap"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output repertory name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'Output path/name option'
		}
		in: {
			description: 'Bam file to analyze.',
			category: 'Required'
		}
		chromLimit: {
			description: 'Paint chromosome limits inside charts [default: true]',
			category: 'Tool option'
		}
		featureFile: {
			description: 'Feature file with regions of interest in GFF/GTF or BED format',
			category: 'Tool option'
		}
		minHomopolymerSize: {
			description: 'Minimum size for a homopolymer to be considered in indel analysis[default: 3]',
			category: 'Tool option'
		}
		collectOverlapPairs: {
			description: 'Activate this option to collect statistics of overlapping paired-end reads [default: false]',
			category: 'Tool option'
		}
		outsideStats: {
			description: 'Report information for the regions outside those defined by feature-file (ignored if no feature file defined) [default: true]',
			category: 'Tool option'
		}
		pdf: {
			description: 'Specify if a pdf report will be generated [default: true]',
			category: 'Tool option'
		}
		skipDuplicates: {
			description: 'Activate this option to skip duplicated alignments from the analysis. [default: false]',
			category: 'Tool option'
		}
		nr: {
			description: 'Sets number of reads analyzed in a chunk [default: 1000]',
			category: 'Tool option'
		}
		nWindows: {
			description: 'Sets number of reads analyzed in a windows [default: 400]',
			category: 'Tool option'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (768M*threads)]',
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
