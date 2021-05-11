version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2021 MoBiDiC
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

task nonOverlappingDesign {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-05-11"
	}

	input {
		String path_exe = "nonOverlappingDesign.py"

		File in
		String? outputPath
		String? name
		String subString = ".bed"
		String subStringReplace = ".nonOverlapped.tsv"

		Int? margin

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
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		${path_exe} \
			~{default="" "--margin " + margin} \
			--input-panel ~{in} \
			--output-design ~{outputFile}

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
			description: 'Path used as executable [default: "nonOverlappingDesign.py"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		in: {
			description: 'Alignement file used as input (sam, bam or cram).',
			category: 'Required'
		}
		subString: {
			description: 'Extension to remove from the input file [default: ".bed"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ".nonOverlapped.tsv"]',
			category: 'Output path/name option'
		}
		margin: {
			description: 'The minimum distance between two areas in same group.',
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

task addAmpliRG {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-05-11"
	}

	input {
		String path_exe = "addAmpliRG.py"

		File in
		File bed
		String? outputPath
		String? name
		String subString = ".bam"
		String subStringReplace = ".RG"

		Boolean summaryTSV = true
		Boolean checkStrand = false
		Boolean singleEnd = false

		Int anchorOffset = 4
		Int minZOI = 10
		String tag = "LB"


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
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}.bam" else "~{baseName}.bam"
	String ext = if summaryTSV then "tsv" else "json"
	String outputFileSummary = if defined(outputPath) then "~{outputPath}/~{baseName}.~{ext}" else "~{baseName}.~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		${path_exe} \
			--summary-format ~{true="tsv" false="json" summaryTSV} \
			~{true="--check-strand" false="" checkStrand} \
			~{true="--single-mode" false="" singleEnd} \
			--anchor-offset ~{anchorOffset} \
			--min-zoi-cov ~{minZOI} \
			--RG-tag ~{tag} \
			--input-aln ~{in} \
			--input-panel ~{bed} \
			--output-aln ~{outputFile} \
			--output-summary ~{outputFileSummary}

	>>>

	output {
		File outputFile = outputFile
		File outputFileSummary = outputFileSummary
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "addAmpliRG.py"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		in: {
			description: 'Alignement file used as input (bam).',
			category: 'Required'
		}
		subString: {
			description: 'Extension to remove from the input file [default: ".bam"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string to create basename output files [default: ".RG."]',
			category: 'Output path/name option'
		}
		summaryTSV: {
			description: '[default: false]',
			category: 'Tool option'
		}
		checkStrand: {
			description: '[default: false]',
			category: 'Tool option'
		}
		singleEnd: {
			description: '[default: false]',
			category: 'Tool option'
		}
		anchorOffset: {
			description: '[default: false]',
			category: 'Tool option'
		}
		minZOI: {
			description: '[default: false]',
			category: 'Tool option'
		}
		tag: {
			description: '[default: false]',
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
