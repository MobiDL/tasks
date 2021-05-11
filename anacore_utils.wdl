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

		~{path_exe} \
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
			description: 'Bed file with panel.',
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

		~{path_exe} \
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
		bed: {
			description: 'Bed file with panel.',
			category: 'Required'
		}
		subString: {
			description: 'Extension to remove from the input file [default: ".bam"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string to create basename output files [default: ".RG"]',
			category: 'Output path/name option'
		}
		summaryTSV: {
			description: 'Ouput format is tsv (false for json) [default: true]',
			category: 'Tool option'
		}
		checkStrand: {
			description: 'With this option the strand of amplicons is checked. [default: false]',
			category: 'Tool option'
		}
		singleEnd: {
			description: 'Process single-end alignments. [default: false]',
			category: 'Tool option'
		}
		anchorOffset: {
			description: ' The alignment of the read can start at N nucleotids after the start of the primer. [default: 4]',
			category: 'Tool option'
		}
		minZOI: {
			description: 'The minimum cumulative length of reads pair in zone of interest.[default: 10]',
			category: 'Tool option'
		}
		tag: {
			description: 'RG tag used to store the area ID. [default: "LB"]',
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

task areaCoverage {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-05-11"
	}

	input {
		String path_exe = "areaCoverage.py"

		Array[File]+ in
		File bed
		String? outputPath
		String? name
		String subString = ".txt"
		String subStringReplace = ".areaCov"

		Boolean outputTSV = false
		Int percentile = 5

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(in[0]),subString,subStringReplace)
	String ext = if outputTSV then ".tsv" else ".json"
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{ext}" else "~{baseName}~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} \
			--percentile ~{percentile} \
			--inputs-depths ~{sep=" " in} \
			--input-regions ~{bed} \
			--output-metrics ~{outputFile}

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
			description: 'Path used as executable [default: "areaCoverage.py"]',
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
			description: 'The path to the depths by position (format: samtools depth output).',
			category: 'Required'
		}
		bed: {
			description: 'Bed file with panel.',
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
		percentile: {
			description: 'Only the depths for this percentile and his multiples are retained.[default: 5]',
			category: 'Tool option'
		}
		outputTSV: {
			description: 'Ouput format is tsv, or json [default: false (ie. json)]',
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

task splitBAMByRG {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.2"
		date: "2021-05-11"
	}

	input {
		String path_exe = "splitBAMByRG.py"

		File in
		File readGroup
		String? outputPath
		String? name
		String subString = ".bam"
		String subStringReplace = "_{GP}.bam"

		Boolean remove = false
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
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"
	String outputGlob = sub(basename(in),"\{GP\}",".*")

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} \
			~{true="--remove-RG" false="" remove} \
			--RG-tag ~{tag} \
			--input-aln ~{in} \
			--input-design ~{readGroup} \
			--output-pattern ~{outputFile}

	>>>

	output {
		Array[File] outputFile = glob(outputPath + "/" + outputGlob)
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "splitBAMByRG.py"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name (pattern for files creation) [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		in: {
			description: 'Alignement file used as input (bam).',
			category: 'Required'
		}
		readGroup: {
			description: 'The path to the file describing RG in each new group (format: TSV).',
			category: 'Required'
		}
		subString: {
			description: 'Extension to remove from the input file [default: ".bam"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: "_{GP}.bam"]',
			category: 'Output path/name option'
		}
		remove: {
			description: 'With this parameter the RG are removed from the outputted alignments files. [default: false]',
			category: 'Tool option'
		}
		tag: {
			description: 'RG tag used to store the area ID. [default: "LB"]',
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

task filterVCFPrimer {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-05-11"
	}

	input {
		String path_exe = "filterVCFPrimer.py"

		File in
		String? outputPath
		String? name
		String subString = ".vcf"
		String subStringReplace = ".filterPrimer.vcf"

		File bed
		File refFasta

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

		~{path_exe} \
			--input-variants ~{in} \
			--input-regions ~{bed} \
			--input-sequences ~{refFasta} \
			--output-variants ~{outputFile}

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
			description: 'Path used as executable [default: "filterVCFPrimer.py"]',
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
			description: 'Bed file with panel.',
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
		bed: {
			description: 'Bed file with panel.',
			category: 'Required'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'Required'
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
