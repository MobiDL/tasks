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

task get_version {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-04-29"
	}

	input {
		String path_exe = "seqkit"

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
		~{path_exe} version
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
			description: 'Path used as executable [default: "seqkit"]',
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

task seq_filter {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-04-29"
	}

	input {
		String path_exe = "seqkit"

		File in
		String? outputPath
		String? name
		String subString = ".fastq"
		String subStringReplace = ".filter.fastq"

		Int maxLen = -1
		Int maxQual = -1
		Int minLen = -1
		Int minQual = -1

		String idRegexp = "^(\\S+)\\s?"
		String seqType = "auto"

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

		~{path_exe} seq \
			--max-len ~{maxLen} \
			--max-qual ~{maxQual} \
			--min-len ~{minLen} \
			--min-qual ~{minQual} \
			--out-file ~{outputFile} \
			--threads ~{threads} \
			~{in}

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
			description: 'Path used as executable [default: "seqkit"]',
			category: 'System'
		}
		in: {
			description: "File used as input",
			category: 'Required'
		}
		outputPath: {
			description: 'Output path',
			category: 'Output path/name option'
		}
		name: {
			description: 'Output name',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Substring to edit from input file to create output [default: ".fastq"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'Substring used to edit input file [default: ".filter.fastq"]',
			category: 'Output path/name option'
		}
		maxLen: {
			description: 'Only print sequences shorter than the maximum length (-1 for no limit) [default: -1]',
			category: 'Option filter'
		}
		maxQual: {
			description: 'Only print sequences with average quality less than this limit (-1 for no limit) [default: -1]',
			category: 'Option filter'
		}
		minLen: {
			description: 'Only print sequences longer than the minimum length (-1 for no limit) [default: -1]',
			category: 'Option filter'
		}
		minQual: {
			description: 'Only print sequences with average quality greater or equal than this limit (-1 for no limit) [default: -1]',
			category: 'Option filter'
		}
		idRegexp: {
			description: 'Regular expression for parsing ID [default: "^(\\S+)\\s?"]',
			category: 'Option global'
		}
		seqType: {
			description: 'Sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) [default: "auto"]',
			category: 'Option global'
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
