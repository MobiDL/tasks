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
		date: "2021-03-24"
	}

	input {
		String path_exe = "rtg"

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
		~{path_exe} --version
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
			description: 'Path used as executable [default: "bcftools"]',
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

task fasta2sdf {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-19"
	}

	input {
		String path_exe = "rtg"

		String? outputPath
		String? sample
		String subString = ".(fa|fasta)"
		String subStringReplace = "_sdf"

		File inputFile
		Boolean protein = false
		Boolean duster = false
		Array[String]? exclude
		Boolean allowDuplicates = false
		Boolean name = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(sample) then sample else sub(basename(inputFile),subString,subStringReplace)
	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	Boolean excludeDefined = defined(exclude)

	command <<<

		if [[ ! -d $(dirname ~{outputRep}) ]]; then
			mkdir -p $(dirname ~{outputRep})
		fi

		~{path_exe} format \
			--format fasta \
			--output ~{outputRep} \
			~{true="--protein" false="" protein} \
			~{true="--duster" false="" duster} \
			~{true="--exclude " false="" excludeDefined}~{sep="--exclude " exclude} \
			~{true="--allow-duplicate-names" false="" allowDuplicates} \
			~{true="" false="--no-names" name} \
			~{inputFile}

	>>>

	output {
		File outputDone = outputRep + "/done"
		File outputLog = outputRep + "/format.log"
		File outputMainIndex = outputRep + "/mainIndex"
		File outputRef = outputRep + "/reference.txt"
		File outputSummary = outputRep + "/summary.txt"
		Array[File]? outputNameIndex = glob(outputPath + "/nameIndex[0-9]+")
		Array[File]? outputNameData = glob(outputPath + "/namedata[0-9]+")
		Array[File]? outputNamePointer = glob(outputPath + "/namepointer[0-9]+")
		Array[File] outputSeqData = glob(outputPath + "/seqdata[0-9]+")
		Array[File] outputSeqPointer = glob(outputPath + "/seqpointer[0-9]+")
		Array[File] outputSequenceIndex = glob(outputPath + "/sequenceIndex[0-9]+")
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "rtg"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where bam file was generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		sample: {
			description: 'Sample name to use for output file name [default: sub(basename(inputFile),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Substring to remove to get sample name [default: ".(fa|fasta)"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: "_sdf"]',
			category: 'Output path/name option'
		}
		inputFile: {
			description: 'Input sequence file.',
			category: 'Required'
		}
		protein: {
			description: 'Input is protein. [default: false]',
			category: 'Options'
		}
		duster: {
			description: 'Treat lower case residues as unknowns [default: false]',
			category: 'Options'
		}
		exclude: {
			description: 'Exclude input sequences based on their name. (array)',
			category: 'Options'
		}
		allowDuplicates: {
			description: 'Disable checking for duplicate sequence names [default: false]',
			category: 'Options'
		}
		name: {
			description: 'Include name data in the SDF output [default: true]',
			category: 'Options'
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
