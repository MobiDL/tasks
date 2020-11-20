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
		String path_exe = "regtools"

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
		~{path_exe} --help | grep "Program\|Version"
	>>>

	output {
		String output = read_string(stdout())
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

task junctionsExtract {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-10-22"
	}

	input {
		String path_exe = "regtools"

		File in
		File inIdx
		String? outputPath
		String? name
		String subString = ".bam"
		String subStringReplace = ""
		String suffix = ".junctions"
		String ext = ".bed"

		Int minAnchor = 8
		Int minIntronSize = 70
		Int maxIntronSize = 500000

		String? region

		Int strandSpec = 0

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
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}~{ext}" else "~{baseName}~{suffix}~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} junctions extract \
			-a ~{minAnchor} \
			-m ~{minIntronSize} \
			-M ~{maxIntronSize} \
			~{default="" "-r " + region} \
			-s ~{strandSpec} \
			-o ~{outputFile} \
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
			description: 'Path used as executable [default: "samtools"]',
			category: 'System'
		}
		in: {
			description: 'Path to the alignement file (bam)',
			category: 'Required'
		}
		inIdx: {
			description: 'Path to the alignement file index (bai)',
			category: 'Required'
		}
		outputPath: {
			description: 'Output path where output file where generated [default: pwd()]',
			category: 'Tool option'
		}
		subString: {
			description: 'Substring to replace from the input file [default: ".bam"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'Substring replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		suffix: {
			description: 'Suffix to add on the output file [default: ".junctions"]',
			category: 'Tool option'
		}
		ext: {
			description: 'Extension to add on the output file [default: ".bed"]',
			category: 'Tool option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Tool option'
		}
		minAnchor: {
			description: 'Minimum anchor length. Junctions which satisfy a minimum anchor length on both sides are reported. [default: 8]',
			category: 'Tool option'
		}
		minIntronSize: {
			description: 'Minimum intron length. [default: 70]',
			category: 'Tool option'
		}
		maxIntronSize: {
			description: 'Maximum intron length. [500000]',
			category: 'Tool option'
		}
		region: {
			description: 'The region to identify junctions in "chr:start-end" format.',
			category: 'Tool option'
		}
		strandSpec: {
			description: 'Strand specificity of RNA library preparation (0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR). [default: 0]',
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

task junctionsAnnotate {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-10-26"
	}

	input {
		String path_exe = "regtools"

		File in
		File refFasta
		File refGTF
		String? outputPath
		String? name
		String subString = ".bed"
		String subStringReplace = ""
		String suffix = ".annotate"
		String ext = ".bed"

		Boolean singleExonGenes = true

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
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}~{ext}" else "~{baseName}~{suffix}~{ext}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} junctions annotate \
			~{true="-S" false="" singleExonGenes} \
			-o ~{outputFile} \
			~{in} \
			~{refFasta} \
			~{refGTF}

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
			description: 'Path used as executable [default: "samtools"]',
			category: 'System'
		}
		in: {
			description: 'Path to the bed file containing junctions (bed)',
			category: 'Required'
		}
		outputPath: {
			description: 'Output path where output file where generated [default: pwd()]',
			category: 'Tool option'
		}
		subString: {
			description: 'Substring to replace from the input file [default: ".bed"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'Substring replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		suffix: {
			description: 'Suffix to add on the output file [default: ".annotate"]',
			category: 'Tool option'
		}
		ext: {
			description: 'Extension to add on the output file [default: ".bed"]',
			category: 'Tool option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Tool option'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'Required'
		}
		refGTF: {
			description: 'Path to the GTF reference file (format: GTF)',
			category: 'Required'
		}
		singleExonGenes: {
			description: 'Include single exon genes [default: true]',
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
