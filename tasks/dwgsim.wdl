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
		date: "2020-11-20"
	}

	input {
		String path_exe = "dwgsim"

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
		~{path_exe} 2>&1 | grep "Program\|Version\|Contact"
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
			description: 'Path used as executable [default: "dwgsim"]',
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

task simulateReadsIllumina {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.2"
		date: "2020-08-05"
	}

	input {
		String path_exe = "dwgsim"

		File in
		String? outputPath
		String name = "dwgsim-illumina"

		Int sizeR1 = 150
		Int sizeR2 = sizeR1

		Boolean useInnerDistance = false

		Int insertSize = 250
		Int meanCov = 100
		Int stdInsert = 10

		Float rateMut = 0.0010
		Float freqMut = 0.5000

		Float rateIndels = 0.1000
		Float probExt = 0.3000
		Int minSizeIndels = 1

		Float probRandomDNA = 0.05
		Int maxN = 0

		File? target
		String? readPrefix

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String outputPrefix = if defined(outputPath) then "~{outputPath}/~{name}" else "~{name}"

	command <<<

		if [[ ! -d $(dirname ~{outputPrefix}) ]]; then
			mkdir -p $(dirname ~{outputPrefix})
		fi

		~{path_exe} \
			-1 ~{sizeR1} \
			-2 ~{sizeR2} \
			-C ~{meanCov} \
			~{true="-i" false="" useInnerDistance} \
			-d ~{insertSize} \
			-s ~{stdInsert} \
			-r ~{rateMut} \
			-F ~{freqMut} \
			-R ~{rateIndels} \
			-X ~{probExt} \
			-I ~{minSizeIndels} \
			-y ~{probRandomDNA} \
			-n ~{maxN} \
			~{default="" "-x " + target} \
			~{default="" "-P " + readPrefix} \
			~{in} \
			~{outputPrefix}
	>>>

	output {
		File bfastFQ = outputPrefix + ".bfast.fastq"
		File fastqR1 = outputPrefix + ".bwa.read1.fastq"
		File fastqR2 = outputPrefix + ".bwa.read2.fastq"
		File mutationsTXT = outputPrefix + ".mutations.txt"
		File mutationsVCF = outputPrefix + ".mutations.vcf"
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "dwgsim"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name [default: dwgsim-illumina]',
			category: 'Output path/name option'
		}
		in: {
			description: 'Fasta file used to create fastq.',
			category: 'Required'
		}
		sizeR1: {
			description: 'Length of the first read [default: 150]',
			category: 'Tool option'
		}
		sizeR2: {
			description: 'Length of the second read [default: sizeR1]',
			category: 'Tool option'
		}
		meanCov: {
			description: 'mean coverage across available positions [default: 100]',
			category: 'Tool option'
		}
		useInnerDistance: {
			description: 'Use the inner distance instead of the outer distance for pairs [default: false]',
			category: 'Tool option'
		}
		insertSize: {
			description: 'Outer distance between the two ends for pairs [default: 500]',
			category: 'Tool option'
		}
		stdInsert: {
			description: 'Standard deviation of the distance for pairs [default: 50]',
			category: 'Tool option'
		}
		rateMut: {
			description: 'Rate of mutations [default: 0.0010]',
			category: 'Tool option'
		}
		freqMut: {
			description: 'Frequency of given mutation to simulate low fequency somatic mutations [default: 0.5000]',
			category: 'Tool option'
		}
		rateIndels: {
			description: 'Fraction of mutations that are indels [default: 0.1000]',
			category: 'Tool option'
		}
		probExt: {
			description: 'Probability an indel is extended [default: 0.3000]',
			category: 'Tool option'
		}
		minSizeIndels: {
			description: 'The minimum length indel [default: 1]',
			category: 'Tool option'
		}
		probRandomDNA: {
			description: 'Probability of a random DNA read [default: 0.05]',
			category: 'Tool option'
		}
		maxN: {
			description: 'Maximum number of Ns allowed in a given read [default: 0]',
			category: 'Tool option'
		}
		target: {
			description: 'Bed of regions to cover [default: WGS]',
			category: 'Tool option'
		}
		readPrefix: {
			description: 'Read prefix to prepend to each read name',
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
