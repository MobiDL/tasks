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
		String path_exe = "bedtools"

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
			description: 'Path used as executable [default: "bedtools"]',
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

# Options optimizations
#	- loj, wao: prefer option -v to have complement
#	- wo: prefer script to count
#	- u : prefer 'uniq'
#	- s and S can be merged together (not defined -> null ; true s; false S)
#	- r and e can be merged together (not defined -> null ; true r; false e)
#	- c and C : prefer uniq -c
#	- sortout and g : prefer bedtools sort
task intersect {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-28"
	}

	input {
		String path_exe = "bedtools"

		File bedA
		Array[File]+ bedB

		String? outputPath
		String? name = "intersect.bed"

		Boolean wa = false
		Boolean wb = false

		Boolean v = false

		Float? f
		Float? F

		Boolean? reciprocal
		Boolean? strandness

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	Boolean filenames = wb

	String outputFile = if defined(outputPath) then "~{outputPath}/~{name}" else "~{name}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} intersect \
			~{true="-wa" false="" wa} ~{true="-wb" false="" wb} \
			~{default="" true="-r" false="-e" reciprocal} ~{default="" true="-s" false="-S" strandness} \
			~{true="-filenames " false="" filenames}\
			~{default="" "-f " + f} \
			~{default="" "-F " + F} \
			-a ~{bedA} \
			-b ~{sep=" " bedB} > ~{outputFile}

	>>>

	output {
		File out = "~{outputFile}"
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bedtools"]',
			category: 'System'
		}
		outputPath: {
			description: 'Path where was generated output. [default: pwd(script)]',
			category: 'Output path/name option'
		}
		bedA: {
			description: 'BAM/BED/GFF/VCF file "A".',
			category: 'Required'
		}
		bedB: {
			description: 'BAM/BED/GFF/VCF files "B" (One or more).',
			category: 'Required'
		}
		wa: {
			description: 'Write the original entry in A for each overlap.',
			category: 'Tool option'
		}
		wb: {
			description: 'Write the original entry in B for each overlap.',
			category: 'Tool option'
		}
		v: {
			description: 'Only report those entries in A that have no overlap in B.',
			category: 'Tool option'
		}
		f: {
			description: 'Minimum overlap required as a fraction of A (e.g 0.1). (default: null = 1bp)',
			category: 'Tool option'
		}
		F: {
			description: 'Minimum overlap required as a fraction of B(e.g 0.9). (default: null = 1bp)',
			category: 'Tool option'
		}
		strandness: {
			description: 'Force "strandedness" (true) or "different strandness" (false). [default: null]',
			category: 'Tool option'
		}
		reciprocal: {
			description: 'true : F = f ; false : fileter OR (f OR F is OK) ; default: need to respect f AND F ',
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

task sort {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-28"
	}

	input {
		String path_exe = "bedtools"

		File in
		String ext = ".bed"

		String? outputPath
		String name = basename(in, ext)
		String suffix = "sort"

		Boolean? sortBySizeAsc
		Boolean? sortByChrSizeAsc
		Boolean? sortByChrScoreAsc

		File? idx

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String outputFile = if defined(outputPath) then "~{outputPath}/~{name}.sort~{ext}" else "~{name}.~{suffix}~{ext}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} sort \
			~{default="" true="-sizeA" false="-sizeD" sortBySizeAsc} \
			~{default="" true="-chrThenSizeA" false="-chrThenSizeD" sortByChrSizeAsc} \
			~{default="" true="-chrThenScoreA" false="-chrThenScoreD" sortByChrScoreAsc} \
			~{"-faidx " + idx} \
			-i ~{in} > ~{outputFile}

	>>>

	output {
		File out = "~{outputFile}"
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bedtools"]',
			category: 'System'
		}
		outputPath: {
			description: 'Path where was generated output. [default: pwd(script)]',
			category: 'Output path/name option'
		}
		in: {
			description: 'BED/GFF/VCF file.',
			category: 'Required'
		}
		ext: {
			description: 'Extension of the input file (BED/GFF/VCF) [default: ".bed"]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Prefix for the output file [default: basename(in, ext)]',
			category: 'Output path/name option'
		}
		suffix: {
			description: 'Suffix for the output file (e.g. name.suffix.ext) [default: "sort"]',
			category: 'Output path/name option'
		}
		sortBySizeAsc: {
			description: 'true: sort by feature size asc; false: sort by feature size desc; default: null',
			category: 'Tool option'
		}
		sortByChrSizeAsc: {
			description: 'true: sort by chrom then by feature size asc; false: sort by chrom then by feature size desc; default: null',
			category: 'Tool option'
		}
		sortByChrScoreAsc: {
			description: 'true: sort by chrom then by score asc; false: sort by chrom then by score desc; default: null',
			category: 'Tool option'
		}
		idx: {
			description: 'sort according to chromosome in file; default: null',
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

task coverage {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.2"
		date: "2021-04-28"
	}

	input {
		String path_exe = "bedtools"

		File bedA
		File bedB

		String? outputPath
		String? name = "coverage.bed"

		String? outMode # values possible => "-hist" "-d" "-counts" or "-mean"
		Boolean d = false

		Boolean v = false

		Float? f
		Float? F

		Boolean? reciprocal
		Boolean? strandness

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String outputFile = if defined(outputPath) then "~{outputPath}/~{name}" else "~{name}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} intersect \
			~{default="" true="-r" false="-e" reciprocal} \
			~{default="" true="-s" false="-S" strandness} \
			~{true="-d" false="" d} \
			~{default="" "-f " + f} \
			~{default="" "-F " + F} \
			-a ~{bedA} \
			-b ~{bedB} > ~{outputFile}

	>>>

	output {
		File out = "~{outputFile}"
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bedtools"]',
			category: 'System'
		}
		outputPath: {
			description: 'Path where was generated output. [default: pwd(script)]',
			category: 'Output path/name option'
		}
		bedA: {
			description: 'BAM/BED/GFF/VCF file "A".',
			category: 'Required'
		}
		bedB: {
			description: 'BAM/BED/GFF/VCF file "B".',
			category: 'Required'
		}
		d: {
			description: 'Report the depth at each position in each A feature. (default: false)',
			category: 'Tool option'
		}
		f: {
			description: 'Minimum overlap required as a fraction of A (e.g 0.1). (default: null = 1bp)',
			category: 'Tool option'
		}
		F: {
			description: 'Minimum overlap required as a fraction of B(e.g 0.9). (default: null = 1bp)',
			category: 'Tool option'
		}
		strandness: {
			description: 'Force "strandedness" (true) or "different strandness" (false). [default: null]',
			category: 'Tool option'
		}
		reciprocal: {
			description: 'true : F = f ; false : fileter OR (f OR F is OK) ; default: need to respect f AND F ',
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
