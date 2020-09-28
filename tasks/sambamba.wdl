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

task markdup {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample
		String suffix = ".markdup"

		Boolean removeDuplicates = false
		Int compressionLevel = 6

		String? tempDir
		Int? hashTableSize
		Int? overflowListSize
		Int? sortBufferSize
		Int? bufferSize

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.sam|\.cram)","")
	String outputBam = if defined(outputPath) then "~{outputPath}/~{sampleName}~{suffix}.bam" else "~{sampleName}~{suffix}.bam"
	String outputBai = if defined(outputPath) then "~{outputPath}/~{sampleName}~{suffix}.bam.bai" else "~{sampleName}~{suffix}.bam.bai"

	command <<<

		if [[ ! -d $(dirname ~{outputBam}) ]]; then
			mkdir -p $(dirname ~{outputBam})
		fi

		~{path_exe} markdup \
			~{true="--remove-duplicates" false="" removeDuplicates} \
			--nthreads ~{threads} \
			--compression-level ~{compressionLevel} \
			~{default="" "--tmpdir " + tempDir} \
			~{default="" "--hash-table-size " + hashTableSize} \
			~{default="" "--overflow-list-size " + overflowListSize} \
			~{default="" "--sort-buffer-size " + sortBufferSize} \
			~{default="" "--io-buffer-size " + bufferSize} \
			~{in} ~{outputBam}

	>>>

	output {
		File outputBam = outputBam
		File outputBai = outputBai
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${totalMemMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "sambamba"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where bam file was generated. [default: pwd()]',
			category: 'optional'
		}
		sample: {
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam file to mark or remove duplicates.',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. sample.suffix.bam) [default: ".markdup"]',
			category: 'optional'
		}
		removeDuplicates: {
			description: 'Remove duplicates instead of just marking them [default: false]',
			category: 'optional'
		}
		compressionLevel: {
			description: 'Specify compression level of the resulting file (from 0 to 9) [default: 6]',
			category: 'optional'
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
		hashTableSize: {
			description: 'Size of hash table for finding read pairs',
			category: 'optional'
		}
		overflowListSize: {
			description: 'Size of the overflow list where reads, thrown from the hash table, get a second chance to meet their pairs',
			category: 'optional'
		}
		sortBufferSize: {
			description: 'Total amount of memory (in *megabytes*) used for sorting purposes',
			category: 'optional'
		}
		bufferSize: {
			description: 'Two buffers of BUFFER_SIZE *megabytes* each are used for reading and writing BAM during the second pass',
			category: 'optional'
		}
	}
}

task index {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample

		Boolean checkBins = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.cram)","")
	String ext = sub(basename(in),"^.*(bam|cram)","$1")
	Boolean cram = if ext=="cram" then true else false
	String extOut = if cram then ".crai" else "" # fix the fact that sambamba write automatically
	String extFile = if cram then "" else ".bai" # extension for output crai but not for bai
	String outputIdx = if defined(outputPath) then "~{outputPath}/~{sampleName}~{extFile}" else "~{sampleName}~{extFile}"

	command <<<

		if [[ ! -d $(dirname ~{outputIdx}) ]]; then
			mkdir -p $(dirname ~{outputIdx})
		fi

		~{path_exe} index \
			--nthreads ~{threads} \
			~{true="--check-bins" false="" checkBins} \
			~{true="--cram-input" false="" cram} \
			~{in} \
			~{outputIdx}

	>>>

	output {
		File outputIdx = outputIdx + extOut
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${totalMemMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "sambamba"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where bam file was generated. [default: pwd()]',
			category: 'optional'
		}
		sample: {
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam or cram file to index.',
			category: 'Required'
		}
		checkBins: {
			description: 'check that bins are set correctly',
			category: 'optional'
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
		date: "2020-07-31"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample
		String suffix = ".sort"

		String? filter
		Boolean? sortByReadName

		Int compressionLevel = 6
		Boolean uncompressedChuncks = false

		String? tempDir

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.sam|\.cram)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{sampleName}~{suffix}.bam" else "~{sampleName}~{suffix}.bam"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} sort \
			--nthreads ~{threads} \
			~{default="" "--tmpdir " + tempDir} \
			~{default="" "--memory-limit " + memory} \
			~{default="" "--filter " + filter} \
			~{default="" true="--sort-by-name" false="--natural-sort" sortByReadName} \
			--compression-level ~{compressionLevel} \
			~{true="--uncompressed-chunks" false="" uncompressedChuncks} \
			--out ~{outputFile} \
			~{in}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${totalMemMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "sambamba"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where bam file was generated. [default: pwd()]',
			category: 'optional'
		}
		sample: {
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam file to sort.',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add on the output file (e.g. sample.suffix.bam) [default: ".sort"]',
			category: 'optional'
		}
		compressionLevel: {
			description: 'Specify compression level of the resulting file (from 0 to 9) [default: 6]',
			category: 'optional'
		}
		tempDir: {
			description: 'Directory for storing intermediate files; default is system directory for temporary files',
			category: 'optional'
		}
		filter: {
			description: 'Keep only reads that satisfy FILTER',
			category: 'optional'
		}
		sortByReadName: {
			description: 'Sort by read name instead of coordinate (true: lexicographical; false: natural) (default: null)',
			category: 'optional'
		}
		uncompressedChuncks: {
			description: 'Write sorted chunks as uncompressed BAM (default is writing with compression level 1), that might be faster in some cases but uses more disk space',
			category: 'optional'
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

task flagstat {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-05"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample
		String ext = ".flagstats"

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.sam|\.cram)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{sampleName}~{ext}" else "~{sampleName}~{ext}"


	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} flagstat \
			--nthreads ~{threads} \
			~{in} > ~{outputFile}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${totalMemMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "sambamba"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where flagstat file was generated. [default: pwd()]',
			category: 'optional'
		}
		sample: {
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam file to sort.',
			category: 'Required'
		}
		ext: {
			description: 'Extension of the output file [default: ".flagstats"]',
			category: 'optional'
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

task view {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-20"
	}

	input {
		String path_exe = "sambamba"

		File in
		String? outputPath
		String? sample
		Boolean cram = false

		String? filters
		String? numFilter

		Boolean header = true
		Boolean valid = false

		File? refFasta

		Int compressionLevel = 6
		File? regions

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String ext = if cram then ".cram" else ".bam"
	String sampleName = if defined(sample) then sample else sub(basename(in),"(\.bam|\.sam|\.cram)","")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{sampleName}~{ext}" else "~{sampleName}~{ext}"


	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} view \
			--nthreads ~{threads} \
			--compression-level ~{compressionLevel} \
			~{default="" "--filter " + filters} \
			--format ~{true="cram" false="bam" cram} \
			~{default="" "--num-filter " + numFilter} \
			~{true="--with-header" false="" header} \
			~{true="--valid" false="" valid} \
			~{default="" "--ref-filename " + refFasta} \
			~{default="" "--regions " + regions} \
			-o ~{outputFile} \
			~{in}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${totalMemMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "sambamba"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where flagstat file was generated. [default: pwd()]',
			category: 'optional'
		}
		sample: {
			description: 'Sample name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'optional'
		}
		in: {
			description: 'Bam.',
			category: 'Required'
		}
		cram: {
			description: 'Output to cram. [default: false]',
			category: 'optional'
		}
		filters: {
			description: 'Set custom filter for alignments',
			category: 'optional'
		}
		numFilter: {
			description: 'Filter flag bits; "i1/i2" corresponds to -f i1 -F i2 samtools arguments; either of the numbers can be omitted',
			category: 'optional'
		}
		header: {
			description: 'Print header before reads (always done for BAM output) [default: true]',
			category: 'optional'
		}
		valid: {
			description: 'Output only valid alignments [default: false]',
			category: 'optional'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'required'
		}
		regions:{
			description: 'Output only reads overlapping one of regions from the BED file',
			category: 'optional'

		}
		compressionLevel: {
			description: 'Specify compression level of the resulting file (from 0 to 9) [default: 6]',
			category: 'optional'
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
