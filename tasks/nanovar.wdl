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
	author: "Thomas GUIGNARD"
		email: "t-guignard(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-19"
	}

	input {
		String path_exe = "nanovar"

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
			description: 'Path used as executable [default: "mytool"]',
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

task nanovar {
	meta {
		author: "Thomas GUIGNARD"
		email: "t-guignard(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-19"
	}

	input {

		String path_exe = "nanovar"
		String path_exe_awk = "awk"

		File inputNanovar

		File refFasta
		File refGenomeIndex
		String gapGenomeBuild
		Float score = 1.0
		String outputPath
		String genomeBuild
		#String? name

		Int minCov = 2
		Int minAlign = 200
		Float splitPct = 0.05

#		String subString = "(substring)"
#		String subStringReplace = "regexp"

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

#	String baseName = if defined(name) then name else sub(basename(fastqR1),subString,subStringReplace)
#	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		~{path_exe} \
	--data_type ont \
	--threads ~{threads} \
	--filter_bed ~{gapGenomeBuild} \
	--score ~{score} \
	--mincov ~{minCov} \
	--minalign ~{minAlign} \
	--splitpct ~{splitPct} \
	~{inputNanovar} \
	~{refFasta} \
	~{outputPath}

	>>>

	output {
	#	File output = outputRep
	File ~{outputPath}/nanovar_run/hsblast_longreads/ALL.hsblast-~{genomeBuild}.tsv
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	# generate documentation
	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "nanovar"]',
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
		minCov: {
			description: 'minimum number of reads required to call a breakend [2]',
			category: 'Option'
		}
		minAlign: {
				description: 'minimum alignment length for single alignment reads [200]',
				category: 'Option'
			}
		splitPct: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'Option'
		}
		gapGenomeBuild: {
			description: 'BED file with genomic regions to be excluded. (e.g. telomeres and centromeres) Either specify name of in-built reference genome filter (i.e. hg38, hg19, mm10) or provide FULL path to own BED file. [None]',
			category: 'Option'
		}
		score: {
			description: 'score threshold for defining PASS/FAIL SVs in VCF. Default score 1.0 was estimated from simulated analysis. [1.0]',
			category: 'Option'
		}
		refFasta: {
			description: 'Path to reference genome in FASTA. Genome indexes created will overwrite indexes created by other aligners (e.g. bwa)',
			category: 'Option'
		}
		inputNanovar: {
			description: 'Path to long reads or mapped BAM file. Formats: fasta/fa/fa.gzip/fa.gz/fastq/fq/fq.gzip/fq.gz or .bam',
			category: 'Input'
		}
		outputPath: {
			description: 'Path to working directory. Directory will be created if it does not exist',
			category: 'Input'
		}
	}
}
