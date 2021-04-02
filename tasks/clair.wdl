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


task callVarBam {
	meta {
		author: "David BAUX"
		email: "d-baux(at)chu-montpellier.fr"
		version: "0.0.5"
		date: "2021-04-02"
	}

	input {
		String path_exe = "clair.py"
		String? outputPath
		String? name
		String subString = "(\.bam)"
		String subStringReplace = ""

		String modelPath
		File refGenome
		File refGenomeIndex
		File bamFile
		File bamFileIndex
		String sampleName = "sample"

		File? bedRegions
		File? candidateVcf

		Float threshold = 0.125000
		Int minCoverage = 4
		Int? qual
		String? contigName
		Int? ctgStart
		Int? ctgEnd
		Boolean stopConsiderLeftEdge = false
		Int dcov = 250
		String samtools = "samtools"
		String pypy = "pypy"
		Boolean delay = false
		Boolean debug = false
		Boolean pysamForAllIndelBases = false
		Boolean haploidPrecision = false
		Boolean haploidSensitive = false
		Boolean activation_only = false
		Int maxPlot = 10
		Boolean logPath = false
		Int parallelLevel = 2
		Boolean fastPlotting = false
		Int workers = 8
		Boolean outputForEnsemble = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(bamFile),subString,subStringReplace)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}.clair.vcf" else "~{baseName}.clair.vcf"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} callVarBam \
			--chkpnt_fn ~{modelPath} \
			--ref_fn ~{refGenome} \
			--bed_fn ~{bedRegions} \
			--bam_fn ~{bamFile} \
			--call_fn ~{outputFile} \
			--vcf_fn ~{candidateVcf} \
			--threshold ~{threshold} \
			--minCoverage ~{minCoverage} \
			--qual ~{qual} \
			--sampleName ~{sampleName} \
			--ctgName ~{contigName} \
			--ctgStart ~{ctgStart} \
			--ctgEnd ~{ctgEnd} \
			--stop_consider_left_edge ~{stopConsiderLeftEdge} \
			--dcov ~{dcov} \
			--samtools ~{samtools} \
			--pypy ~{pypy} \
			--delay ~{delay} \
			--debug ~{debug} \
			--pysam_for_all_indel_bases ~{pysamForAllIndelBases} \
			--haploid_precision ~{haploidPrecision} \
			--haploid_sensitive ~{haploidSensitive} \
			--activation_only ~{activation_only} \
			--max_plot ~{maxPlot} \
			--log_path ~{logPath} \
			--parallel_level ~{parallelLevel} \
			--fast_plotting ~{fastPlotting} \
			--workers ~{workers} \
			--output_for_ensemble ~{outputForEnsemble} \
			--threads ~{threads}

	>>>

	output {
		File outputFile = "~{outputFile}"
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	# generate documentation
	parameter_meta {
		path_exe: {
			description: 'Path to CLAIR python script [default: "clair.py"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where vcf file was generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Sample name to use for output file name [default: sub(basename(bamFile),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Substring to remove to get sample name [default: "(\.bam)"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		modelPath: {
			description: 'Path to the model folder',
			category: 'input'
		}
		refGenome: {
			description: 'Reference fasta file input, [default: ref.fa]',
			category: 'input'
		}
		refGenomeIndex: {
			description: 'Index of the FASTA reference.',
			category: 'input'
		}
		bedRegions: {
			description: 'Call variant only in these regions, works in intersection with ctgName, ctgStart and ctgEnd, optional, [default: as defined by ctgName, ctgStart and ctgEnd]',
			category: 'Option'
		}
		bamFile: {
			description: 'BAM file input',
			category: 'input'
		}
		bamFileIndex: {
			description: 'Index of the bam file',
			category: 'input'
		}
		candidateVcf: {
			description: 'Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file, [default: None]',
			category: 'Option'
		}
		threshold: {
			description: 'Minimum allele frequence of the 1st non-reference allele for a site to be considered as a condidate site, [default: 0.125000]',
			category: 'Option'
		}
		minCoverage: {
			description: 'Minimum coverage required to call a variant, [default: 4]',
			category: 'Option'
		}
		qual: {
			description: 'If set, variant with equal or higher quality will be marked PASS, or LowQual otherwise, optional',
			category: 'Option'
		}
		sampleName: {
			description: 'Define the sample name to be shown in the VCF file',
			category: 'input'
		}
		contigName: {
			description: 'The name of sequence to be processed, [default: None]',
			category: 'Option'
		}
		ctgStart: {
			description: 'The 1-based starting position of the sequence to be processed',
			category: 'Option'
		}
		ctgEnd: {
			description: 'The 1-based inclusive ending position of the sequence to be processed',
			category: 'Option'
		}
		stopConsiderLeftEdge: {
			description: 'If not set, would consider left edge only. That is, count the left-most base-pairs of a read for coverage even if the starting position of a read is after the starting position of a tensor',
			category: 'Option'
		}
		dcov: {
			description: 'Cap depth per position at 250, [default: 250]',
			category: 'Option'
		}
		samtools: {
			description: 'Path to the "samtools", [default: samtools]',
			category: 'Option'
		}
		pypy: {
			description: 'Path to the "pypy", [default: pypy3]',
			category: 'Option'
		}
		delay: {
			description: 'Wait a short while for no more than 10 to start the job. This is to avoid starting multiple jobs simultaneously that might use up the maximum number of threads allowed, because Tensorflow will create more threads than needed at the beginning of running the program.',
			category: 'Option'
		}
		debug: {
			description: 'Debug mode, optional',
			category: 'Option'
		}
		pysamForAllIndelBases: {
			description: 'Always using pysam for outputting indel bases, optional',
			category: 'Option'
		}
		haploidPrecision: {
			description: 'call haploid instead of diploid (output homo-variant only)',
			category: 'Option'
		}
		haploidSensitive: {
			description: 'call haploid instead of diploid (output non-multi-variant only)',
			category: 'Option'
		}
		activation_only: {
			description: 'Output activation only, no prediction',
			category: 'Option'
		}
		maxPlot: {
			description: 'The maximum number of plots output, negative number, means no limit (plot all), [default: 10]',
			category: 'Option'
		}
		logPath: {
			description: 'The path for tensorflow logging, [default: None]',
			category: 'Option'
		}
		parallelLevel: {
			description: 'The level of parallelism in plotting (currently available: 0, 2), [default: 2]',
			category: 'Option'
		}
		fastPlotting: {
			description: 'Enable fast plotting.',
			category: 'Option'
		}
		workers: {
			description: 'The number of workers in plotting, [default: 8]',
			category: 'Option'
		}
		outputForEnsemble: {
			description: 'Output for ensemble',
			category: 'Option'
		}
		threads: {
			description: 'Number of threads, optional',
			category: 'System'
		}
		memoryByThreads: {
			description: '',
			category: 'System'
		}
		memory: {
			description: '',
			category: 'System'
		}
	}
}
