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

task vcfEval {
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
		String subString = ".vcf.gz$"
		String subStringReplace = ""

		File legacyVCF
		File legacyVCFidx
		File testVCF
		File testVCFidx
		File genomeTemplate

		File? bedRegion

		Boolean allRecords = false
		Boolean decompose = false
		Boolean refOverlap = false
		Int ploidy = 2
		Boolean squashPloidy = false

		Float? precision
		Float? sensitivity
		Boolean roc = true
		String outMode = "split"
		Boolean sortAsc = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(sample) then sample else sub(basename(testVCF),subString,subStringReplace)
	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d $(dirname ~{outputRep}) ]]; then
			mkdir -p $(dirname ~{outputRep})
		fi

		~{path_exe} vcfeval \
			--baseline ~{legacyVCF} \
			--calls ~{testVCF} \
			--template ~{genomeTemplate} \
			~{default="" "--bed-regions " + bedRegion} \
			~{true="--all-records" false="" allRecords} \
			~{true="--decompose" false="" decompose} \
			~{true="--ref-overlap" false="" refOverlap} \
			--sample-ploid ~{ploidy} \
			~{true="--squash-ploidy" false="" squashPloidy} \
			~{default="" "--at-precision" + precision} \
			~{default="" "--at-sensitivity" + sensitivity} \
			~{true="" false="--no-roc" roc} \
			--output-mode ~{outMode} \
			--sort-order ~{true="ascending" false="descending" sortAsc} \
			--threads ~{threads} \
			--output ~{outputRep}

	>>>

	output {
		File falseNegVCF = outputRep + "/fn.vcf.gz"
		File falseNegVCFidx = outputRep + "/fn.vcf.gz.tbi"
		File falsePosVCF = outputRep + "/fp.vcf.gz"
		File falsePosVCFidx = outputRep + "/fp.vcf.gz.tbi"
		File truePosBaseVCF = outputRep + "/tp-baseline.vcf.gz"
		File truePosBaseVCFidx = outputRep + "/tp-baseline.vcf.gz.tbi"
		File truePosVCF = outputRep + "/tp.vcf.gz"
		File truePosVCFidx = outputRep + "/tp.vcf.gz.tbi"
		File phasing = outputRep + "/phasing.txt"
		File progress = outputRep + "/progress"
		File summary = outputRep + "/summary.txt"
		File log = outputRep + "/vcfeval.log"
		File done = outputRep + "/done"
		File? rocSNP = outputRep + "/snp_roc.tsv.gz"
		File? rocNonSNP = outputRep + "/non_snp_roc.tsv.gz"
		File? rocWeighted = outputRep + "/weighted_roc.tsv.gz"
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
			description: 'Substring to remove to get sample name [default: ".vcf"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		legacyVCF: {
			description: 'Legacy VCF to use as baseline.',
			category: 'Required'
		}
		legacyVCFidx: {
			description: 'Index for the legacy vcf input.',
			category: 'Required'
		}
		testVCF: {
			description: 'VCF to check against legacy.',
			category: 'Required'
		}
		testVCFidx: {
			description: 'Index for the test vcf input.',
			category: 'Required'
		}
		genomeTemplate: {
			description: 'SDF of the reference genome the variants are called against',
			category: 'Required'
		}
		bedRegion: {
			description: 'Only read VCF records that overlap the ranges contained in the specified BED file',
			category: 'Option: Filtering'
		}
		allRecords: {
			description: 'Use all records regardless of FILTER status (Default is to only process records where FILTER is "." or "PASS") [default: false]',
			category: 'Option: Filtering'
		}
		decompose: {
			description: 'Decompose complex variants into smaller constituents to allow partial credit [default: false]',
			category: 'Option: Filtering'
		}
		refOverlap: {
			description: 'Allow alleles to overlap where bases of either allele are same-as-ref (Default is to only allow VCF anchor base overlap) [default: false]',
			category: 'Option: Filtering'
		}
		ploidy: {
			description: 'Expected ploidy of samples [default: 2]',
			category: 'Option: Filtering'
		}
		squashPloidy: {
			description: 'Treat heterozygous genotypes as homozygous ALT in both baseline and calls, to allow matches that ignore zygosity differences [default: false]',
			category: 'Option: Filtering'
		}
		precision: {
			description: 'Output summary statistics where precision >= supplied value (Default is to summarize at maximum F-measure)',
			category: 'Option: Reporting'
		}
		sensitivity: {
			description: 'Output summary statistics where sensitivity >= supplied value (Default is to summarize at maximum F-measure)',
			category: 'Option: Reporting'
		}
		roc: {
			description: 'Produce ROCs data',
			category: 'Option: Reporting'
		}
		outMode: {
			description: 'output reporting mode. Allowed values are [split, annotate, combine, ga4gh, roc-only] [Default: "split"]',
			category: 'Option: Reporting'
		}
		sortAsc: {
			description: 'Order in which to sort the ROC scores so that "good" scores come before "bad" scores. [deafault: false i.e descending]',
			category: 'Option: Reporting'
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
