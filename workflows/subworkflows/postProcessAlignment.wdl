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

import "../../tasks/utilities.wdl" as utilities
import "../../tasks/crumble.wdl" as crumble
import "../../tasks/GATK4.wdl" as GATK4
import "../../tasks/sambamba.wdl" as sambamba

workflow postProcessAlignment {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.4"
		date: "2021-04-12"
	}

	input {
		File inBam

		Array[File] splittedIntervals

		File refFasta
		File refFai
		File refDict

		Array[File] knownSites
		Array[File] knownSitesIdx

		Boolean amplicons = false
		Boolean rna = false

		String outputPath = "./"
		String subString = "\.bam?"
		String? name

		Int SortThreads = 1
		Int MarkdupThreads = 1
		Int IndexThreads = 1
		Int ViewThreads = 1

		String SortMemory = "1G"
		String MarkdupMemory = "1G"
		String IndexMemory = "1G"
		String ViewMemory = "1G"
	}

	String sampleName = if defined(name) then "~{name}" else sub(basename(inBam),subString,"")

################################################################################

################################################################################
## Markduplicate or index if amplicons

	call sambamba.sort as SORT {
		input :
			in = inBam,
			outputPath = outputPath,

			threads = SortThreads,
			memory = SortMemory
	}

	if (! amplicons) {
		call sambamba.markdup as Markdup {
			input :
				in = SORT.outputBam,
				outputPath = outputPath,

				threads = MarkdupThreads,
				memory = MarkdupMemory
		}
	}

	if (rna) {
		call GATK4.splitNCigarReads as SPLITNCIGARREADS {
			input :
				in = select_first([Markdup.outputBam, SORT.outputBam]),

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "/Alignment/"
		}
	}

	File bamUsed = select_first([SPLITNCIGARREADS.outputBam, Markdup.outputBam, SORT.outputBam])
	File baiUsed = select_first([SPLITNCIGARREADS.outputIdx, Markdup.outputBai, SORT.outputBai])

################################################################################

################################################################################
## Base recalibrator & LeftAlign

	scatter (intervals in splittedIntervals) {
		call GATK4.baseRecalibrator as BR {
			input :
				in = bamUsed,
				bamIdx = baiUsed,

				intervals = intervals,

				knownSites = knownSites,
				knownSitesIdx = knownSitesIdx,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "/BQSR"
		}
	}

	call GATK4.gatherBQSRReports as GBR {
		input :
			in = BR.outputFile,
			outputPath = outputPath + "/BQSR"
	}

	scatter (intervals in splittedIntervals) {
		call GATK4.applyBQSR as ABR {
			input :
				in = bamUsed,
				bamIdx = baiUsed,
				bqsrReport = GBR.outputFile,
				intervals = intervals,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "/BQSR"
		}

		call GATK4.leftAlignIndels as LAI {
			input :
				in = ABR.outputBam,
				bamIdx = ABR.outputBai,
				intervals = intervals,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "/BQSR"
		}
	}

	call GATK4.gatherBamFiles {
		input :
			in = LAI.outputBam,
			bamIdx = LAI.outputBai,

			outputPath = outputPath
	}

################################################################################

################################################################################
## Postprocess (sort/index/compression)

	call sambamba.sort as SortBamProcessed {
		input :
			in = gatherBamFiles.outputBam,
			outputPath = outputPath,

			threads = SortThreads,
			memory = SortMemory
	}

	call sambamba.view as Bam2Cram{
		input :
			in = SortBamProcessed.outputBam,
			refFasta = refFasta,
			refFai = refFai,
			cram = true,
			outputPath = outputPath,

			threads = ViewThreads,
			memory = ViewMemory
	}

	call sambamba.index as CramIdx {
		input :
			in = Bam2Cram.outputFile,
			outputPath = outputPath,

			threads = IndexThreads,
			memory = IndexMemory
	}

	## TODO: fix crumble
	/* call crumble.crumble as Crumble {
		input :
			in = Bam2Cram.outputFile,
			outputPath = outputPath
	}

	call sambamba.index as CrumbleIdx {
		input :
			in = Crumble.outputFile,
			outputPath = outputPath
	} */

################################################################################

################################################################################

	output {
		File cram = Bam2Cram.outputFile
		File crai = CramIdx.outputFile
		File bam = SortBamProcessed.outputBam
		File bai = SortBamProcessed.outputBai
	}

################################################################################

	parameter_meta {
		inBam: {
			description: 'Input BAM.',
			category: 'Input (required)'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'Input (required)'
		}
		refFai: {
			description: 'Path to the reference file index (format: fai)',
			category: 'Input (required)'
		}
		refDict: {
			description: 'Path to the reference file dict (format: dict)',
			category: 'Input (required)'
		}
		knownSites : {
			description: 'Path to the knownSites vcf files (format: vcf.gz)',
			category: 'Input (required)'
		}
		knownSitesIdx : {
			description: 'Path to the knownSites vcf indexes (format: vcf.gz.tbi)',
			category: 'Input (required)'
		}
		amplicons : {
			description: 'Mark duplicates if set to false [default: false]',
			category: 'Option'
		}
		outputPath : {
			description: 'Path where the output repository will be created [default: ./]',
			category: 'Output (optional)'
		}
		subString : {
			description: 'The regexp substring to remove from fastq R1 file to create sampleName, and used by bwa-mem [default: "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"]',
			category: 'Output (optional)'
		}
		name : {
			description: 'The name used as sampleName',
			category: 'Output (optional)'
		}
		SortThreads : {
			description : 'Number of threads for "Sort" step [default : 1]',
			category : 'System'
		}
		MarkdupThreads : {
			description : 'Number of threads for "Markdup" step [default : 1]',
			category : 'System'
		}
		IndexThreads : {
			description : 'Number of threads for "Index" step [default : 1]',
			category : 'System'
		}
		ViewThreads : {
			description : 'Number of threads for "View" step [default : 1]',
			category : 'System'
		}
		SortMemory : {
			description : 'Memory allocated for "Sort" step [default : "1G"]',
			category : 'System'
		}
		MarkdupMemory : {
			description : 'Memory allocated for "Markdup" step [default : "1G"]',
			category : 'System'
		}
		IndexMemory : {
			description : 'Memory allocated for "Index" step [default : "1G"]',
			category : 'System'
		}
		ViewMemory : {
			description : 'Memory allocated for "View" step [default : "1G"]',
			category : 'System'
		}
	}
}
