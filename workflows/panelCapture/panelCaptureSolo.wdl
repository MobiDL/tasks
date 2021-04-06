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
import "../../tasks/GATK4.wdl" as GATK4
import "../../tasks/vep.wdl" as vep

import "../subworkflows/alignmentShortReadsDNA.wdl" as alignment
import "../subworkflows/postProcessVCF.wdl" as postProcessVCF

workflow panelCaptureSolo {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.3"
		date: "2021-02-25"
	}

	input {
		File fastqR1
		File? fastqR2

		File? intervalBedFile

		File refFasta
		File refFai
		File refDict
		File refAmb
		File refAnn
		File refBwt
		File refPac
		File refSa

		Array[File?] knownSites
		Array[File?] knownSitesIdx

		File dbsnp
		File dbsnpIdx

		Int scatterCount=24

		String outputPath = "./"
		String subString = "(_S[0-9]+)?(_L[0-9]+)?([_.]?R[12])?(_[0-9]+)?.(fastq|fq)(.gz)?"
		String? name

		## primary systems options
		Int BWAThreads = 1
		Int SortThreads = 1
		Int MarkdupThreads = 1
		Int IndexThreads = 1
		Int ViewThreads = 1

		String BWAMemory = "4G"
		String SortMemory = "1G"
		String MarkdupMemory = "4G"
		String IndexMemory = "1G"
		String ViewMemory = "1G"
	}

	String sampleName = if defined(name) then "~{name}" else sub(basename(fastqR1),subString,"")

################################################################################
## Preprocessing

	call utilities.fai2bed as bedGenome {
		input :
			in = refFai
	}

	File bed2use = select_first([intervalBedFile,bedGenome.outputFile])

	call utilities.convertBedToIntervals as Bed2Intervals {
		input :
			in = bed2use,
			outputPath = outputPath + "/regionOfInterest/"
	}

	call GATK4.splitIntervals as SplitIntervals {
		input :
			in = Bed2Intervals.outputFile,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,

			subdivisionMode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",

			outputPath = outputPath + "/regionOfInterest/",

			scatterCount = scatterCount
	}

	scatter (intervals in SplitIntervals.splittedIntervals) {
		call GATK4.intervalListToBed as IL2B {
			input :
				in = intervals,
				outputPath = outputPath + "regionOfInterest/split-bed/"
		}
	}

################################################################################

################################################################################
## Alignement

	call alignment.alignDNA as primary {
		input :
			fastqR1 = fastqR1,
			fastqR2 = fastqR2,

			splittedIntervals = SplitIntervals.splittedIntervals,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,
			refAmb = refAmb,
			refAnn = refAnn,
			refBwt = refBwt,
			refPac = refPac,
			refSa = refSa,

			knownSites = flatten([select_all(knownSites), [dbsnp]]),
			knownSitesIdx = flatten([select_all(knownSitesIdx), [dbsnpIdx]]),

			name = sampleName,
			outputPath = outputPath + "/alignment/",

			BWAThreads = BWAThreads,
			SortThreads = SortThreads,
			MarkdupThreads = MarkdupThreads,
			IndexThreads = IndexThreads,
			ViewThreads = ViewThreads,

			BWAMemory = BWAMemory,
			SortMemory = SortMemory,
			MarkdupMemory = MarkdupMemory,
			IndexMemory = IndexMemory,
			ViewMemory = ViewMemory
	}

################################################################################

################################################################################
## Variant Calling

	### HaplotypeCaller
	scatter (intervals in SplitIntervals.splittedIntervals) {
		call GATK4.haplotypeCaller as GHC {
			input :
				in = primary.cram,
				idx = primary.idx,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				name = sampleName,

				dbsnp = dbsnp,
				dbsnpIdx = dbsnpIdx,

				intervals = intervals,

				outputPath = outputPath + "variant-calling/HaplotypeCaller/"
		}
	}

	call GATK4.gatherVcfFiles as GGV {
		input :
			in = GHC.outputFile,
			outputPath = outputPath + "variant-calling/HaplotypeCaller/"
	}

	call postProcessVCF.normalization as normVCF {
		input :
			vcf = GGV.outputFile,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,

			outputPath = outputPath + "variant-calling/HaplotypeCaller/"
	}
	###

################################################################################

################################################################################
## Annotation

	### VEP
	call vep.vep_cache as annotVEP {
		input :
			in = normVCF.vcfOut,
			outputPath = outputPath + "variant-annot/",
			name = sampleName,
			refFasta = refFasta,
			refFai = refFai
	}
	###

################################################################################

################################################################################

	parameter_meta {
		fastqR1 : {
			description: 'Input file with reads 1 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		fastqR2 : {
			description: 'Input file with reads 2 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		intervalBedFile : {
			description: 'Input bed to convert to intervals (chr:start-end).',
			category: 'Required'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'Required'
		}
		refFai: {
			description: 'Path to the reference file index (format: fai)',
			category: 'Required'
		}
		refDict: {
			description: 'Path to the reference file dict (format: dict)',
			category: 'Required'
		}
		refAmb : {
			description: 'Path to the reference Amb file (generate by BWA index)',
			category: 'Required'
		}
		refAnn : {
			description: 'Path to the reference Ann file (generate by BWA index)',
			category: 'Required'
		}
		refBwt : {
			description: 'Path to the reference Bwt file (generate by BWA index)',
			category: 'Required'
		}
		refPac : {
			description: 'Path to the reference Pac file (generate by BWA index)',
			category: 'Required'
		}
		refSa : {
			description: 'Path to the reference Sa file (generate by BWA index)',
			category: 'Required'
		}
		knownSites : {
			description: 'Path to the knownSites vcf files (format: vcf.gz)',
			category: 'Required'
		}
		knownSitesIdx : {
			description: 'Path to the knownSites vcf indexes (format: vcf.gz.tbi)',
			category: 'Required'
		}
		dbsnp : {
			description: 'Path to the dbsnp vcf file (format: vcf.gz)',
			category: 'Required'
		}
		dbsnpIdx : {
			description: 'Path to the dbsnp vcf index (format: vcf.gz.tbi)',
			category: 'Required'
		}
		outputPath : {
			description: 'Path where the output repository will be created [default: ./]',
			category: 'Option'
		}
		subString : {
			description: 'The regexp substring to remove from fastq R1 file to create sampleName, and used by bwa-mem [default: "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"]',
			category: 'Option'
		}
		name : {
			description: 'The name used as sampleName',
			category: 'Option'
		}
		BWAThreads : {
			description : 'Number of threads for "BWA" step [default : 1]',
			category : 'System'
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
		BWAMemory : {
			description : 'Memory allocated for "BWA" step [default : "4G"]',
			category : 'System'
		}
		SortMemory : {
			description : 'Memory allocated for "Sort" step [default : "1G"]',
			category : 'System'
		}
		MarkdupMemory : {
			description : 'Memory allocated for "Markdup" step [default : "4G"]',
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
