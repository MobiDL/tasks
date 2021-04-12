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

import "../../tasks/bwa.wdl" as bwa
import "../../tasks/utilities.wdl" as utilities
import "../../tasks/GATK4.wdl" as GATK4

import "../subworkflows/postProcessAlignment.wdl" as postProcessAlignment

workflow ampliconsSoma {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-02-17"
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
		String subString = "(_S[0-9]+)?(_L[0-9]+)?(_R[12])?(_[0-9]+)?.(fastq|fq)(.gz)?"
		String? name

		## primary systems options
		Int SortThreads = 1
		Int MarkdupThreads = 1
		Int IndexThreads = 1
		Int ViewThreads = 1

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

	call utilities.convertBedToIntervals as Bed2Intervals {
		input :
			in = select_first([intervalBedFile,bedGenome.outputFile]),
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

################################################################################

################################################################################
## Alignement

	call bwa.mem as Alignment {
		input :
			fastqR1 = fastqR1,
			fastqR2 = fastqR2,

			subString = subString,
			sample = sampleName,

			refFasta = refFasta,
			refFai = refFai,
			refAmb = refAmb,
			refAnn = refAnn,
			refBwt = refBwt,
			refPac = refPac,
			refSa = refSa,

			outputPath = outputPath
	}

	call postProcessAlignment.postProcessAlignment as postProcessAlignment {
		input :
			inBam = Alignment.outputFile,

			splittedIntervals = SplitIntervals.splittedIntervals,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,

			knownSites = flatten([select_all(knownSites), [dbsnp]]),
			knownSitesIdx = flatten([select_all(knownSitesIdx), [dbsnpIdx]]),

			amplicons = true,

			name = sampleName,
			outputPath = outputPath + "/alignment/",

			SortThreads = SortThreads,
			MarkdupThreads = MarkdupThreads,
			IndexThreads = IndexThreads,
			ViewThreads = ViewThreads,

			SortMemory = SortMemory,
			MarkdupMemory = MarkdupMemory,
			IndexMemory = IndexMemory,
			ViewMemory = ViewMemory
	}

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
