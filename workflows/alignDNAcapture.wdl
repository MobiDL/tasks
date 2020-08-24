version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2020  MoBiDiC
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import "../tasks/bash.wdl" as bash
import "../tasks/bwa.wdl" as bwa
import "../tasks/crumble.wdl" as crumble
import "../tasks/fastqc.wdl" as fastqc
import "../tasks/GATK4.wdl" as GATK4
import "../tasks/sambamba.wdl" as sambamba
import "../tasks/samtools.wdl" as samtools
import "../tasks/qualimap.wdl" as qualimap

workflow alignDNAcapture {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
	}

	input {
		File fastqR1
		File fastqR2

		File intervalBedFile

		File refFasta
		File refFai
		File refDict
		File refAmb
		File refAnn
		File refBwt
		File refPac
		File refSa

		File Indels1000G
		File IndelsMills
		File dbsnp
		File Indels1000GIdx
		File IndelsMillsIdx
		File dbsnpIdx

		String outputRep = "./"
		String subString = "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"
		String? name

		Int threads = 1
		Int maxThreads = threads
		Int minThreads = threads
	}

	String sampleName = sub(basename(fastqR1),subString,"")
	String outputPath = if defined(name) then "~{outputRep}/~{name}/" else"~{outputRep}/~{sampleName}/"

	## Preprocessing

	call bash.convertBedToIntervals as Bed2Intervals {
		input :
			in = intervalBedFile,
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

			scatterCount = maxThreads,
			threads = minThreads
	}

	scatter (intervals in SplitIntervals.splittedIntervals) {
		call GATK4.intervalListToBed as IL2B {
			input :
				in = intervals,
				outputPath = outputPath + "regionOfInterest/split-bed/"
		}
	}

	## Alignement

	call bwa.mem as Alignement {
		input :
			fastqR1 = fastqR1,
			fastqR2 = fastqR2,
			subString = subString,

			refFasta = refFasta,
			refFai = refFai,
			refAmb = refAmb,
			refAnn = refAnn,
			refBwt = refBwt,
			refPac = refPac,
			refSa = refSa,

			outputPath = outputPath + "/alignement/",

			threads = maxThreads
	}

	call sambamba.markdup as Markdup {
		input :
			in = Alignement.outputFile,
			outputPath = outputPath + "/alignement/"
	}

	scatter (intervals in SplitIntervals.splittedIntervals) {
		call GATK4.baseRecalibrator as BR {
			input :
				in = Markdup.outputBam,
				bamIdx = Markdup.outputBai,

				intervals = intervals,

				knownSites = [Indels1000G, IndelsMills, dbsnp],
				knownSitesIdx = [Indels1000GIdx, IndelsMillsIdx, dbsnpIdx],

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "alignement/BQSR",

				threads = threads
		}
	}

	call GATK4.gatherBQSRReports as GBR {
		input :
			in = BR.outputFile,
			outputPath = outputPath + "alignement/BQSR"
	}

	scatter (intervals in SplitIntervals.splittedIntervals) {
		call GATK4.applyBQSR as ABR {
			input :
				in = Markdup.outputBam,
				bamIdx = Markdup.outputBai,
				bqsrReport = GBR.outputFile,
				intervals = intervals,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "alignement/BQSR"
		}

		call GATK4.leftAlignIndels as LAI {
			input :
				in = ABR.outputBam,
				bamIdx = ABR.outputBai,
				intervals = intervals,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "alignement/BQSR"
		}
	}

	call GATK4.gatherBamFiles {
		input :
			in = LAI.outputBam,
			bamIdx = LAI.outputBai,

			outputPath = outputPath + "alignement/"
	}

	call samtools.sort as SortBamProcessed {
		input :
			in = gatherBamFiles.outputBam,
			outputPath = outputPath + "alignement/",

			threads = maxThreads
	}

	call sambamba.index as IdxBamGatherSort {
		input :
			in = SortBamProcessed.outputFile,
			outputPath = outputPath + "alignement/",
			threads = maxThreads
	}

	call samtools.view as Bam2Cram{
		input :
			in = SortBamProcessed.outputFile,
			refFasta = refFasta,
			faidx = refFai,
			cram = true,
			outputPath = outputPath + "alignement/",
			threads = minThreads
	}

	call samtools.index as CramIdx {
		input :
			in = Bam2Cram.outputFile,
			outputPath = outputPath + "alignement/",
			threads = minThreads
	}

	call crumble.crumble as Crumble {
		input :
			in = Bam2Cram.outputFile,
			outputPath = outputPath + "alignement/",
			threads = minThreads
	}

	call samtools.index as CrumbleIdx {
		input :
			in = Crumble.outputFile,
			outputPath = outputPath + "alignement/",
			threads = minThreads
	}

	## Quality

	### Raw reads

	call fastqc.fastqc as FastQC {
		input :
			in = [fastqR1, fastqR2],
			outputPath = outputPath + "/quality/"
	}

	### Alignement

	scatter (bed in IL2B.outputFile) {
		call samtools.bedcov as BC {
			input :
				inBam = Crumble.outputFile,
				inBamIdx = CrumbleIdx.outputFile,
				qualityThreshold = 30,
				inBed = bed,
				name = sampleName + "_" + sub(basename(bed),".bed",""),
				outputPath = outputPath + "/quality/scatter/"
		}
	}

	call bash.concatenateFiles as BCT {
		input :
			in = BC.outputFile,
			subString = sampleName + "_[0-9]+",
			outputPath = outputPath + "/quality/"
	}

	call samtools.flagstat as FS {
		input :
			in = Crumble.outputFile,
			outputPath = outputPath + "/quality/"
	}

	call GATK4.collectMultipleMetrics as CMM {
		input :
			in = Crumble.outputFile,
			bamIdx = CrumbleIdx.outputFile,
			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,
			outputPath = outputPath + "/quality/GATK4Metrics/"
	}

	call qualimap.bamqc as QM {
		input :
			in = SortBamProcessed.outputFile,
			featureFile = intervalBedFile,
			threads = 4,
			outsideStats = false,
			outputPath = outputPath + "/quality/qualimap/"
	}

	output {
		File cram = Crumble.outputFile
		File idx = CrumbleIdx.outputFile
		Array[File] fastQCreport = FastQC.outHTML
		File bedcov = BCT.outputFile
		File flagstatReport = FS.outputFile
		Array[File] qualimap = QM.reports
		Array[File] collectMultipleMetrics = CMM.collectMultipleMetrics
	}

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
			category: 'required'
		}
		refFai: {
			description: 'Path to the reference file index (format: fai)',
			category: 'required'
		}
		refDict: {
			description: 'Path to the reference file dict (format: dict)',
			category: 'required'
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
		Indels1000G : {
			description: 'Path to the Indels1000G vcf file (format: vcf.gz)',
			category: 'Required'
		}
		IndelsMills : {
			description: 'Path to the IndelsMills vcf file (format: vcf.gz)',
			category: 'Required'
		}
		dbsnp : {
			description: 'Path to the dbsnp vcf file (format: vcf.gz)',
			category: 'Required'
		}
		Indels1000GIdx : {
			description: 'Path to the Indels1000G vcf index file (format: vcf.gz.tbi)',
			category: 'Required'
		}
		IndelsMillsIdx : {
			description: 'Path to the IndelsMills vcf index file (format: vcf.gz.tbi)',
			category: 'Required'
		}
		dbsnpIdx : {
			description: 'Path to the dbsnp vcf index file (format: vcf.gz.tbi)',
			category: 'Required'
		}
		outputRep : {
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
		threads : {
			description: 'Sets the number of threads to use by default [default: 1]',
			category: 'System'
		}
		maxThreads : {
			description: 'Sets the number of threads to use for high computing jobs [default: threads]',
			category: 'System'
		}
		minThreads : {
			description: 'Sets the number of threads to use for low computing jobs [default: threads]',
			category: 'System'
		}
	}
}
