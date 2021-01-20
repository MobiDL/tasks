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

import "../tasks/bash.wdl" as bash
import "../tasks/bcftools.wdl" as bcftools
import "../tasks/GATK4.wdl" as GATK4
import "../tasks/tabix.wdl" as tabix

workflow VCCaptureHC {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
	}

	input {
		File bam
		File bamIdx

		File intervalBedFile

		File refFasta
		File refFai
		File refDict

		File? dbsnp
		File? dbsnpIdx

		String outputPath = "./"
		String subString = "\.(sam|bam|cram)"
		String? name

		Int LowQualByDepth = 2
		Int HomopolymerRegion = 7
		Int LowCoverage = 20

		Int FSStrandBiasSNP = 60
		Int LowReadPosRankSumSNP = -4
		Int SORStrandBiasSNP = 3
		Int LowMappingQualitySNP = 40
		Int LowMappingQualityRankSumSNP = -3

		Int FSStrandBiasIndels = 200
		Int LowReadPosRankSumIndels = -5
		Int SORStrandBiasIndels = 10

		Int threads = 1
		Int maxThreads = threads
		Int minThreads = threads
	}

	String sampleName = if defined(name) then "~{name}" else sub(basename(bam),subString,"")

	call bash.convertBedToIntervals as BB2I {
		input :
			in = intervalBedFile,
			outputPath = outputPath + "regionOfInterest/"
	}

	call GATK4.splitIntervals as GSI {
		input :
			in = BB2I.outputFile,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,

			subdivisionMode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",

			outputPath = outputPath + "regionOfInterest/",

			scatterCount = maxThreads,
			threads = minThreads
	}

	scatter (intervals in GSI.splittedIntervals) {
		call GATK4.haplotypeCaller as GHC {
			input :
				in = bam,
				idx = bamIdx,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				name = sampleName,

				dbsnp = dbsnp,
				dbsnpIdx = dbsnpIdx,

				intervals = intervals,

				outputPath = outputPath + "variant-calling/split/"
		}
	}

	call GATK4.gatherVcfFiles as GGV{
		input :
			in = GHC.outputFile,
			outputPath = outputPath + "variant-calling/"
	}

	call GATK4.splitVcfs as GSpV {
		input :
			in = GGV.outputFile,
			outputPath = outputPath + "variant-calling/"
	}

	call GATK4.variantFiltration as GVF_SPN {
		input :
			in = GSpV.outputFileSnp,
			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,
			LowQualByDepth = LowQualByDepth,
			FSStrandBias = FSStrandBiasSNP,
			LowReadPosRankSum = LowReadPosRankSumSNP,
			SORStrandBias = SORStrandBiasSNP,
			HomopolymerRegion = HomopolymerRegion,
			LowCoverage = LowCoverage,
			LowMappingQuality = LowMappingQualitySNP,
			LowMappingQualityRankSum = LowMappingQualityRankSumSNP,
			outputPath = outputPath + "variant-calling/"
	}

	call GATK4.variantFiltration as GVF_Indels {
		input :
			in = GSpV.outputFileIndel,
			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,
			LowQualByDepth = LowQualByDepth,
			FSStrandBias = FSStrandBiasIndels,
			LowReadPosRankSum = LowReadPosRankSumIndels,
			SORStrandBias = SORStrandBiasIndels,
			HomopolymerRegion = HomopolymerRegion,
			LowCoverage = LowCoverage,
			outputPath = outputPath + "variant-calling/"
	}

	call GATK4.mergeVcfs as GMV {
		input :
			in = [GVF_SPN.outputFile,GVF_Indels.outputFile],
			subString = "(snps|indels)\.filter\.(vcf|bcf)$",
			subStringReplace = "filter",
			outputPath = outputPath + "variant-calling/"
	}

	call GATK4.sortVcf as GSoV {
		input :
			in = GMV.outputFile,
			outputPath = outputPath + "variant-calling/"
	}

	call bcftools.norm as BN {
		input :
			in = GSoV.outputFile,
			splitMA = true,
			outputType = "z",
			refFasta = refFasta,
			outputPath = outputPath + "variant-calling/"
	}

	call tabix.index as TI {
		input :
			in = BN.outputFile,
			outputPath = outputPath + "variant-calling/"
	}

	call bcftools.stats as BS {
		input :
			in = BN.outputFile,
			refFasta = refFasta,
			refFai = refFai,
			outputPath = outputPath + "variant-calling/quality/"
	}

	output {
		File vcf = TI.outputFile
		File vcfIdx = BS.outputFile
		File vcfStats = BS.outputFile
	}

	parameter_meta {
		bam : {
			description: 'Input file to process variant calling (bam, cram).',
			category: 'Required'
		}
		bamIdx : {
			description: 'Index of the input file (bai, crai).',
			category: 'Required'
		}
		intervalBedFile : {
			description: 'Input bed of the design.',
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
		dbsnp : {
			description: 'Path to the dbsnp vcf file (format: vcf.gz)',
			category: 'Required'
		}
		dbsnpIdx : {
			description: 'Path to the dbsnp vcf index file (format: vcf.gz.tbi)',
			category: 'Required'
		}
		outputPath : {
			description: 'Path where the output repository will be created [default: ./]',
			category: 'Option'
		}
		subString : {
			description: 'The regexp substring to remove from the input file to create sampleName [default: "\.(sam|bam|cram)"]',
			category: 'Option'
		}
		name : {
			description: 'The name used as sampleName (at the place of subString)',
			category: 'Option'
		}
		LowQualByDepth: {
			description : 'Threshold below which QD the variant will be tagged as LowQualByDepth',
			category: 'Option'
		}
		HomopolymerRegion: {
			description : 'Threshold above which POLYX the variant will be tagged as HomopolymerRegion',
			category: 'Option'
		}
		LowCoverage: {
			description : 'Threshold below which DP the variant will be tagged as LowCoverage',
			category: 'Option'
		}
		FSStrandBiasSNP: {
			description : 'Threshold above which FS the variant will be tagged as FSStrandBias for SNPs',
			category: 'Option'
		}
		LowReadPosRankSumSNP: {
			description : 'Threshold below which ReadPosRankSum the variant will be tagged as LowReadPosRankSum for SNPs',
			category: 'Option'
		}
		SORStrandBiasSNP: {
			description : 'Threshold above which SOR the variant will be tagged as SORStrandBias for SNPs',
			category: 'Option'
		}
		LowMappingQualitySNP: {
			description : 'Threshold below which MQ the variant will be tagged as LowMappingQuality for SNPs',
			category: 'Option'
		}
		LowMappingQualityRankSumSNP: {
			description : 'Threshold below which MQRankSum the variant will be tagged as LowMappingQualityRankSum for SNPs',
			category: 'Option'
		}
		FSStrandBiasIndels: {
			description : 'Threshold above which FS the variant will be tagged as FSStrandBias for INDELs',
			category: 'Option'
		}
		LowReadPosRankSumIndels: {
			description : 'Threshold below which ReadPosRankSum the variant will be tagged as LowReadPosRankSum for INDELs',
			category: 'Option'
		}
		SORStrandBiasIndels: {
			description : 'Threshold above which SOR the variant will be tagged as SORStrandBias for INDELs',
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
