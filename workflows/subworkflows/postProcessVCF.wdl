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

import "../../tasks/bcftools.wdl" as bcftools
import "../../tasks/GATK4.wdl" as GATK4
import "../../tasks/tabix.wdl" as tabix

workflow normalization {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-01-20"
	}

	input {
		File vcf

		File refFasta
		File refFai
		File refDict

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

		String outputPath = "./postProcessVCF/"
	}

################################################################################
## Post process

	call bcftools.norm as BNHC {
		input :
			in = vcf,
			splitMA = true,
			outputType = "z",
			refFasta = refFasta,
			outputPath = outputPath
	}

	call tabix.index as TBIBNHC {
		input :
			in = BNHC.outputFile,
			outputPath = outputPath
	}

	call GATK4.splitVcfs as GSpV {
		input :
			in = BNHC.outputFile,
			outputPath = outputPath
	}

	call GATK4.variantFiltration as GVF_SPN {
		input :
			in = GSpV.outputFileSnp,
			inIdx = GSpV.outputFileSnpIdx,
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
			outputPath = outputPath
	}

	call GATK4.variantFiltration as GVF_Indels {
		input :
			in = GSpV.outputFileIndel,
			inIdx = GSpV.outputFileIndelIdx,
			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,
			LowQualByDepth = LowQualByDepth,
			FSStrandBias = FSStrandBiasIndels,
			LowReadPosRankSum = LowReadPosRankSumIndels,
			SORStrandBias = SORStrandBiasIndels,
			HomopolymerRegion = HomopolymerRegion,
			LowCoverage = LowCoverage,
			outputPath = outputPath
	}

	call GATK4.mergeVcfs as GMV {
		input :
			in = [GVF_SPN.outputFile,GVF_Indels.outputFile],
			subString = "\.(snps|indels)\.filter\.(vcf|bcf|vcf\.gz)$",
			subStringReplace = ".filter",
			outputPath = outputPath
	}

	call GATK4.sortVcf as GSoV {
		input :
			in = GMV.outputFile,
			outputPath = outputPath
	}

	call bcftools.norm as BN {
		input :
			in = GSoV.outputFile,
			outputType = "z",
			multiallelicType = "any",
			refFasta = refFasta,
			outputPath = outputPath
	}

	call tabix.index as TI {
		input :
			in = BN.outputFile,
			outputPath = outputPath
	}

################################################################################

################################################################################

	output {
		File vcf = BN.outputFile
		File vcfIdx = TI.outputFile
	}

	parameter_meta {
		vcf : {
			description: 'Input file to process (vcf).',
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
		outputPath : {
			description: 'Path where the output repository will be created [default: ./]',
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
		outputPath : {
			description: 'Path where the output repository will be created [default: ./]',
			category: 'Output (optional)'
		}
	}
}
