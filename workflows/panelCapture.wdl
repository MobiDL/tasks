version 1.0

#MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
#Copyright (C) 2020  MoBiDiC
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

import "prepareGenome.wdl" as prepareGenome
import "alignDNAcapture.wdl" as alignDNAcapture
import "variantCallingCaptureHC.wdl" as variantCallingCaptureHC
import "../tasks/tabix.wdl" as tabix


workflow panelCapture {
	meta {
		author: "Mobidic"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
	}

	input {
		# Input
		File fastqR1
		File fastqR2
		String fastqSubString = "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"
		String fastqSubStringReplace = ""
		String? name

		# Output
		String outputRep = "./"

		# Targets
		File intervalBedFile

		# Reference Genome
		File refFasta
		File? refFai
		File? refDict
		File? refAmb
		File? refAnn
		File? refBwt
		File? refPac
		File? refSa

		# Known sites
		File Indels1000G
		File IndelsMills
		File dbsnp
		File? Indels1000GIdx
		File? IndelsMillsIdx
		File? dbsnpIdx

		# Threshold for filters
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

		# System
		Int threads = 1
		Int maxThreads = threads
		Int minThreads = threads
	}

	String sampleName = if defined(name) then "~{name}" else sub(basename(fastqR1),fastqSubString,fastqSubStringReplace)
	String outputPath = "~{outputRep}/~{sampleName}/"

	if (! defined(refFai) || ! defined(refDict) || ! defined(refAmb) || ! defined(refAnn) || ! defined(refBwt) || ! defined(refPac) || ! defined(refSa)) {
		call prepareGenome.prepareGenome as PG {
			input :
				fasta = refFasta,
				outputPath = outputPath + "Genome/",

				threads = threads,
				maxThreads = maxThreads,
				minThreads = minThreads
		}
	}

	File refFaiUsed = if defined(PG.refFai) then "~{PG.refFai}" else "~{refFai}"
	File refDictUsed = if defined(PG.refDict) then "~{PG.refDict}" else "~{refDict}"
	File refAmbUsed = if defined(PG.refAmb) then "~{PG.refAmb}" else "~{refAmb}"
	File refAnnUsed = if defined(PG.refAnn) then "~{PG.refAnn}" else "~{refAnn}"
	File refBwtUsed = if defined(PG.refBwt) then "~{PG.refBwt}" else "~{refBwt}"
	File refPacUsed = if defined(PG.refPac) then "~{PG.refPac}" else "~{refPac}"
	File refSaUsed = if defined(PG.refSa) then "~{PG.refSa}" else "~{refSa}"

	if (! defined(Indels1000GIdx)) {
		call tabix.index as TI_Indels1000GIdx {
			input :
				in = Indels1000G,
				outputPath = outputPath + "KnownSitesIdx/"
		}
	}
	File Indels1000GIdxUsed = if defined(Indels1000GIdx) then "~{Indels1000GIdx}" else "~{TI_Indels1000GIdx.outputFile}"
	if (! defined(IndelsMillsIdx)) {
		call tabix.index as TI_IndelsMillsIdx {
			input :
				in = IndelsMills,
				outputPath = outputPath + "KnownSitesIdx/"
		}
	}
	File IndelsMillsIdxUsed = if defined(IndelsMillsIdx) then "~{IndelsMillsIdx}" else "~{TI_IndelsMillsIdx.outputFile}"
	if (! defined(dbsnpIdx)) {
		call tabix.index as TI_dbsnpIdx {
			input :
				in = dbsnp,
				outputPath = outputPath + "KnownSitesIdx/"
		}
	}
	File dbsnpIdxUsed = if defined(dbsnpIdx) then "~{dbsnpIdx}" else "~{TI_dbsnpIdx.outputFile}"

	call alignDNAcapture.alignDNAcapture as ADC {
		input :
			fastqR1 = fastqR1,
			fastqR2 = fastqR2,
			intervalBedFile = intervalBedFile,

			refFasta = refFasta,
			refFai = refFaiUsed,
			refDict = refDictUsed,
			refAmb = refAmbUsed,
			refAnn = refAnnUsed,
			refBwt = refBwtUsed,
			refPac = refPacUsed,
			refSa = refSaUsed,

			Indels1000G = Indels1000G,
			IndelsMills = IndelsMills,
			dbsnp = dbsnp,

			Indels1000GIdx = Indels1000GIdxUsed,
			IndelsMillsIdx = IndelsMillsIdxUsed,
			dbsnpIdx = dbsnpIdxUsed,

			outputRep = outputPath + "Alignement/",
			subString = fastqSubString,
			name = sampleName,

			threads = threads,
			maxThreads = maxThreads,
			minThreads = minThreads
	}

	call variantCallingCaptureHC.VCCaptureHC as VCHC {
		input :
			bam = ADC.cram,
			bamIdx = ADC.idx,
			intervalBedFile = intervalBedFile,
			refFasta = refFasta,
			refFai = refFaiUsed,
			refDict = refDictUsed,
			dbsnp = dbsnp,
			dbsnpIdx = dbsnpIdxUsed,
			outputRep =  outputPath + "VariantCallingHC/",
			name = sampleName,
			LowQualByDepth = LowQualByDepth,
			HomopolymerRegion = HomopolymerRegion,
			LowCoverage = LowCoverage,
			FSStrandBiasSNP = FSStrandBiasSNP,
			LowReadPosRankSumSNP = LowReadPosRankSumSNP,
			SORStrandBiasSNP = SORStrandBiasSNP,
			LowMappingQualitySNP = LowMappingQualitySNP,
			LowMappingQualityRankSumSNP = LowMappingQualityRankSumSNP,
			FSStrandBiasIndels = FSStrandBiasIndels,
			LowReadPosRankSumIndels = LowReadPosRankSumIndels,
			SORStrandBiasIndels = SORStrandBiasIndels,
			threads = threads,
			maxThreads = maxThreads,
			minThreads = minThreads
	}

	output {
		File cram = ADC.cram
		File idx = ADC.idx
		Array[File] fastQCreport = ADC.fastQCreport
		File bedcov = ADC.bedcov
		File flagstatReport = ADC.flagstatReport
		Array[File] qualimap = ADC.qualimap
		Array[File] collectMultipleMetrics = ADC.collectMultipleMetrics
		File vcf = VCHC.vcf
		File vcfIdx = VCHC.vcfIdx
		File vcfStats = VCHC.vcfStats

	}

}
