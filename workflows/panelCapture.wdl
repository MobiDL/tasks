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

			outputPath = outputPath + "Alignement/",
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

	parameter_meta {
		fastqR1 : {
			description: 'Input file with reads 1 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		fastqR2 : {
			description: 'Input file with reads 2 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		fastqSubString : {
			description: 'The regexp substring to remove from fastq R1 file to create sampleName, and used by bwa-mem [default: "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"]',
			category: 'Option'
		}
		fastqSubStringReplace : {
			description: 'A string to replace regex catch by fastqSubString to create sampleName [default: ""]',
			category: 'Option'
		}
		name : {
			description: 'The name used as sampleName',
			category: 'Option'
		}
		outputRep: {
			description: 'Output path of the repertory where files will be generated. [default: pwd()]',
			category: 'Option'
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
			category: 'Option'
		}
		refDict: {
			description: 'Path to the reference file dict (format: dict)',
			category: 'Option'
		}
		refAmb : {
			description: 'Path to the reference Amb file (generate by BWA index)',
			category: 'Option'
		}
		refAnn : {
			description: 'Path to the reference Ann file (generate by BWA index)',
			category: 'Option'
		}
		refBwt : {
			description: 'Path to the reference Bwt file (generate by BWA index)',
			category: 'Option'
		}
		refPac : {
			description: 'Path to the reference Pac file (generate by BWA index)',
			category: 'Option'
		}
		refSa : {
			description: 'Path to the reference Sa file (generate by BWA index)',
			category: 'Option'
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
			category: 'Option'
		}
		IndelsMillsIdx : {
			description: 'Path to the IndelsMills vcf index file (format: vcf.gz.tbi)',
			category: 'Option'
		}
		dbsnpIdx : {
			description: 'Path to the dbsnp vcf index file (format: vcf.gz.tbi)',
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
