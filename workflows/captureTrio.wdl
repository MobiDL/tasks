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

import "panelCapture.wdl" as panelCapture
import "prepareGenome.wdl" as prepareGenome
import "../tasks/samtools.wdl" as samtools
import "../tasks/GATK4.wdl" as GATK4
import "../tasks/tabix.wdl" as tabix

workflow captureTrio {
  meta {
    author: "Mobidic"
    email: "o-ardouin<arobase>chu-montpellier.fr"
    version: "0.0.1"
  }

  input {
    # Input files
    # Fastq files
    File dadFastqR1
    File dadFastqR2
    File mumFastqR1
    File mumFastqR2
    File casFastqR1
    File casFastqR2

    # ref files
    File IntervalBedFile
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

		File dbsnp

    # Samples name
    String dadName
    String mumName
    String casName

    # Output
		String OutputRep = "./"

    # Find a way to specify threads and memoryByThreads
    # System
		Int threads = 1
		Int maxThreads = threads
		Int minThreads = threads
  }
  # define folders
  String OutputRepCas = "~{OutputRep}/~{casName}"
  String OutputRepDad = "~{OutputRep}/~{dadName}"
  String OutputRepMum = "~{OutputRep}/~{mumName}"
  String outputPath = OutputRep

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

  call panelCapture.panelCapture as panelCase {
    input :
      fastqR1 = casFastqR1,
      fastqR2 = casFastqR2,
      name = casName,
      outputRep = OutputRepCas,
      intervalBedFile = IntervalBedFile,
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
      dbsnpIdx = dbsnpIdxUsed
  }

  call panelCapture.panelCapture as panelDad {
    input :
      fastqR1 = dadFastqR1,
      fastqR2 = dadFastqR2,
      name = dadName,
      outputRep = OutputRepDad,
      intervalBedFile = IntervalBedFile,
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
      dbsnpIdx = dbsnpIdxUsed
  }

  call panelCapture.panelCapture as panelMum {
    input :
      fastqR1 = mumFastqR1,
      fastqR2 = mumFastqR2,
      name = mumName,
      outputRep = OutputRepDad,
      intervalBedFile = IntervalBedFile,
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
      dbsnpIdx = dbsnpIdxUsed
  }

  # merge CRAMs
  call samtools.merge {
    input :
      inputPaths = [panelCase.cram, panelDad.cram, panelMum.cram],
      outputPath = "~{OutputRep}/merge",
      name = "~{casName}_Trio",
  }

  # merge VCFs
  call GATK4.mergeVcfs {
    input :
      in = [panelCase.vcf, panelDad.vcf, panelMum.vcf],
      outputPath = "~{OutputRep}/merge",
      name = "~{casName}_Trio"

  }

}
