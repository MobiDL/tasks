version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2020 MoBiDiC
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

import "../tasks/bwa.wdl" as bwa
import "../tasks/bash.wdl" as bash
import "../tasks/samtools.wdl" as samtools
import "../tasks/star.wdl" as star

workflow prepareGenome {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.3"
	}

	input {
		String fastaLink = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz"
		String gtfLink = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz"
		String outputPath

		Boolean localeFasta = false
		Boolean localeGTF = false

		Boolean prepareStar = true

		Int memoryByThreads = 768
		Int memoryByThreadsHigh = memoryByThreads
		Int memoryByThreadsLow = memoryByThreads

		Int threads = 1
		Int maxThreads = threads
		Int minThreads = threads

		String memory = threads * memoryByThreads
		String memoryHigh = maxThreads * memoryByThreadsHigh
		String memoryLow = minThreads * memoryByThreadsLow
	}

################################################################################
####### STEP 1 : Fasta file

### 1.1 Download and extract fasta
	if (! localeFasta) {
		call bash.wget as fastaWget {
			input :
				in = fastaLink,
				outputPath = outputPath,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	if (localeFasta) {
		call bash.makeLink as fastaLn {
			input :
				in = fastaLink,
				outputPath = outputPath,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	File refFastaTemp = select_first([fastaWget.outputFile, fastaLn.outputFile])

### 1.2 Unzip if fasta is gzipped
	Boolean isFastaGzip = if sub(refFastaTemp, "(.*)(\.gz)$","$2") == ".gz" then true else false
	if (isFastaGzip) {
		call bash.gzip as gunzipFasta {
			input :
				in = refFastaTemp,
				outputPath = outputPath,
				decompress = true,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	File refFasta = select_first([gunzipFasta.outputFile, refFastaTemp])

### 1.3 Create indexes and dict for fasta
	call bwa.index as BwaIndexGenome {
		input :
			in = refFasta,
			outputPath = outputPath,
			threads = maxThreads
	}

	call samtools.faidx as SamtoolsIndexGenome {
		input :
			in = refFasta,
			outputPath = outputPath,
			threads = maxThreads
	}

	call samtools.dict as SamtoolsDictGenome {
		input :
			in = refFasta,
			outputPath = outputPath,
			threads = maxThreads
	}

####### END STEP 1
################################################################################

################################################################################
####### STEP 2 : Download annotations files

### 2.1 Download GTF
	if (! localeGTF) {
		call bash.wget as GTFWget {
			input :
				in = gtfLink,
				outputPath = outputPath,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	if (localeGTF) {
		call bash.makeLink as GTFLn {
			input :
				in = gtfLink,
				outputPath = outputPath,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	File refGTFTemp = select_first([GTFWget.outputFile, GTFLn.outputFile])

### 2.1 Unzip if GTF is gzipped
	Boolean isGTFGzip = if sub(refGTFTemp, "(.*)(\.gz)$","$2") == ".gz" then true else false
	if (isGTFGzip) {
		call bash.gzip as gunzipGTF {
			input :
				in = refGTFTemp,
				outputPath = outputPath,
				decompress = true,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	File refGTF = select_first([gunzipGTF.outputFile, refGTFTemp])

####### END STEP 2
################################################################################

################################################################################
####### STEP 3 : Prepare RNAseq files for STAR

### 3.1 Prepare Genome for RNAseq
	if (prepareStar) {
		call star.genomeGenerate as SGG {
			input :
				refFasta = refFasta,
				refGTF = refGTF,
				outputPath = outputPath,
				threads = maxThreads
		}
	}

####### END STEP 3
################################################################################

	output {
		File fasta = refFasta
		File refFai = SamtoolsIndexGenome.outputFile
		File refDict = SamtoolsDictGenome.outputFile
		File refAmb = BwaIndexGenome.refAmb
		File refAnn = BwaIndexGenome.refAnn
		File refBwt = BwaIndexGenome.refBwt
		File refPac = BwaIndexGenome.refPac
		File refSa = BwaIndexGenome.refSa
		File GTF = refGTF
		File? chrLength = SGG.chrLength
		File? chrNameLength = SGG.chrNameLength
		File? chrName = SGG.chrName
		File? chrStart = SGG.chrStart
		File? exonGeTrInfo = SGG.exonGeTrInfo
		File? exonInfo = SGG.exonInfo
		File? geneInfo = SGG.geneInfo
		File? Genome = SGG.Genome
		File? genomeParameters = SGG.genomeParameters
		File? SA = SGG.SA
		File? SAindex = SGG.SAindex
		File? sjdbInfo = SGG.sjdbInfo
		File? sjdbListGTF = SGG.sjdbListGTF
		File? sjdbList = SGG.sjdbList
		File? transcriptInfo = SGG.transcriptInfo
	}

	parameter_meta {
		fastaLink: {
			description : 'Link from genome to download [default: download GRCh37 from ensembl]',
			category : 'Optional'
		}
		outputPath: {
			description : 'Output path where fasta will be download and indexes created.',
			category : 'Required'
		}
		localeFasta: {
			description : 'Defined if fasta file is locale (make hard link) or to download (wget) [default: false]',
			category : 'Optional'
		}
		localeGTF: {
			description : 'Defined if GTF file is locale (make hard link) or to download (wget) [default: false]',
			category : 'Optional'
		}
		prepareStar: {
			description : 'Should launch star.genomeGenerate (needs at least 32G of RAM) [default: false]',
			category : 'Optional'
		}
		memoryByThreads : {
			description: 'Sets the memory by threads used by default [default: 768]',
			category: 'System'
		}
		memoryByThreadsHigh : {
			description: 'Sets the number of threads to use for high computing jobs [default: memoryByThreads]',
			category: 'System'
		}
		memoryByThreadsLow : {
			description: 'Sets the number of threads to use for low computing jobs [default: memoryByThreads]',
			category: 'System'
		}
		memory : {
			description: 'Sets the memory to use by default [default: memoryByThreads * threads]',
			category: 'System'
		}
		memoryHigh : {
			description: 'Sets the number of threads to use for high computing jobs [default: memoryByThreadsHigh * maxThreads]',
			category: 'System'
		}
		memoryLow : {
			description: 'Sets the number of threads to use for low computing jobs [default: memoryByThreadsLow * minThreads]',
			category: 'System'
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
