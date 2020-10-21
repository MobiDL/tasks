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

workflow prepareGenome {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.2"
	}

	input {
		String genomeLink = "ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.fa.gz"
		String outputPath

		Boolean locale = false

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

	## Download and extract genome

	if (! locale) {
		call bash.wget as getGenome {
			input :
				in = genomeLink,
				outputPath = outputPath,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	if (locale) {
		call bash.makeLink as GenomeLn {
			input :
				in = genomeLink,
				outputPath = outputPath,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	File genome = select_first([GenomeLn.outputFile, getGenome.outputFile])

	Boolean isGzip = if sub(genomeLink, "(.*)(\.gz)$","$2") == ".gz" then true else false
	if (isGzip) {
		call bash.gzip as gzippedGenome {
			input :
				in = genome,
				outputPath = outputPath,
				decompress = true,
				memoryByThreads = memoryByThreadsLow,
				threads = minThreads
		}
	}
	File fastaRef = select_first([gzippedGenome.outputFile, genome])

	## Index and create dict from fasta
	call bwa.index as BwaIndexGenome {
		input :
			in = fastaRef,
			outputPath = outputPath,
			threads = maxThreads
	}

	call samtools.faidx as SamtoolsIndexGenome {
		input :
			in = fastaRef,
			outputPath = outputPath,
			threads = maxThreads
	}

	call samtools.dict as SamtoolsDictGenome {
		input :
			in = fastaRef,
			outputPath = outputPath,
			threads = maxThreads
	}

	output {
		File refFasta = fastaRef
		File refFai = SamtoolsIndexGenome.outputFile
		File refDict = SamtoolsDictGenome.outputFile
		File refAmb = BwaIndexGenome.refAmb
		File refAnn = BwaIndexGenome.refAnn
		File refBwt = BwaIndexGenome.refBwt
		File refPac = BwaIndexGenome.refPac
		File refSa = BwaIndexGenome.refSa
	}

	parameter_meta {
		genomeLink: {
			description : 'Link from genome to download [default: download GRCh37 from ensembl]',
			category : 'Optional'
		}
		outputPath: {
			description : 'Output path where fasta will be download and indexes created.',
			category : 'Required'
		}
		locale: {
			description : 'Defined if genome is locale (make hard link) or to download (wget) [default: false]',
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
