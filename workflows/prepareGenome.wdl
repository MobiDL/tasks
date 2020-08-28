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

import "../tasks/bwa.wdl" as bwa
import "../tasks/bash.wdl" as bash
import "../tasks/samtools.wdl" as samtools

workflow prepareGenome {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
	}

	input {
		File fasta
		String outputPath

		Int threads = 1
		Int maxThreads = threads
		Int minThreads = threads
	}

	call bash.makeLink as GenomeLn {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = minThreads
	}

	call bwa.index as BwaIndexGenome {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = maxThreads
	}

	call samtools.faidx as SamtoolsIndexGenome {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = maxThreads
	}

	call samtools.dict as SamtoolsDictGenome {
		input :
			in = fasta,
			outputPath = outputPath,
			threads = maxThreads
	}

	output {
		File refFasta = GenomeLn.outputFile
		File refFai = SamtoolsIndexGenome.outputFile
		File refDict = SamtoolsDictGenome.outputFile
		File refAmb = BwaIndexGenome.refAmb
		File refAnn = BwaIndexGenome.refAnn
		File refBwt = BwaIndexGenome.refBwt
		File refPac = BwaIndexGenome.refPac
		File refSa = BwaIndexGenome.refSa
	}

	parameter_meta {
		fasta: {
			description : 'Path to the reference file (format: fasta)',
			category : 'Required'
		}
		outputPath: {
			description : 'Output path where bam file was generated.',
			category : 'Required'
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
