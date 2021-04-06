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

import "../../tasks/nanovar.wdl" as nanovar
import "../../tasks/utilities.wdl" as utlities
import "../../tasks/minimap2.wdl" as minimap2
import "../../tasks/sambamba.wdl" as sambamba
import "../../tasks/clair.wdl" as clair

workflow allVariantCalling {
	meta {
		author: "Thomas GUIGNARD"
		email: "t-guignard(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-19"
	}
	input {

		Boolean SVcalling = true
		Boolean SNVcalling = true

		String fastqPath = "/mnt/RS_IURC/MinION/DPNNI/Analyses_bioinfo_Anaïs/202103160950_basecalling_703-DIA1-EXP/Guppy_703-DIA1-EXP/"
		String outputRep = "/home/bioinfo/MobiDL2.0-new/TEST/"

		File refFa = "/mnt/RS_IURC/MobiDL/Datasets/genomes/hg19/hg19.fa"
		File refFai = "/mnt/RS_IURC/MobiDL/Datasets/genomes/hg19/hg19.fa.fai"

		String modelPath = "/home/bioinfo/clair/ont/"

		Int qual = 748
		String? name


		File InputNanovar
		File RefFasta
		File RefGenomeIndex
		String OutputPath
		Float Score
		String GenomeBuild
		Int Threads
		String GapGenomeBuild

	}
	call utilities.findFiles as FF {
		input :
			path = fastqPath,
			regexpName = "*.fastq",
			maxDepth = 1
	}
	call utilities.concatenateFiles as ConcFQ {
		input :
			in = FF.files,
			name = sampleName + ".fastq",
			outputPath = outputRep + "/fastq_concatenate/"
	}
	call minimap2.mapOnt as align {
		input :
			fastq = ConcFQ.outputFile,
			refFasta = refFa,
			sample = sampleName,
			outputPath = outputRep + "/Alignment/"
	}
	call sambamba.index as idxBam {
		input :
			in = align.outputFile,
			outputPath = outputRep + "/Alignment/"
	}

	if (SVcalling == true){

		call nanovar.nanovar {
			input :
				inputNanovar = InputNanovar,
				refFasta = RefFasta,
				outputPath = OutputPath,
				score = Score,
				gapGenomeBuild = GapGenomeBuild,
				threads = Threads,
				genomeBuild = GenomeBuild,
				refGenomeIndex = RefGenomeIndex
		}
		call utilities.nanoVar2Bed {
			input :


		}

	}
	if (SNVcalling == true){

		call clair.callVarBam {
			input :
				modelPath = modelPath,
				refGenome = refFa,
				refGenomeIndex = refFai,
				bamFile = align.outputFile,
				bamFileIndex = idxBam.outputFile,
				sampleName = sampleName,
				qual = qual,
				outputPath = outputRep + "/Variant-Calling/"
			}
			/* call bgzip.compress{
				input :

			}
			call tabix.index{
				input :
			}
			call bcftools.view {
				input :
				bcftools
				} */
			}



}
