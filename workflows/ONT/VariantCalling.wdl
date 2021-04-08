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

import "../../tasks/utilities.wdl" as utilities
import "../../tasks/minimap2.wdl" as minimap2
import "../../tasks/sambamba.wdl" as sambamba
import "../../tasks/clair.wdl" as clair
import "../../tasks/GATK4.wdl" as GATK4

workflow variantCallingONT {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-04-02"
	}

	input {
		String fastqPath
		String outputPath

		File refFa
		File refFai

		String modelPath

		Int qual = 748
		File? bedRegions
		String? includeFilter

		String? name
	}

	String sampleName = if defined(name) then "~{name}" else "sample"

################################################################################
## Alignment

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
			outputPath = outputPath + "/fastq_concatenate/"
	}

	call minimap2.mapOnt as align {
		input :
			fastq = ConcFQ.outputFile,
			refFasta = refFa,
			sample = sampleName,
			outputPath = outputPath + "/Alignment/"
	}

	call sambamba.index as idxBam {
		input :
			in = align.outputFile,
			outputPath = outputPath + "/Alignment/"
	}

################################################################################

################################################################################
## Variant Calling

	call utilities.fai2bed as F2B {
		input :
			in = refFai
	}

	call utilities.bed2Array as B2A {
		input :
			bed = select_first([bedRegions,F2B.outputFile])
	}

	scatter (region in B2A.bedObj) {
		call clair.callVarBam as CALLVARBAM {
			input :
				modelPath = modelPath,
				refGenome = refFa,
				refGenomeIndex = refFai,
				bamFile = align.outputFile,
				bamFileIndex = idxBam.outputFile,
				name = "~{sampleName}-~{region['chrom']}_~{region['start']}_~{region['end']}",
				sampleName = sampleName,
				qual = qual,
				outputPath = "~{outputPath}/Variant-Calling/clair/scatter/",
				contigName = region["chrom"],
				ctgStart = region["start"],
				ctgEnd = region["end"],
				delay = 10
		}
	}

	call GATK4.gatherVcfFiles as GATHERVCFFILES {
		input :
			in = CALLVARBAM.outputFile,
			outputPath = outputPath + "/Variant-Calling/clair/",
			subString = "-(chr)?[0-9WXYMT]+_[0-9]+_[0-9]+\.clair\.vcf$"
	}

################################################################################

################################################################################
## Outputs

	output {
		File bam = align.outputFile
		File bai = idxBam.outputFile
		File vcf = GATHERVCFFILES.outputFile
		File vcfIdx = GATHERVCFFILES.outputFileIdx
	}

################################################################################
}
