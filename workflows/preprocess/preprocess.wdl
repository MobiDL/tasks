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

import "getAnnotationsFiles.wdl" as getAnnotationsFiles
import "getGenomeDir.wdl" as getGenomeDir
import "../../tasks/star.wdl" as star

workflow preprocess {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		date: "2020-12-04"
		version: "0.0.1"
	}

	input {
		String outputPathFasta
		String outputPathAnnot
		String outputPathStar = outputPathFasta

		String linkFa
		String linkGTF
		Array[String] vcfs_bundle_broad

		Int readLengthRNA = 100
	}

################################################################################

	call getGenomeDir.getGenomeDir as GGD {
		input :
			outputPath = outputPathFasta,
			linkFa = linkFa
	}

	call getAnnotationsFiles.getAnnotationsFiles as GAF {
		input :
			outputPath = outputPathAnnot,
			linkGTF = linkGTF,
			vcfs_bundle_broad = vcfs_bundle_broad
	}

	call star.genomeGenerate as SGG {
		input :
			outputPath = outputPathStar,
			refFasta = GGD.refFasta,
			refGTF = GAF.gtf,
			readLength = readLengthRNA
	}

################################################################################

	output {
		File refFasta = GGD.refFasta
		File refFai = GGD.refFai
		File refDict = GGD.refDict
		File refAmb = GGD.refAmb
		File refAnn = GGD.refAnn
		File refBwt = GGD.refBwt
		File refPac = GGD.refPac
		File refSa = GGD.refSa

		File gtf = GAF.gtf
		Array[File] vcfs = GAF.vcfs
		Array[File] vcfsgz = GAF.vcfsgz
		Array[File] vcfsgztbi = GAF.vcfsgztbi
		File vep_infos = GAF.vep_infos

		File chrLength = SGG.chrLength
		File chrNameLength = SGG.chrNameLength
		File chrName = SGG.chrName
		File chrStart = SGG.chrStart
		File exonGeTrInfo = SGG.exonGeTrInfo
		File exonInfo = SGG.exonInfo
		File geneInfo = SGG.geneInfo
		File Genome = SGG.Genome
		File genomeParameters = SGG.genomeParameters
		File SA = SGG.SA
		File SAindex = SGG.SAindex
		File sjdbInfo = SGG.sjdbInfo
		File sjdbListGTF = SGG.sjdbListGTF
		File sjdbList = SGG.sjdbList
		File transcriptInfo = SGG.transcriptInfo

	}

	parameter_meta {
		outputPathFasta : {
			outputPath: 'Path where fatsa files will be written.',
			category: 'Output'
		}
		outputPathAnnot : {
			description: 'Path where annotations files will be written.',
			category: 'Required'
		}
		outputPathStar : {
			description: 'Path where STAR files will be written. [default: outputPathFasta]',
			category: 'Option'
		}
		linkFa : {
			description: 'Link to a fasta to download (format: *.fa.gz)',
			category: 'Input'
		}
		linkGTF : {
			description: 'Link to a gtf to download (format: *.gtf.gz)',
			category: 'Required'
		}
		vcfs_bundle_broad : {
			description: 'Link to vcfs in the broad bundle (ftp://ftp.broadinstitute.org/bundle/hg19/, format: *.vcf.gz)',
			category: 'Required'
		}
		readLength: {
			description: 'Read length of the sequencing [default: 100]',
			category: 'Option'
		}
	}
}
