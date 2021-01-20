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

import "../../tasks/rsync.wdl" as rsync
import "../../tasks/bash.wdl" as bash
import "../../tasks/bgzip.wdl" as bgzip
import "../../tasks/tabix.wdl" as tabix
import "../../tasks/vep.wdl" as vep

workflow getAnnotationsFiles {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		date: "2020-12-04"
		last_date: "2021-01-20"
		version: "0.0.2"
	}

	input {
		String outputPath
		String linkGTF
		Array[String] vcfs_bundle_broad
	}

################################################################################

	call bash.wget as getGTF {
		input :
			in = linkGTF,
			outputPath = outputPath
	}

	call bash.gzip as gunzipGTF {
		input :
			in = getGTF.outputFile,
			outputPath = outputPath,
			decompress = true
	}

	call bash.sortgtf as sortGTF {
		input :
			in = gunzipGTF.outputFile,
			outputPath = outputPath
	}

	call bgzip.compress as bgzipGTF {
		input :
			in = sortGTF.outputFile,
			outputPath = outputPath
	}

	call tabix.index as indexGTF {
		input :
			in = bgzipGTF.outputFile,
			outputPath = outputPath
	}

################################################################################

################################################################################

	scatter (vcfs in select_all(vcfs_bundle_broad)) {
		call bash.wget as getVCF_bundle {
			input :
				in = vcfs,
				user = "gsapubftp-anonymous",
				outputPath = outputPath
		}

		call bash.gzip as gunzipVCF_bundle {
			input :
				in = getVCF_bundle.outputFile,
				decompress = true,
				outputPath = outputPath
		}

		call bgzip.compress as bgzipVCF_bundle {
			input :
				in = gunzipVCF_bundle.outputFile,
				outputPath = outputPath
		}

		call tabix.index as indexVCF_bundle {
			input :
				in = bgzipVCF_bundle.outputFile,
				outputPath = outputPath
		}
	}
################################################################################

################################################################################

	call vep.install as vep {
		input :
			destDir = outputPath
	}

################################################################################

################################################################################

	call snpEff.install as snpEff {
		input :
			destDir = outputPath,
			genome = "hg19"
	}

################################################################################

	output {
		File gtf = gunzipGTF.outputFile
		Array[File] vcfs = gunzipVCF_bundle.outputFile
		Array[File] vcfsgz = bgzipVCF_bundle.outputFile
		Array[File] vcfsgztbi = indexVCF_bundle.outputFile
		File vep_infos = vep.info
	}

	parameter_meta {
		outputPath : {
			description: 'Path where files will be written.',
			category: 'Required'
		}
		linkGTF : {
			description: 'Link to a gtf to download (format: *.gtf.gz)',
			category: 'Required'
		}
		vcfs_bundle_broad : {
			description: 'Link to vcfs in the broad bundle (ftp://ftp.broadinstitute.org/bundle/hg19/, format: *.vcf.gz)',
			category: 'Required'
		}
	}
}
