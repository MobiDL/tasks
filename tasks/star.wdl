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

task genomeGenerate {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-13"
	}

	input {
		String path_exe = "STAR"

		File refFasta
		String? outputPath

		File refGTF

		Int readLength = 100

		Int threads = 1
		String memory = "32G"
	}

	command <<<

		if [[ ! -d ~{outputPath} ]]; then
			mkdir -p ~{outputPath}
		fi

		~{path_exe} --runMode genomeGenerate \
			--genomeDir ~{outputPath} \
			--genomeFastaFiles ~{refFasta} \
			--sjdbGTFfile ~{refGTF} \
			--sjdbOverhang ~{readLength - 1}

	>>>

	output {
		 File chrLength = outputPath + "chrLength.txt"
		 File chrNameLength = outputPath + "chrNameLength.txt"
		 File chrName = outputPath + "chrName.txt"
		 File chrStart = outputPath + "chrStart.txt"
		 File exonGeTrInfo = outputPath + "exonGeTrInfo.tab"
		 File exonInfo = outputPath + "exonInfo.tab"
		 File geneInfo = outputPath + "geneInfo.tab"
		 File Genome = outputPath + "Genome"
		 File genomeParameters = outputPath + "genomeParameters.txt"
		 File SA = outputPath + "SA"
		 File SAindex = outputPath + "SAindex"
		 File sjdbInfo = outputPath + "sjdbInfo.txt"
		 File sjdbListGTF = outputPath + "sjdbList.fromGTF.out.tab"
		 File sjdbList = outputPath + "sjdbList.out.tab"
		 File transcriptInfo = outputPath + "transcriptInfo.tab"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "samtools"]',
			category: 'optional'
		}
		outputPath: {
			description: 'Output path where files will be generated. [default: pwd()]',
			category: 'optional'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'required'
		}
		refGTF: {
			description: 'Path to the GTF reference file (format: GTF)',
			category: 'required'
		}
		readLength: {
			description: 'Read length of the sequencing [default: 100]',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: 32G]',
			category: 'optional'
		}
	}
}
