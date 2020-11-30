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

task get_version {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-11-20"
	}

	input {
		String path_exe = "star"

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	command <<<
		~{path_exe} --version
	>>>

	output {
		String version = read_string(stdout())
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "star"]',
			category: 'System'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}

task genomeGenerate {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-10-21"
	}

	input {
		String path_exe = "STAR"

		File refFasta
		File refGTF
		String? outputPath
		String? name
		String subString = "(.fa(sta)?)?(.gtf)?"
		String subStringReplace = ""


		Int readLength = 100

		Int threads = 1
		Int memoryByThreads = 768
		String memory = "32G"
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String basenameFa = sub(basename(refFasta),subString,subStringReplace)
	String basenameGTF = sub(basename(refGTF),subString,subStringReplace)
	String baseName = if defined(name) then name else "~{basenameFa}_~{basenameGTF}_~{readLength}"
	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d ~{outputRep} ]]; then
			mkdir -p ~{outputRep}
		fi

		~{path_exe} --runMode genomeGenerate \
			--genomeDir ~{outputRep} \
			--genomeFastaFiles ~{refFasta} \
			--sjdbGTFfile ~{refGTF} \
			--sjdbOverhang ~{readLength - 1} \
			--runThreadN ~{threads}

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

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "samtools"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files will be generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'Tool option'
		}
		subString: {
			description: 'Substring to remove to create name file (used for fasta file name and gtf) [default: "(.fa(sta)?)?(.gtf)?"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'Required'
		}
		refGTF: {
			description: 'Path to the GTF reference file (format: GTF)',
			category: 'Required'
		}
		readLength: {
			description: 'Read length of the sequencing [default: 100]',
			category: 'Tool option'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}
