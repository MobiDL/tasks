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
		String path_exe = "fastqc"

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
		String output = read_string(stdout())
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bcftools"]',
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

task fastqc {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-29"
	}

	input {
		String path_exe = "fastqc"
		String? path_java

		Array[File]+ in
		String outputPath = "."

		Boolean extract = false
		Boolean nogroup = false
		Int? minLength
		String? format
		Int kmers = 7

		File? contaminants
		File? adapters
		File? limits

		String? tempDir

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

		if [[ ! -d ~{outputPath} ]]; then
			mkdir -p ~{outputPath}
		fi

		~{path_exe} ~{default="" "--java " + path_java} \
			~{true="--extract" false="--noextract" extract} \
			~{true="--nogroup" false="" nogroup} \
			~{default="" "--min_length " + minLength} \
			~{default="" "--format " + format} \
			--kmers ~{kmers} \
			~{default="" "--contaminants " + contaminants} \
			~{default="" "--adapters " + adapters} \
			~{default="" "--limits " + limits} \
			~{default="" "--dir " + tempDir} \
			--threads ~{threads} \
			--outdir ~{outputPath} \
			~{sep=" " in}

	>>>

	output {
		Array[File] outHTML = glob(outputPath + "/*.html")
		Array[File] outZIP = glob(outputPath + "/*.zip")
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "fastqc"]',
			category: 'System'
		}
		path_java: {
			description: 'Provides the full path to the java binary you want to use to launch fastqc. [default: assuming is in the path]',
			category: 'System'
		}
		in: {
			description: 'A set of sequence files (one or more)',
			category: 'Required'
		}
		extract: {
			description: 'The zip file will be uncompressed [default: false]',
			category: 'Tool option'
		}
		nogroup: {
			description: 'Disable grouping of bases for reads >50bp. [default: false]',
			category: 'Tool option'
		}
		minLength: {
			description: 'Sets an artificial lower limit on the length of the sequence to be shown in the report.',
			category: 'Tool option'
		}
		format: {
			description: 'Bypasses the normal sequence file format detection and forces the program to use the specified format.',
			category: 'Tool option'
		}
		kmers: {
			description: 'Specifies the length of Kmer to look for in the Kmer content module (between 2 and 10). [default: 7]',
			category: 'Tool option'
		}
		contaminants: {
			description: 'Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against.',
			category: 'Tool option'
		}
		adapters: {
			description: 'Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library.',
			category: 'Tool option'
		}
		limits: {
			description: 'Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules.',
			category: 'Tool option'
		}
		tempDir: {
			description: 'Selects a directory to be used for temporary files written when generating report images.',
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

task fastqcNano {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-29"
	}

	input {
		String path_exe = "fastqc"
		String? path_java

		Array[File]+ in
		String outputPath = "."

		Boolean extract = false
		Int? minLength
		String? format
		Int kmers = 7

		File? contaminants
		File? adapters
		File? limits

		String? tempDir

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

		if [[ ! -d ~{outputPath} ]]; then
			mkdir -p ~{outputPath}
		fi

		~{path_exe} ~{default="" "--java " + path_java} --nano \
			~{true="--extract" false="--noextract" extract} \
			~{default="" "--min_length " + minLength} \
			~{default="" "--format " + format} \
			--kmers ~{kmers} \
			~{default="" "--contaminants " + contaminants} \
			~{default="" "--adapters " + adapters} \
			~{default="" "--limits " + limits} \
			~{default="" "--dir " + tempDir} \
			--threads ~{threads} \
			--outdir ~{outputPath} \
			~{sep=" " in}

	>>>

	output {
		Array[File] outHTML = glob(outputPath + "/*.html")
		Array[File] outZIP = glob(outputPath + "/*.zip")
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "fastqc"]',
			category: 'System'
		}
		path_java: {
			description: 'Provides the full path to the java binary you want to use to launch fastqc. [default: assuming is in the path]',
			category: 'System'
		}
		in: {
			description: 'A set of sequence files (one or more)',
			category: 'Required'
		}
		extract: {
			description: 'The zip file will be uncompressed [default: false]',
			category: 'Tool option'
		}
		minLength: {
			description: 'Sets an artificial lower limit on the length of the sequence to be shown in the report.',
			category: 'Tool option'
		}
		format: {
			description: 'Bypasses the normal sequence file format detection and forces the program to use the specified format.',
			category: 'Tool option'
		}
		kmers: {
			description: 'Specifies the length of Kmer to look for in the Kmer content module (between 2 and 10). [default: 7]',
			category: 'Tool option'
		}
		contaminants: {
			description: 'Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against.',
			category: 'Tool option'
		}
		adapters: {
			description: 'Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library.',
			category: 'Tool option'
		}
		limits: {
			description: 'Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules.',
			category: 'Tool option'
		}
		tempDir: {
			description: 'Selects a directory to be used for temporary files written when generating report images.',
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

task fastqcCasava {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-29"
	}

	input {
		String path_exe = "fastqc"
		String? path_java

		Array[File]+ in
		String outputPath = "."

		Boolean extract = false
		Boolean nogroup = false
		Boolean nofilter = false
		Int? minLength
		String? format
		Int kmers = 7

		File? contaminants
		File? adapters
		File? limits

		String? tempDir

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

		if [[ ! -d ~{outputPath} ]]; then
			mkdir -p ~{outputPath}
		fi

		~{path_exe} ~{default="" "--java " + path_java} --casava \
			~{true="--extract" false="--noextract" extract} \
			~{true="--nogroup" false="" nogroup} \
			~{true="--nofilter" false="" nofilter} \
			~{default="" "--min_length " + minLength} \
			~{default="" "--format " + format} \
			--kmers ~{kmers} \
			~{default="" "--contaminants " + contaminants} \
			~{default="" "--adapters " + adapters} \
			~{default="" "--limits " + limits} \
			~{default="" "--dir " + tempDir} \
			--threads ~{threads} \
			--outdir ~{outputPath} \
			~{sep=" " in}

	>>>

	output {
		Array[File] outHTML = glob(outputPath + "/*.html")
		Array[File] outZIP = glob(outputPath + "/*.zip")
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "fastqc"]',
			category: 'System'
		}
		path_java: {
			description: 'Provides the full path to the java binary you want to use to launch fastqc. [default: assuming is in the path]',
			category: 'System'
		}
		in: {
			description: 'A set of sequence files (one or more)',
			category: 'Required'
		}
		extract: {
			description: 'The zip file will be uncompressed [default: false]',
			category: 'Tool option'
		}
		nogroup: {
			description: 'Disable grouping of bases for reads >50bp. [default: false]',
			category: 'Tool option'
		}
		nofilter: {
			description: "Don't remove read flagged by casava as poor quality when performing the QC analysis. [default: false]",
			category: 'Tool option'
		}
		minLength: {
			description: 'Sets an artificial lower limit on the length of the sequence to be shown in the report.',
			category: 'Tool option'
		}
		format: {
			description: 'Bypasses the normal sequence file format detection and forces the program to use the specified format.',
			category: 'Tool option'
		}
		kmers: {
			description: 'Specifies the length of Kmer to look for in the Kmer content module (between 2 and 10). [default: 7]',
			category: 'Tool option'
		}
		contaminants: {
			description: 'Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against.',
			category: 'Tool option'
		}
		adapters: {
			description: 'Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library.',
			category: 'Tool option'
		}
		limits: {
			description: 'Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules.',
			category: 'Tool option'
		}
		tempDir: {
			description: 'Selects a directory to be used for temporary files written when generating report images.',
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
