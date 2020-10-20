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

task findFiles {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-24"
	}

	input {
		String path

		String? regexpName
		String? regexpPath

		Int? maxDepth
		Int? minDepth

		Int? uid
		Int? gid

		Boolean? readable
		Boolean? writable
		Boolean? executable

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String regexpNameOpt = if defined(regexpName) then "-name \"~{regexpName}\" " else ""
	String regexpPathOpt = if defined(regexpPath) then "-path \"~{regexpPath}\" " else ""

	String maxDepthOpt = if defined(maxDepth) then "-maxdepth ~{maxDepth} " else ""
	String minDepthOpt = if defined(minDepth) then "-mindepth ~{minDepth} " else ""

	String uidOpt = if defined(uid) then "-uid ~{uid} " else ""
	String gidOpt = if defined(gid) then "-gid ~{gid} " else ""

	String readableOpt = if (defined(readable) && readable) then "-readable " else ""
	String writableOpt = if (defined(writable) && writable) then "-writable " else ""
	String executableOpt = if (defined(executable) && executable) then "-executable " else ""

	command <<<

		find ~{path} ~{maxDepthOpt}~{minDepthOpt}~{uidOpt}~{gidOpt}~{readableOpt}~{writableOpt}~{executableOpt}~{regexpPathOpt}~{regexpNameOpt}-fprint files.txt

	>>>

	output {
		Array[File] files = read_lines("files.txt")
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path: {
			description: "Path where find will work on.",
			category: 'Required'
		}
		regexpName: {
			description: "Base of file name (the path with the leading directories removed) matches shell pattern pattern.",
			category: 'Tool option'
		}
		regexpPath: {
			description: "File name matches shell pattern pattern. The metacharacters do not treat `/' or `.' specially.",
			category: 'Tool option'
		}
		maxDepth: {
			description: "Descend at most levels (a non-negative integer) levels of directories below the starting-points.",
			category: 'Tool option'
		}
		minDepth: {
			description: "Do not apply any tests or actions at levels less than levels (a non-negative integer).",
			category: 'Tool option'
		}
		uid: {
			description: "File's numeric user ID is n.",
			category: 'Tool option'
		}
		gid: {
			description: "File's numeric group ID is n.",
			category: 'Tool option'
		}
		readable: {
			description: "Matches files which are readable.",
			category: 'Tool option'
		}
		writable: {
			description: "Matches files which are writable.",
			category: 'Tool option'
		}
		executable: {
			description: "Matches files which are executable and directories which are searchable (in a file name resolution sense).",
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

task convertBedToIntervals {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-30"
	}

	input {
		File in

		String? outputPath
		String? name

		String ext =".intervals"

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String outputName = if defined(name) then name else sub(basename(in),".bed", "")
	String outputFile = if defined(outputPath) then outputPath + "/" + outputName + ext else outputName + ext

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		grep -E -v "^track|^browser|^#" ~{in} | \
		awk -F "\t" '{
			if($2 == $3){
				$3++
			}
			$2++;
			print $1":"$2"-"$3;
		}' > ~{outputFile}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		in: {
			description: 'Input bed to convert.',
			category: 'Required'
		}
		outputPath: {
			description: 'Path where was generated output. [default: "."]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Basename of the output file. [default: sub(basename(in), ".bed")]',
			category: 'Output path/name option'
		}
		ext: {
			description: 'Extension of the output file. [default: ".intervals"]',
			category: 'Output path/name option'
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

task makeLink {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-10"
	}

	input {
		File in

		String outputPath
		String? name

		Boolean softLink = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String outputName = if defined(name) then name else basename(in)
	String outputFile = outputPath + "/" + outputName

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		if [[ ! -f ~{outputFile} ]]; then
			ln \
				~{true="-s" false="" softLink} \
				~{in} \
				~{outputFile}
		else
			1>&2 echo "[WARN][$(date +'%Y-%m-%d')] Link of ~{outputFile} failed : file already exist !"
		fi

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		in: {
			description: 'Input bed to convert.',
			category: 'Required'
		}
		outputPath: {
			description: 'Path where was generated output.',
			category: 'Output path/name option'
		}
		name: {
			description: 'Basename of the output file. [default: basename(in)]',
			category: 'Output path/name option'
		}
		softLink: {
			description: 'Make soft link (-s). [default: false]',
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

task concatenateFiles {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-08-19"
	}

	input {
		Array[File] in

		String? outputPath
		String? name
		String subString = "^[0-9]+\-"

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(in[0]),subString,"")
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		cat ~{sep=" " in} > ~{outputFile}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		in: {
			description: 'Array of files.',
			category: 'Required'
		}
		outputPath: {
			description: 'Path where output will be generated.',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name of the output file. [default: sub(basename(in[0]),subString,"")]',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Substring to remove to create name file [default: "^[0-9]+\-"]',
			category: 'Output path/name option'
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

task wget {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-10-05"
	}

	input {
		String path_exe = "wget"

		String in
		String? outputPath
		String? name
		String subString = ".*/(.*)$"
		String subStringReplace = "$1"

		Boolean forceHTML = false
		File? configuration
		Boolean timestamping = false
		String? user
		String? password
		Boolean recursive = false
		Int? recursiveLVL
		Boolean convertLinks = false
		Boolean mirror = false
		Boolean pageRequisites = false

		Array[String]? extAccept
		Array[String]? extReject
		Array[String]? repAccept
		Array[String]? repReject

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(in,subString,subStringReplace)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	Boolean extAcceptDefined = defined(extAccept)
	Boolean extRejectDefined = defined(extReject)
	Boolean repAcceptDefined = defined(repAccept)
	Boolean repRejectDefined = defined(repReject)

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} \
			~{true="--force-html" false="" forceHTML} \
			~{default="" "--config " + configuration} \
			~{true="--timestamping" false="" timestamping} \
			~{default="" "--user " + user} \
			~{default="" "--password " + password} \
			~{true="--recursive" false="" recursive} \
			~{default="" "--level " + recursiveLVL} \
			~{true="--convert-links" false="" convertLinks} \
			~{true="--mirror" false="" mirror} \
			~{true="--page-requisites" false="" pageRequisites} \
			~{true="--accept" false="" extAcceptDefined} ~{default="" sep="," extAccept} \
			~{true="--reject" false="" extRejectDefined} ~{default="" sep="," extReject} \
			~{true="--include-directories" false="" repAcceptDefined} ~{default="" sep="," repAccept} \
			~{true="--exclude-directories" false="" repRejectDefined} ~{default="" sep="," repReject} \
			~{in} \
			-O ~{outputFile}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "wget"]',
			category: 'System'
		}
		in: {
			description: 'Url of the input file.',
			category: 'Required'
		}
		outputPath: {
			description: 'Path where output will be generated.',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name of the output file. [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Substring to remove to create name file [default: ".*/(.*)$"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: "$1"]',
			category: 'Output path/name option'
		}
		forceHTML: {
			description: 'Force to resolve inputs as HTML [default: false]',
			category: 'Tool option'
		}
		configuration: {
			description: 'Use a specific config file.',
			category: 'Tool option'
		}
		timestamping: {
			description: 'Do not download files if they are not more recent than locally [default: false]',
			category: 'Tool option'
		}
		user: {
			description: 'Define the user for FTP and HTTP (use with caution because it will appear in clear)',
			category: 'Tool option'
		}
		password: {
			description: 'Define the password for FTP and HTTP (use with caution because it will appear in clear)',
			category: 'Tool option'
		}
		recursive: {
			description: 'Turn on/off recursive mode [default: false]',
			category: 'Tool option'
		}
		recursiveLVL: {
			description: 'Level of maximal recursivity',
			category: 'Tool option'
		}
		convertLinks: {
			description: 'Convert links locally into HTML and CSS [default: false]',
			category: 'Tool option'
		}
		mirror: {
			description: 'Activate timestamping and recursive to infinity [default: false]',
			category: 'Tool option'
		}
		pageRequisites: {
			description: 'Get all images and others to display HTML pages [default: false]',
			category: 'Tool option'
		}
		extAccept: {
			description: 'Array of accepted extensions.',
			category: 'Tool option'
		}
		extReject: {
			description: 'Array of rejected extensions.',
			category: 'Tool option'
		}
		repAccept: {
			description: 'Array of accepted directories.',
			category: 'Tool option'
		}
		repReject: {
			description: 'Array of rejected directories.',
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

task gzip {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-10-06"
	}

	input {
		String path_exe = "gzip"

		String in
		String? outputPath
		String? name
		String subString = ""
		String subStringReplace = ""
		String ext = ".gz"

		Boolean decompress = false
		Int levelCompression = 6
		Boolean recursive = false


		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseNameTemp = if defined(name) then name else sub(basename(in),subString,subStringReplace)
	String baseName = if decompress then sub(baseNameTemp,ext,"") else "~{baseNameTemp}.~{ext}"
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} \
			~{true="--decompress" false="" decompress} \
			-~{levelCompression} \
			--stdout \
			~{in} \
			> ~{outputFile}

	>>>

	output {
		File outputFile = outputFile
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "gzip"]',
			category: 'System'
		}
		in: {
			description: 'Path of the input file.',
			category: 'Required'
		}
		outputPath: {
			description: 'Path where output will be generated.',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name of the output file. [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Substring to remove to create name file [default: ""]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		ext: {
			description: 'Extension of the output (compress) or input (decompress) file [default: ".gz"]',
			category: 'Tool option'
		}
		decompress: {
			description: 'Enable decompression mode [default: false]',
			category: 'System'
		}
		levelCompression: {
			description: 'Set compression level [default: 6]',
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
