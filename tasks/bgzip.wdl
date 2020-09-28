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

task bgzip {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-24"
	}

	input {
		String path_exe = "bgzip"

		File in
		String outputPath

		Boolean decompress = false
		Boolean force = false
		Boolean index = false
		Boolean keepFile = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String outputFile = if ! decompress then "~{outputPath}/" + basename(in) + ".gz" else "~{outputPath}/" + basename(in,".gz")
	String outputIndex = if (index && ! decompress) then "~{outputFile}.gzi" else "~{outputFile}"
	String indexOpt = if (index && ! decompress) then "--index --index-name ~{outputIndex}" else ""

	command <<<

		~{path_exe} \
			~{true="--stdout" false="" keepFile} \
			~{true="--decompress" false="" decompress} \
			~{true="--force" false="" force}\
			~{indexOpt} \
			~{in} \
			--threads ~{threads} > ~{outputFile}

	>>>

	output {
		File out = "~{outputFile}"
		File index = "~{outputIndex}"
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${totalMemMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "bgzip"]',
			category: 'System'
		}
		in: {
			description: "File to compres/decompress.",
			category: "required"
		}
		outputPath: {
			description: "Path where was generated output.",
			category: "required"
		}
		decompress: {
			description: "Decompress file (incompatible with index) [default: false]",
			category: 'optional'
		}
		force: {
			description: "Overwrite files without asking [default: false]",
			category: 'optional'
		}
		index: {
			description: "Compress and create BGZF index [default: false]",
			category: 'optional'
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
