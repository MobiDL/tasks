version 1.0

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
	}

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

    parameter_meta {

        path: {
			description: "Path used as executable [default: 'bcftools']",
			category: "required"
		}
        regexpName: {
			description: "Base of file name (the path with the leading directories removed) matches shell pattern pattern.",
			category: "tests"
		}
        regexpPath: {
			description: "File name matches shell pattern pattern.  The metacharacters do not treat `/' or `.' specially.",
			category: "tests"
		}
        maxDepth: {
			description: "Descend at most levels (a non-negative integer) levels of directories below the starting-points.",
			category: "global"
		}
        minDepth: {
			description: "Do not apply any tests or actions at levels less than levels (a non-negative integer).",
			category: "global"
		}
        uid: {
			description: "File's numeric user ID is n.",
			category: "tests"
		}
        gid: {
			description: "File's numeric group ID is n.",
			category: "tests"
		}
        readable: {
			description: "Matches  files  which  are  readable.",
			category: "global"
		}
        writable: {
			description: "Matches  files  which  are  writable.",
			category: "tests"
		}
        executable: {
			description: "Matches files which are executable and directories which are searchable (in a file name resolution sense).",
			category: "tests"
		}

	}
}
