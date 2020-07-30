version 1.0

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
	}

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

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "fastqc"]',
			category: 'optional'
		}
		path_java: {
			description: 'Provides the full path to the java binary you want to use to launch fastqc. [default: assuming is in the path]',
			category: 'optional'
		}
		in: {
			description: 'A set of sequence files (one or more)',
			category: 'Required'
		}
		extract: {
			description: 'The zip file will be uncompressed [default: false]',
			category: 'optional'
		}
		nogroup: {
			description: 'Disable grouping of bases for reads >50bp. [default: false]',
			category: 'optional'
		}
		minLength: {
			description: 'Sets an artificial lower limit on the length of the sequence to be shown in the report.',
			category: 'optional'
		}
		format: {
			description: 'Bypasses the normal sequence file format detection and forces the program to use the specified format.',
			category: 'optional'
		}
		kmers: {
			description: 'Specifies the length of Kmer to look for in the Kmer content module (between 2 and 10). [default: 7]',
			category: 'optional'
		}
		contaminants: {
			description: 'Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against.',
			category: 'optional'
		}
		adapters: {
			description: 'Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library.',
			category: 'optional'
		}
		limits: {
			description: 'Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules.',
			category: 'optional'
		}
		tempDir: {
			description: 'Selects a directory to be used for temporary files written when generating report images.',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
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
	}

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

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "fastqc"]',
			category: 'optional'
		}
		path_java: {
			description: 'Provides the full path to the java binary you want to use to launch fastqc. [default: assuming is in the path]',
			category: 'optional'
		}
		in: {
			description: 'A set of sequence files (one or more)',
			category: 'Required'
		}
		extract: {
			description: 'The zip file will be uncompressed [default: false]',
			category: 'optional'
		}
		minLength: {
			description: 'Sets an artificial lower limit on the length of the sequence to be shown in the report.',
			category: 'optional'
		}
		format: {
			description: 'Bypasses the normal sequence file format detection and forces the program to use the specified format.',
			category: 'optional'
		}
		kmers: {
			description: 'Specifies the length of Kmer to look for in the Kmer content module (between 2 and 10). [default: 7]',
			category: 'optional'
		}
		contaminants: {
			description: 'Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against.',
			category: 'optional'
		}
		adapters: {
			description: 'Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library.',
			category: 'optional'
		}
		limits: {
			description: 'Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules.',
			category: 'optional'
		}
		tempDir: {
			description: 'Selects a directory to be used for temporary files written when generating report images.',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
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
	}

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

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "fastqc"]',
			category: 'optional'
		}
		path_java: {
			description: 'Provides the full path to the java binary you want to use to launch fastqc. [default: assuming is in the path]',
			category: 'optional'
		}
		in: {
			description: 'A set of sequence files (one or more)',
			category: 'Required'
		}
		extract: {
			description: 'The zip file will be uncompressed [default: false]',
			category: 'optional'
		}
		nogroup: {
			description: 'Disable grouping of bases for reads >50bp. [default: false]',
			category: 'optional'
		}
		nofilter: {
			description: "Don't remove read flagged by casava as poor quality when performing the QC analysis. [default: false]",
			category: 'optional'
		}
		minLength: {
			description: 'Sets an artificial lower limit on the length of the sequence to be shown in the report.',
			category: 'optional'
		}
		format: {
			description: 'Bypasses the normal sequence file format detection and forces the program to use the specified format.',
			category: 'optional'
		}
		kmers: {
			description: 'Specifies the length of Kmer to look for in the Kmer content module (between 2 and 10). [default: 7]',
			category: 'optional'
		}
		contaminants: {
			description: 'Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against.',
			category: 'optional'
		}
		adapters: {
			description: 'Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library.',
			category: 'optional'
		}
		limits: {
			description: 'Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules.',
			category: 'optional'
		}
		tempDir: {
			description: 'Selects a directory to be used for temporary files written when generating report images.',
			category: 'optional'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'optional'
		}
	}
}