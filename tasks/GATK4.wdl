version 1.0

task reorderSam {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-07-27"
	}

	input {
		String path_exe = "gatk"
		Array[String]? javaOptions

		File in
		String outputPath
		String? prefix
		String ext = ".bam"
		String suffix = "reorder"
		File sequenceDict

		Boolean allowIncompleteLengthDiscordance = false
		Boolean allowIncompleteDictDiscordance = false
		Int compressionLevel = 2
		Boolean createIndex = false
		Boolean createMD5 = false
		File? GA4GHClientSecrets
		Int maxRecordsInRam = 500000
		Boolean quiet = false
		File? referenceSequence
		Array[File]? tmpDir
		Boolean useJDKDeflater = false
		Boolean useJDKInflater = false
		String validationStringency = "STRICT"
		String verbosity = "INFO"
		Boolean showHidden = false
	}

	String GA4GHClientSecretsOpt = if defined(GA4GHClientSecrets) then "--GA4GH_CLIENT_SECRETS ~{GA4GHClientSecrets} " else ""
	String referenceSequenceOpt = if defined(referenceSequence) then "--REFERENCE_SEQUENCE ~{referenceSequence} " else ""
	String tmpDirOpt = if defined(tmpDir) then "~{prefix='--TMP_DIR ' sep=' ' tmpDir} " else ""
	String prefixTmpDirOpt = if defined(tmpDir) then "--TMP_DIR " else ""

	String outputName = if defined(prefix) then "~{prefix}.~{suffix}~{ext}" else basename(in,ext) + ".~{suffix}~{ext}"
	String outputFile = "~{outputPath}/~{outputName}"

	command <<<

		if [[ ! -f ~{outputFile} ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} ReorderSam ~{GA4GHClientSecrets}~{referenceSequenceOpt}~{prefixTmpDirOpt}~{sep=' --TMP_DIR ' tmpDir} \
			--java-options '~{sep=" " javaOptions}' \
			--ALLOW_CONTIG_LENGTH_DISCORDANCE '~{allowIncompleteLengthDiscordance}' \
			--ALLOW_INCOMPLETE_DICT_CONCORDANCE '~{allowIncompleteDictDiscordance}' \
			--COMPRESSION_LEVEL ~{compressionLevel} \
			--CREATE_INDEX '~{createIndex}' \
			--CREATE_MD5_FILE '~{createMD5}' \
			--MAX_RECORDS_IN_RAM ~{maxRecordsInRam} \
			--QUIET '~{quiet}' \
			--USE_JDK_DEFLATER '~{useJDKDeflater}' \
			--USE_JDK_INFLATER '~{useJDKInflater}' \
			--VALIDATION_STRINGENCY '~{validationStringency}' \
			--VERBOSITY ~{verbosity} \
			--showHidden ~{showHidden} \
			--INPUT ~{in} \
			--SEQUENCE_DICTIONARY ~{sequenceDict} \
			--OUTPUT ~{outputFile}

	>>>

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "gatk"]',
			category: 'optional'
		}
		in: {
			description: 'Output path where file (SAM or BAM) were generated.',
			category: 'Required'
		}
		outputPath: {
			description: 'Output path where file (SAM or BAM) were generated.',
			category: 'Required'
		}
		ext: {
			description: 'Extension of the input file (".sam" or ".bam") [default: ".bam"]',
			category: 'optional'
		}
		prefix: {
			description: 'Prefix for the output file [default: basename(in, ext)]',
			category: 'optional'
		}
		suffix: {
			description: 'Suffix for the output file (e.g. sample.suffix.bam) [default: "reorder"]',
			category: 'optional'
		}
		sequenceDict: {
			description: 'Sequence Dictionary for the OUTPUT file (can be read from one of the following file types (SAM, BAM, VCF, BCF, Interval List, Fasta, or Dict)',
			category: 'Required'
		}
		allowIncompleteLengthDiscordance: {
			description: 'If true, then permits mapping from a read contig to a new reference contig with the same name but a different length.  Highly dangerous, only use if you know what you are doing. [Default: false]',
			category: 'optional'
		}
		allowIncompleteDictDiscordance: {
			description: 'If true, allows only a partial overlap of the original contigs with the new reference sequence contigs.  By default, this tool requires a corresponding contig in the new reference for each read contig. [Default: false]',
			category: 'optional'
		}
		compressionLevel: {
			description: 'Compression level for all compressed files created (e.g. BAM and VCF). [default: 2]',
			category: 'optional'
		}
		createIndex: {
			description: 'Whether to create a BAM index when writing a coordinate-sorted BAM file. [Default: false]',
			category: 'optional'
		}
		createMD5: {
			description: 'Whether to create an MD5 digest for any BAM or FASTQ files created. [Default: false]',
			category: 'optional'
		}
		GA4GHClientSecrets: {
			description: 'Google Genomics API client_secrets.json file path. [Default: null]',
			category: 'optional'
		}
		maxRecordsInRam: {
			description: 'When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed. [Default: 500000]',
			category: 'optional'
		}
		referenceSequence: {
			description: 'Reference sequence file. [Default: null]',
			category: 'optional'
		}
		tmpDir: {
			description: 'One or more directories with space available to be used by this program for temporary storage of working files. [Default: null]',
			category: 'optional'
		}
		useJDKDeflater: {
			description: 'Use the JDK Deflater instead of the Intel Deflater for writing compressed output. [Default: false]',
			category: 'optional'
		}
		useJDKInflater: {
			description: 'Use the JDK Inflater instead of the Intel Inflater for reading compressed input. [Default: false]',
			category: 'optional'
		}
		validationStringency: {
			description: ' Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. [Default: "STRIC"]',
			category: 'optional'
		}
		verbosity: {
			description: 'Control verbosity of logging. [Default: "INFO"]',
			category: 'optional'
		}
		showHidden: {
			description: 'Display hidden arguments. [Default: false]',
			category: 'optional'
		}
	}
}
