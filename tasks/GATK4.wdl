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
			description: 'Input file (SAM or BAM) to extract reads from..',
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

### WARNING : depthOfCoverage is on BETA
## 		--calculate-coverage-over-genes,-gene-list: option seems not working
##			properly : https://github.com/broadinstitute/gatk/issues/6714
##		--count-type: only "COUNT_READS" is supported
task depthOfCoverage {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1b"
		date: "2020-07-27"
	}

	input {
		String path_exe = "gatk"
		Array[String]? javaOptions

		File in
		File intervals
		String outputPath
		String? prefix
		String ext = ".bam"
		String suffix = "DoC"
		File referenceFasta
		File referenceFai
		File referenceDict

		File? geneList
		String countType = "COUNT_READS"
		Boolean disableBamIndexCaching = false
		Boolean disableSequenceDictionaryValidation = false
		String intervalMergingRule = "ALL"
		Int maxBaseQuality = 127
		Int minBaseQuality = 0
		Int maxDepthPerSample = 0
		String outputFormat = "CSV"
		String partitionType = "sample"
		Boolean printBaseCounts = false

		Boolean quiet = false
		File? tmpDir
		Boolean useJDKDeflater = false
		Boolean useJDKInflater = false
		String validationStringency = "SILENT"
		String verbosity = "INFO"
		Boolean showHidden = false

		Array[Int] summaryCoverageThreshold = [5]
	}

	Array[String] summaryCoverageThresholdOpt = prefix("--summary-coverage-threshold ", summaryCoverageThreshold)

	String outputName = if defined(prefix) then "~{prefix}.~{suffix}" else basename(in,ext) + ".~{suffix}"
	String outputFile = "~{outputPath}/~{outputName}"

	command <<<
		echo """
			~{path_exe} DepthOfCoverage \
				--java-options '~{sep=" " javaOptions}' \
				--calculate-coverage-over-genes ~{default='null' geneList} \
				--count-type ~{countType} \
				--disable-bam-index-caching ~{disableBamIndexCaching} \
				--disable-sequence-dictionary-validation ~{disableSequenceDictionaryValidation} \
				--interval-merging-rule ~{intervalMergingRule} \
				--max-base-quality ~{maxBaseQuality} \
				--max-depth-per-sample ~{maxDepthPerSample} \
				--min-base-quality ~{minBaseQuality} \
				--output-format ~{outputFormat} \
				--partition-type ~{partitionType} \
				--print-base-counts ~{printBaseCounts} \
				--QUIET '~{quiet}' \
				--use-jdk-deflater '~{useJDKDeflater}' \
				--use-jdk-inflater '~{useJDKInflater}' \
				--read-validation-stringency '~{validationStringency}' \
				--verbosity ~{verbosity} \
				--showHidden ~{showHidden} \
				--input ~{in} \
				--intervals ~{intervals} \
				--reference ~{referenceFasta} \
				--output ~{outputFile} \
				~{sep=" " summaryCoverageThresholdOpt}
		"""
	>>>

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "gatk"]',
			category: 'optional'
		}
		in: {
			description: 'BAM/SAM/CRAM file containing reads.',
			category: 'Required'
		}
		intervals: {
			description: 'Path to a file containing genomic intervals over which to operate. (format intervals list: chr1:1000-2000)',
			category: 'Required'
		}
		outputPath: {
			description: 'Output path where files were generated.',
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
			description: 'Suffix for the output file (e.g. prefix.suffix.bam) [default: "DoC"]',
			category: 'optional'
		}
		referenceFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'required'
		}
		referenceFai: {
			description: 'Path to the reference file index (format: fai)',
			category: 'required'
		}
		referenceDict: {
			description: 'Path to the reference file dict (format: dict)',
			category: 'required'
		}
		geneList: {
			description: 'Calculate coverage statistics over this list of genes. (refseq format)',
			category: 'optional'
		}
		countType: {
			description: 'How should overlapping reads from the same fragment be handled? (Possible values: {COUNT_READS, COUNT_FRAGMENTS, COUNT_FRAGMENTS_REQUIRE_SAME_BASE}) [default: COUNT_READS]',
			category: 'optional'
		}
		disableBamIndexCaching: {
			description: "If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified. [default: false]",
			category: 'optional'
		}
		disableSequenceDictionaryValidation: {
			description: 'If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk! [default: false]',
			category: 'optional'
		}
		intervalMergingRule: {
			description: 'Interval merging rule for abutting intervals. (possible values: ALL, OVERLAPPING_ONLY) [default: ALL]',
			category: 'optional'
		}
		maxBaseQuality: {
			description: 'Maximum quality of bases to count towards depth. [default: 127]',
			category: 'optional'
		}
		minBaseQuality: {
			description: 'Minimum quality of bases to count towards depth. [default: 0]',
			category: 'optional'
		}
		maxDepthPerSample: {
			description: 'Maximum number of reads to retain per sample per locus. Reads above this threshold will be downsampled. Set to 0 to disable. [default: 0]',
			category: 'optional'
		}
		outputFormat: {
			description: 'The format of the output file. (possible values: CSV, TABLE) [default: CSV]',
			category: 'optional'
		}
		partitionType: {
			description: 'Partition type for depth of coverage. (possbile values: sample, readgroup and/or library) [default: sample]',
			category: 'optional'
		}
		printBaseCounts: {
			description: 'Add base counts to per-locus output. [default: false]',
			category: 'optional'
		}
		quiet: {
			description: 'Whether to suppress job-summary info on System.err. [Default: false]',
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
