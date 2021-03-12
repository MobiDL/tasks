version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2021 MoBiDiC
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
			description: 'Path used as executable [default: "star"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files will be generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'Output path/name option'
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

task alignReads {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-12"
	}

	input {
		String path_exe = "STAR"

		File fastqR1
		File? fastqR2

		String? outputPath
		String? name
		String subString = "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"
		String subStringReplace = ""

		# Genome Parameters
		File genomeDir
		Boolean genomeLoad = false
		File? genomeFastaFiles

		# Splice Junctions Database
		File? sjdbGTF
		String? sjdbGTFchrPrefix
		String sjdbGTFfeatureExon = "exon"
		String sjdbGTFtagExonParentTranscript = "transcript_id"
		String sjdbGTFtagExonParentGene = "gene_id"
		Array[String] sjdbGTFtagExonParentGeneName = ["gene_name"]
		Array[String] sjdbGTFtagExonParentGeneType = ["gene_type","gene_biotype"]
		Int sjdbOverhang = 100
		Int sjdbScore = 2
		Boolean sjdbInsertSaveAll = false
		File? sjdbBED

		# Read Parameters
		String readFilesCommand = "zcat"
		Int readQualityScoreBase = 33

		# Read Clipping
		Boolean clipAdapterHamming = true
		Array[Int] clip3pNbases = [0]
		Array[String] clip3pAdapterSeq = ["-"]
		Array[Float] clip3pAdapterMMp = [0.1]
		Array[Int] clip3pAfterAdapterNbases = [0]
		Array[Int] clip5pNbases = [0]

		# Limits
		Int limitIObufferSize = 150000000
		Int limitOutSAMoneReadBytes = 100000
		Int limitOutSJoneRead = 1000
		Int limitOutSJcollapsed = 1000000
		Int limitBAMsortRAM = 0
		Int limitSjdbInsertNsj = 1000000
		Int limitNreadsSoft = -1

		# Output: general
		Boolean outReadsUnmapped = false

		# Output: SAM and BAM
		Boolean outBam = true
		Boolean sorted = true
		Boolean outSAMmodeQuality = true
		Boolean outSAMstrandFieldIntron = true
		Array[String] outSAMattributes = ["NH","HI","AS","nM"]

		# BAM processing
		Boolean markDup = true
		Int bamRemoveDuplicatesMate2basesN = 0

		# Output Filtering
		Boolean outFilterBySJ = true
		Int outFilterMultimapScoreRange = 1
		Int outFilterMultimapNmax = 20
		Int outFilterMismatchNmax = 10
		Float outFilterMismatchNoverLmax = 0.3
		Float outFilterMismatchNoverReadLmax = 0.04
		Int outFilterScoreMin = 0
		Float outFilterScoreMinOverLread = 0.66
		Int outFilterMatchNmin = 0
		Float outFilterMatchNminOverLread = 0.66
		Boolean outFilterNonCanonical = false
		Boolean removeInconsistentStrands = true

		# Output Filtering: Splice Junctions
		Boolean outSJfilterReadsUniq = false
		Int outSJfilterOverhangMin1 = 30
		Int outSJfilterOverhangMin2 = 12
		Int outSJfilterOverhangMin3 = 12
		Int outSJfilterOverhangMin4 = 12
		Int outSJfilterCountUniqueMin1 = 3
		Int outSJfilterCountUniqueMin2 = 1
		Int outSJfilterCountUniqueMin3 = 1
		Int outSJfilterCountUniqueMin4 = 1
		Int outSJfilterCountTotalMin1 = 3
		Int outSJfilterCountTotalMin2 = 1
		Int outSJfilterCountTotalMin3 = 1
		Int outSJfilterCountTotalMin4 = 1
		Int outSJfilterDistToOtherSJmin1 = 10
		Int outSJfilterDistToOtherSJmin2 = 0
		Int outSJfilterDistToOtherSJmin3 = 5
		Int outSJfilterDistToOtherSJmin4 = 10
		Array[Int] outSJfilterIntronMaxVsReadN = [50000,100000,200000]

		# Scoring
		Int scoreGap = 0
		Int scoreGapNoncan = -8
		Int scoreGapGCAG = -4
		Int scoreGapATAC = -8
		Float scoreGenomicLengthLog2scale = -0.25
		Int scoreDelOpen = -2
		Int scoreDelBase = -2
		Int scoreInsOpen = -2
		Int scoreInsBase = -2
		Int scoreStitchSJshift = 1

		# Alignments and Seeding
		## Seed
		Int seedSearchStartLmax = 50
		Float seedSearchStartLmaxOverLread = 1.0
		Int seedSearchLmax = 0
		Int seedMultimapNmax = 10000
		Int seedPerReadNmax = 1000
		Int seedPerWindowNmax = 50
		Int seedNoneLociPerWindow = 10
		Int seedSplitMin = 12
		Int seedMapMin = 5
		## Alignments
		Int alignIntronMin = 21
		Int alignIntronMax = 0
		Int alignMatesGapMax = 0
		Int alignSJoverhangMin = 5
		Int alignSJstitchMismatchNmax1 = 0
		Int alignSJstitchMismatchNmax2 = -1
		Int alignSJstitchMismatchNmax3 = 0
		Int alignSJstitchMismatchNmax4 = 0
		Int alignSJDBoverhangMin = 3
		Int alignSplicedMateMapLmin = 0
		Float alignSplicedMateMapLminOverLmate = 0.66
		Int alignWindowsPerReadNmax = 10000
		Int alignTranscriptsPerWindowNmax = 100
		Int alignTranscriptsPerReadNmax = 10000
		Boolean alignEndsTypeE2E = false
		Int alignEndsProtrudeMax = 0
		Boolean alignEndsProtrudeConc = true
		Boolean alignSoftClipAtReferenceEnds = true
		Int alignInsertionFlushRight = false

		# Paired-End reads
		Int peOverlapNbasesMin = 0
		Float peOverlapMMp = 0.01

		# Windows, Anchors, Binning
		Int winAnchorMultimapNmax = 50
		Int winBinNbits = 16
		Int winAnchorDistNbins = 9
		Int winFlankNbins = 4

		# Quantification of Annotations
		Boolean quantModeSAM = true
		Int quantTranscriptomeBAMcompression = -1
		Int quantTranscriptomeBanSingleend = false

		# 2-pass Mapping
		Boolean twopassMode = true
		Int twopass1readsN = -1

		# WASP
		Boolean waspOutputMode = false

		Int threads = 1
		Int memoryByThreads = 768
		String memory = "32G"
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(fastqR1),subString,subStringReplace)
	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d ~{outputRep} ]]; then
			mkdir -p ~{outputRep}
		fi

		~{path_exe} --runMode alignReads \
			--runThreadN ~{threads}
			--genomeDir ~{genomeDir} \
			--genomeLoad ~{true="LoadAndKeep" false="NoSharedMemory" genomeLoad} \
			--readFilesIn ~{fastqR1} ~{default="" fastqR2} \
			~{default="" "--sjdbGTFfile " + sjdbGTF} \
			~{default="" "--sjdbFileChrStartEnd " + sjdbBED} \
			~{default="" "--sjdbGTFchrPrefix " + sjdbGTFchrPrefix} \
			--sjdbGTFfeatureExon ~{sjdbGTFfeatureExon} \
			--sjdbGTFtagExonParentTranscript ~{sjdbGTFtagExonParentTranscript} \
			--sjdbGTFtagExonParentGene ~{sjdbGTFtagExonParentGene} \
			--sjdbGTFtagExonParentGeneName ~{sep=" " sjdbGTFtagExonParentGeneName} \
			--sjdbGTFtagExonParentGeneType ~{sep=" " sjdbGTFtagExonParentGeneType} \
			--sjdbOverhang ~{sjdbOverhang} \
			--sjdbScore ~{sjdbScore} \
			--sjdbInsertSave ~{true="All" false="Basic" sjdbInsertSaveAll} \
			--readQualityScoreBase ~{readQualityScoreBase} \
			--clipAdapterType ~{true="Hammer" false="None" clipAdapterHamming} \
			--clip3pNbases ~{sep=" " clip3pNbases} \
			--clip3pAdapterSeq ~{sep=" " clip3pAdapterSeq} \
			--clip3pAdapterMMp ~{sep=" " clip3pAdapterMMp} \
			--clip3pAfterAdapterNbases ~{sep=" " clip3pAfterAdapterNbases} \
			--clip5pNbases ~{sep=" " clip5pNbases} \
			--limitIObufferSize ~{limitIObufferSize} \
			--limitOutSAMoneReadBytes ~{limitOutSAMoneReadBytes} \
			--limitOutSJoneRead ~{limitOutSJoneRead} \
			--limitOutSJcollapsed ~{limitOutSJcollapsed} \
			--limitBAMsortRAM ~{limitBAMsortRAM} \
			--limitSjdbInsertNsj ~{limitSjdbInsertNsj} \
			--limitNreadsSoft ~{limitNreadsSoft} \
			~{true="--outReadsUnmapped" false="" outReadsUnmapped} \
			--outSAMtype ~{true="BAM" false="SAM" outBam} ~{true="SortedByCoordinate" false="Unsorted" sorted} \
			--outSAMmode ~{true="Full" false="NoQS" outSAMmodeQuality} \
			--outSAMstrandField ~{true="IntronMotif" false="None" outSAMstrandFieldIntron} \
			--outSAMattributes ~{sep=" " outSAMattributes} \
			~{true="--bamRemoveDuplicatesType UniqueIdenticalNotMulti" false="" markDup} \
			--bamRemoveDuplicatesMate2basesN ~{bamRemoveDuplicatesMate2basesN} \
			--outFilterType ~{true="BySJout" false="Normal" outFilterBySJ} \
			--outFilterMultimapScoreRange ~{outFilterMultimapScoreRange} \
			--outFilterMultimapNmax ~{outFilterMultimapNmax} \
			--outFilterMismatchNmax ~{outFilterMismatchNmax} \
			--outFilterMismatchNoverLmax ~{outFilterMismatchNoverLmax} \
			--outFilterMismatchNoverReadLmax ~{outFilterMismatchNoverReadLmax} \
			--outFilterScoreMin ~{outFilterScoreMin} \
			--outFilterScoreMinOverLread ~{outFilterScoreMinOverLread} \
			--outFilterMatchNmin ~{outFilterMatchNmin} \
			--outFilterMatchNminOverLread ~{outFilterMatchNminOverLread} \
			--outFilterIntronMotifs ~{true="RemoveNoncanonicalUnannotated" false="None" outFilterNonCanonical} \
			--outFilterIntronStrands ~{true="RemoveInconsistentStrands" false="None" removeInconsistentStrands} \
			--outSJfilterReads ~{true="Unique" false="All" outSJfilterReadsUniq} \
			--outSJfilterOverhangMin ~{outSJfilterOverhangMin1} ~{outSJfilterOverhangMin2} ~{outSJfilterOverhangMin3} ~{outSJfilterOverhangMin4} \
			--outSJfilterCountUniqueMin ~{outSJfilterCountUniqueMin1} ~{outSJfilterCountUniqueMin2} ~{outSJfilterCountUniqueMin3} ~{outSJfilterCountUniqueMin4} \
			--outSJfilterCountTotalMin ~{outSJfilterCountTotalMin1} ~{outSJfilterCountTotalMin2} ~{outSJfilterCountTotalMin3} ~{outSJfilterCountTotalMin4} \
			--outSJfilterDistToOtherSJmin ~{outSJfilterDistToOtherSJmin1} ~{outSJfilterDistToOtherSJmin2} ~{outSJfilterDistToOtherSJmin3} ~{outSJfilterDistToOtherSJmin4} \
			--outSJfilterIntronMaxVsReadN ~{sep=" " outSJfilterIntronMaxVsReadN} \
			--scoreGap ~{scoreGap} \
			--scoreGapNoncan ~{scoreGapNoncan} \
			--scoreGapGCAG ~{scoreGapGCAG} \
			--scoreGapATAC ~{scoreGapATAC} \
			--scoreGenomicLengthLog2scale ~{scoreGenomicLengthLog2scale} \
			--scoreDelOpen ~{scoreDelOpen} \
			--scoreDelBase ~{scoreDelBase} \
			--scoreInsOpen ~{scoreInsOpen} \
			--scoreInsBase ~{scoreInsBase} \
			--scoreStitchSJshift ~{scoreStitchSJshift} \
			--seedSearchStartLmax ~{seedSearchStartLmax} \
			--seedSearchStartLmaxOverLread ~{seedSearchStartLmaxOverLread} \
			--seedSearchLmax ~{seedSearchLmax} \
			--seedMultimapNmax ~{seedMultimapNmax} \
			--seedPerReadNmax ~{seedPerReadNmax} \
			--seedPerWindowNmax ~{seedPerWindowNmax} \
			--seedNoneLociPerWindow ~{seedNoneLociPerWindow} \
			--seedSplitMin ~{seedSplitMin} \
			--seedMapMin ~{seedMapMin} \
			--alignIntronMin ~{alignIntronMin} \
			--alignIntronMax ~{alignIntronMax} \
			--alignMatesGapMax ~{alignMatesGapMax} \
			--alignSJoverhangMin ~{alignSJoverhangMin} \
			--alignSJstitchMismatchNmax ~{alignSJstitchMismatchNmax1} ~{alignSJstitchMismatchNmax2} ~{alignSJstitchMismatchNmax3} ~{alignSJstitchMismatchNmax4} \
			--alignSJDBoverhangMin ~{alignSJDBoverhangMin} \
			--alignSplicedMateMapLmin ~{alignSplicedMateMapLmin} \
			--alignSplicedMateMapLminOverLmate ~{alignSplicedMateMapLminOverLmate} \
			--alignWindowsPerReadNmax ~{alignWindowsPerReadNmax} \
			--alignTranscriptsPerWindowNmax ~{alignTranscriptsPerWindowNmax} \
			--alignTranscriptsPerReadNmax ~{alignTranscriptsPerReadNmax} \
			--alignEndsType ~{true="EndToEnd" false="Local" alignEndsTypeE2E} \
			--alignEndsProtrude ~{alignEndsProtrudeMax} ~{true="ConcordantPair" false="DiscordantPair" alignEndsProtrudeConc} \
			--alignSoftClipAtReferenceEnds ~{true="Yes" false="No" alignSoftClipAtReferenceEnds} \
			--alignInsertionFlush ~{true="Right" false="None" alignInsertionFlushRight} \
			--peOverlapNbasesMin ~{peOverlapNbasesMin} \
			--peOverlapMMp ~{peOverlapMMp} \
			--winAnchorMultimapNmax ~{winAnchorMultimapNmax} \
			--winBinNbits ~{winBinNbits} \
			--winAnchorDistNbins ~{winAnchorDistNbins} \
			--winFlankNbins ~{winFlankNbins} \
			--quantMode ~{true="TranscriptomeSAM" false="GeneCounts" quantModeSAM} \
			--quantTranscriptomeBAMcompression ~{quantTranscriptomeBAMcompression} \
			--quantTranscriptomeBan ~{true="Singleend" false="IndelSoftclipSingleend" quantTranscriptomeBanSingleend} \
			--twopassMode ~{true="Basic" false="None" twopassMode} \
			--twopass1readsN ~{twopass1readsN} \
			~{true="--waspOutputMode SAMtag" false="" waspOutputMode} \
			--outFileNamePrefix ~{outputRep} \
			--readFilesCommand ~{readFilesCommand}

	>>>

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "star"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files will be generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),"(\.bam|\.sam|\.cram)","")]',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Substring to remove to get sample name [default: "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		fastqR1: {
			description: 'Input file with reads 1 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		fastqR2: {
			description: 'Input file with reads 2 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Input - not required'
		}
		genomeDir: {
			description: 'Path to the directory where genome files are stored',
			category: 'Required'
		}
		genomeLoad: {
			description: 'Mode of shared memory usage for the genome files (active LoadAndKeep Mode)',
			category: 'Genome Parameters'
		}
		genomeFastaFiles: {
			description: 'Path to the fasta files with the genome sequences (add extra (new) sequences tothe genome)',
			category: 'Genome Parameters'
		}
		sjdbGTF: {
			description: 'Path to the GTF file with annotations',
			category: 'Splice Junctions Database'
		}
		sjdbGTFchrPrefix: {
			description: 'Prefix for chromosome names in a GTF file (e.g. "chr" for usingENSMEBL annotations with UCSC genomes)',
			category: 'Splice Junctions Database'
		}
		sjdbGTFfeatureExon: {
			description: 'Feature type in GTF file to be used as exons for building transcripts',
			category: 'Splice Junctions Database'
		}
		sjdbGTFtagExonParentTranscript: {
			description: 'GTF attribute name for parent transcript ID (default ”transcriptid”works for GTF files)',
			category: 'Splice Junctions Database'
		}
		sjdbGTFtagExonParentGene: {
			description: 'GTF attribute name for parent gene ID (default ”geneid” works forGTF files)',
			category: 'Splice Junctions Database'
		}
		sjdbGTFtagExonParentGeneName: {
			description: 'GTF attribute name for parent gene name',
			category: 'Splice Junctions Database'
		}
		sjdbGTFtagExonParentGeneType: {
			description: 'GTF attribute name for parent gene type',
			category: 'Splice Junctions Database'
		}
		sjdbOverhang: {
			description: 'Length of the donor/acceptor sequence on each side of the junctions, ideally = (mate length - 1)',
			category: 'Splice Junctions Database'
		}
		sjdbScore: {
			description: 'Extra alignment score for alignments that cross database junctions',
			category: 'Splice Junctions Database'
		}
		sjdbInsertSaveAll: {
			description: 'Save all files when sjdb junctions are inserted on the fly at the mapping step',
			category: 'Splice Junctions Database'
		}
		sjdbBED: {
			description: 'Path to the file with genomic coordinates or the splice junction introns.',
			category: 'Splice Junctions Database'
		}
		readFilesCommand: {
			description: 'Command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout',
			category: 'Read Parameters'
		}
		readQualityScoreBase: {
			description: 'Number to be subtracted from the ASCII code to get Phred quality score',
			category: 'Read Parameters'
		}
		clipAdapterHamming: {
			description: 'Adapter clipping type in Hamming mode',
			category: 'Read Clipping'
		}
		clip3pNbases: {
			description: 'Number(s) of bases to clip from 3p of each mate. If one value is given, it will be assumed the same for both mates.',
			category: 'Read Clipping'
		}
		clip3pAdapterSeq: {
			description: 'Adapter sequences to clip from 3p of each mate. If one value is given, it will be assumed the same for both mates.',
			category: 'Read Clipping'
		}
		clip3pAdapterMMp: {
			description: 'Max proportion of mismatches for 3p adapter clipping for each mate. If one value is given, it will be assumed the same for both mates.',
			category: 'Read Clipping'
		}
		clip3pAfterAdapterNbases: {
			description: 'Number of bases to clip from 3p of each mate after the adapter clipping. If one value is given, it will be assumed the same for both mates.',
			category: 'Read Clipping'
		}
		clip5pNbases: {
			description: 'Number(s) of bases to clip from 5p of each mate. If one value is given, it will be assumed the same for both mates.',
			category: 'Read Clipping'
		}
		limitIObufferSize: {
			description: 'Max available buffers size (bytes) for input/output, per thread',
			category: 'Limits'
		}
		limitOutSAMoneReadBytes: {
			description: 'Max size of the SAM record (bytes) for one read. Recommended value: >(2*(LengthMate1+LengthMate2+100)*outFilterMultimapNmax',
			category: 'Limits'
		}
		limitOutSJoneRead: {
			description: 'Max number of junctions for one read (including all multi-mappers)',
			category: 'Limits'
		}
		limitOutSJcollapsed: {
			description: 'Max number of collapsed junctions',
			category: 'Limits'
		}
		limitBAMsortRAM: {
			description: 'Maximum available RAM (bytes) for sorting BAM.',
			category: 'Limits'
		}
		limitSjdbInsertNsj: {
			description: 'Maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run',
			category: 'Limits'
		}
		limitNreadsSoft: {
			description: 'Soft limit on the number of reads',
			category: 'Limits'
		}
		outReadsUnmapped: {
			description: 'Output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s).',
			category: 'Output: general'
		}
		outBam: {
			description: 'Output in bam or sam file',
			category: 'Output: SAM and BAM'
		}
		sorted: {
			description: 'Sorted output file.',
			category: 'Output: SAM and BAM'
		}
		outSAMmodeQuality: {
			description: 'Sets the quality mode for SAM/BAM output',
			category: 'Output: SAM and BAM'
		}
		outSAMstrandFieldIntron: {
			description: 'Cufflinks-like strand field flag.',
			category: 'Output: SAM and BAM'
		}
		outSAMattributes: {
			description: 'Desired SAM attributes, in the order desired for the output SAM.',
			category: 'Output: SAM and BAM'
		}
		markDup: {
			description: 'Mark duplicates in the BAM file, for now only works with (i) sorted BAM fed with inputBAMfile, and (ii) for paired-end alignments only',
			category: 'BAM processing'
		}
		bamRemoveDuplicatesMate2basesN: {
			description: 'Number of bases from the 5’ of mate 2 to use in collapsing',
			category: 'BAM processing'
		}
		outFilterBySJ: {
			description: 'Keep only those reads that contain junctions that passed filtering into SJ.out.tab',
			category: 'Output filtering'
		}
		outFilterMultimapScoreRange: {
			description: 'Score range below the maximum score for multimapping alignments',
			category: 'Output filtering'
		}
		outFilterMultimapNmax: {
			description: 'Maximum number of loci the read is allowed to map to.',
			category: 'Output filtering'
		}
		outFilterMismatchNmax: {
			description: 'Alignment will be output only if it has no more mismatches than this value.',
			category: 'Output filtering'
		}
		outFilterMismatchNoverLmax: {
			description: 'Alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value.',
			category: 'Output filtering'
		}
		outFilterMismatchNoverReadLmax: {
			description: 'Alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value.',
			category: 'Output filtering'
		}
		outFilterScoreMin: {
			description: 'Alignment will be output only if its score is higher than or equal to this',
			category: 'Output filtering'
		}
		outFilterScoreMinOverLread: {
			description: 'Same as outFilterScoreMin, but normalized to read length (sum of mates’ lengths for paired-end reads)',
			category: 'Output filtering'
		}
		outFilterMatchNmin: {
			description: 'Alignment will be output only if the number of matched bases is higher than or equal to this value.',
			category: 'Output filtering'
		}
		outFilterMatchNminOverLread: {
			description: 'Same as outFilterMatchNmin, but normalized to the read length (sum of mates’ lengths for paired-end reads).',
			category: 'Output filtering'
		}
		outFilterNonCanonical: {
			description: 'Filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. The annotated non-canonical junctions will be kept.',
			category: 'Output filtering'
		}
		removeInconsistentStrands: {
			description: 'Remove alignments that have junctions with inconsistent strands.',
			category: 'Output filtering'
		}
		outSJfilterReadsUniq: {
			description: 'Uniquely mapping reads consider for collapsed splice junctions output.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterOverhangMin1: {
			description: 'Minimum overhang length for splice junctions on both sides for non-canonical motifs.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterOverhangMin2: {
			description: 'Minimum overhang length for splice junctions on both sides for GT/AG and CT/AC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterOverhangMin3: {
			description: 'Minimum overhang length for splice junctions on both sides for GC/AG and CT/GC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterOverhangMin4: {
			description: 'Minimum overhang length for splice junctions on both sides for AT/AC and GT/AT motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountUniqueMin1: {
			description: 'Minimum uniquely mapping read count per junction for non-canonical motifs',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountUniqueMin2: {
			description: 'Minimum uniquely mapping read count per junction for GT/AG and CT/AC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountUniqueMin3: {
			description: 'Minimum uniquely mapping read count per junction for GC/AG and CT/GC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountUniqueMin4: {
			description: 'Minimum uniquely mapping read count per junction for AT/AC and GT/AT motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountTotalMin1: {
			description: 'Minimum total (multi-mapping+unique) read count per junction for non-canonical motifs',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountTotalMin2: {
			description: 'Minimum total (multi-mapping+unique) read count per junction for GT/AG and CT/AC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountTotalMin3: {
			description: 'Minimum total (multi-mapping+unique) read count per junction for GC/AG and CT/GC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterCountTotalMin4: {
			description: 'Minimum total (multi-mapping+unique) read count per junction for AT/AC and GT/AT motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterDistToOtherSJmin1: {
			description: 'Minimum allowed distance to other junctions’ donor/acceptor for non-canonical motifs',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterDistToOtherSJmin2: {
			description: 'Minimum allowed distance to other junctions’ donor/acceptor for GT/AG and CT/AC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterDistToOtherSJmin3: {
			description: 'Minimum allowed distance to other junctions’ donor/acceptor for GC/AG and CT/GC motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterDistToOtherSJmin4: {
			description: 'Minimum allowed distance to other junctions’ donor/acceptor for AT/AC and GT/AT motif.',
			category: 'Output filtering (Splice Junctions)'
		}
		outSJfilterIntronMaxVsReadN: {
			description: 'Maximum gap allowed for junctions supported by 1,2,3,,,N reads',
			category: 'Output filtering (Splice Junctions)'
		}
		scoreGap: {
			description: 'Splice junction penalty (independent on intron motif)',
			category: 'Scoring'
		}
		scoreGapNoncan: {
			description: 'non-canonical junction penalty (in addition to scoreGap)',
			category: 'Scoring'
		}
		scoreGapGCAG: {
			description: 'GC/AG and CT/GC junction penalty (in addition to scoreGap)',
			category: 'Scoring'
		}
		scoreGapATAC: {
			description: 'AT/AC and GT/AT junction penalty (in addition to scoreGap)',
			category: 'Scoring'
		}
		scoreGenomicLengthLog2scale: {
			description: 'Extra score logarithmically scaled with genomic length of the alignment: scoreGenomicLengthLog2scale*log2(genomicLength)',
			category: 'Scoring'
		}
		scoreDelOpen: {
			description: 'Deletion open penalty',
			category: 'Scoring'
		}
		scoreDelBase: {
			description: 'Deletion extension penalty per base (in addition to scoreDelOpen)',
			category: 'Scoring'
		}
		scoreInsOpen: {
			description: 'Insertion open penalty',
			category: 'Scoring'
		}
		scoreInsBase: {
			description: 'Insertion extension penalty per base (in addition to scoreInsOpen)',
			category: 'Scoring'
		}
		scoreStitchSJshift: {
			description: 'Maximum score reduction while searching for SJ boundaries in the stitching step',
			category: 'Scoring'
		}
		seedSearchStartLmax: {
			description: 'Defines the search start point through the read - the read is split into pieces no longer than this value',
			category: 'Alignments and Seeding'
		}
		seedSearchStartLmaxOverLread: {
			description: 'SeedSearchStartLmax normalized to read length (sum of mates’ lengths for paired-end reads)',
			category: 'Alignments and Seeding'
		}
		seedSearchLmax: {
			description: 'Defines the maximum length of the seeds, if =0 seed length is not limited',
			category: 'Alignments and Seeding'
		}
		seedMultimapNmax: {
			description: 'Only pieces that map fewer than this value are utilized in the stitching procedure',
			category: 'Alignments and Seeding'
		}
		seedPerReadNmax: {
			description: 'Max number of seeds per read',
			category: 'Alignments and Seeding'
		}
		seedPerWindowNmax: {
			description: 'Max number of seeds per window',
			category: 'Alignments and Seeding'
		}
		seedNoneLociPerWindow: {
			description: 'Max number of one seed loci per window',
			category: 'Alignments and Seeding'
		}
		seedSplitMin: {
			description: 'Min length of the seed sequences split by Ns or mate gap',
			category: 'Alignments and Seeding'
		}
		seedMapMin: {
			description: 'Min length of seeds to be mapped',
			category: 'Alignments and Seeding'
		}
		alignIntronMin: {
			description: 'Minimum intron size: genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion',
			category: 'Alignments and Seeding'
		}
		alignIntronMax: {
			description: 'Maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins',
			category: 'Alignments and Seeding'
		}
		alignMatesGapMax: {
			description: 'Maximum gap between two mates, if 0, max intron gap will be determined by (2ˆwinBinNbits)*winAnchorDistNbins',
			category: 'Alignments and Seeding'
		}
		alignSJoverhangMin: {
			description: 'Minimum overhang (i.e. block size) for spliced alignments',
			category: 'Alignments and Seeding'
		}
		alignSJstitchMismatchNmax1: {
			description: 'Maximum number of mismatches for stitching of the splice junctions (-1: no limit) for non-canonical motifs.',
			category: 'Alignments and Seeding'
		}
		alignSJstitchMismatchNmax2: {
			description: 'Maximum number of mismatches for stitching of the splice junctions (-1: no limit) for GT/AG and CT/AC motif.',
			category: 'Alignments and Seeding'
		}
		alignSJstitchMismatchNmax3: {
			description: 'Maximum number of mismatches for stitching of the splice junctions (-1: no limit) for GC/AG and CT/GC motif.',
			category: 'Alignments and Seeding'
		}
		alignSJstitchMismatchNmax4: {
			description: 'Maximum number of mismatches for stitching of the splice junctions (-1: no limit) for AT/AC and GT/AT motif.',
			category: 'Alignments and Seeding'
		}
		alignSJDBoverhangMin: {
			description: 'Minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments',
			category: 'Alignments and Seeding'
		}
		alignSplicedMateMapLmin: {
			description: 'Minimum mapped length for a read mate that is spliced',
			category: 'Alignments and Seeding'
		}
		alignSplicedMateMapLminOverLmate: {
			description: 'AlignSplicedMateMapLmin normalized to mate length',
			category: 'Alignments and Seeding'
		}
		alignWindowsPerReadNmax: {
			description: 'Max number of windows per read',
			category: 'Alignments and Seeding'
		}
		alignTranscriptsPerWindowNmax: {
			description: 'Max number of transcripts per window.',
			category: 'Alignments and Seeding'
		}
		alignTranscriptsPerReadNmax: {
			description: 'Max number of different alignments per read to consider.',
			category: 'Alignments and Seeding'
		}
		alignEndsTypeE2E: {
			description: 'Set type of read ends alignment to EndToEnd mode.',
			category: 'Alignments and Seeding'
		}
		alignEndsProtrudeMax: {
			description: 'Maximum number of protrusion bases allowed.',
			category: 'Alignments and Seeding'
		}
		alignEndsProtrudeConc: {
			description: 'Set report alignments with non-zero protrusion as concordant pairs.',
			category: 'Alignments and Seeding'
		}
		alignSoftClipAtReferenceEnds: {
			description: 'Allow the soft-clipping of the alignments past the end of the chromosomes',
			category: 'Alignments and Seeding'
		}
		alignInsertionFlushRight: {
			description: 'Insertions are flushed to the right',
			category: 'Alignments and Seeding'
		}
		peOverlapNbasesMin: {
			description: 'Minimum number of overlap bases to trigger mates merging and realignment',
			category: 'Paired-End reads'
		}
		peOverlapMMp: {
			description: 'Maximum proportion of mismatched bases in the overlap area',
			category: 'Paired-End reads'
		}
		winAnchorMultimapNmax: {
			description: 'max number of loci anchors are allowed to map to.',
			category: 'Windows, Anchors, Binning'
		}
		winBinNbits: {
			description: 'log2(winBin), where winBin is the size of the bin for the windows/clustering, each window will occupy an integer number of bins.',
			category: 'Windows, Anchors, Binning'
		}
		winAnchorDistNbins: {
			description: 'max number of bins between two anchors that allows aggregation of anchors into one window.',
			category: 'Windows, Anchors, Binning'
		}
		winFlankNbins: {
			description: 'og2(winFlank), where win Flank is the size of the left and right flanking regions for each window.',
			category: 'Windows, Anchors, Binning'
		}
		quantModeSAM: {
			description: 'Set quantMode to sam mode',
			category: 'Quantification of Annotations'
		}
		quantTranscriptomeBAMcompression: {
			description: 'Transcriptome BAM compression level.',
			category: 'Quantification of Annotations'
		}
		quantTranscriptomeBanSingleend: {
			description: 'prohibit indels, soft clipping and single-end alignments - compatible with RSEM.',
			category: 'Quantification of Annotations'
		}
		twopassMode: {
			description: 'Active 2-pass mapping mode.',
			category: '2-pass Mapping'
		}
		twopass1readsN: {
			description: 'number of reads to process for the 1st step. Use very large number (or default -1) to map all reads in the first step.',
			category: '2-pass Mapping'
		}
		waspOutputMode: {
			description: 'WASP allele-specific output type. https://www.nature.com/articles/nmeth.3582.',
			category: 'WASP'
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
