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

task alignReads {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-05"
	}

	input {
		String path_exe = "STAR"

		File fastqR1
		File? fastqR2

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
		Int readQualityScoreBase = 33

		# Read Clipping
		Boolean clipAdapterHamming = true
		Int clip3pNbases = 0
		Boolean clip3pAdapterSeq = false
		Float clip3pAdapterMMp = 0.1
		Int clip3pAfterAdapterNbases = 0
		Int clip5pNbases = 0

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

		# Output Filtering:  Splice Junctions
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
		Int alignIntronMax = O
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
		Float winReadCoverageRelativeMin = 0.5
		Int winReadCoverageBasesMin = 0

		# Quantification of Annotations
		Boolean quantModeSAM = true
		Int quantTranscriptomeBAMcompression = -1
		Int quantTranscriptomeBanSingleend = false

		# 2-pass Mapping
		Boolean twopassMode = true
		Int twopass1readsN = -1

		# WASP
		Boolean waspOutputMode = false

		String readFilesCommand = "zcat"

		String? outputPath
		String? name
		String subString = "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"
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

	String baseName = if defined(sample) then sample else sub(basename(fastqR1),subString,subStringReplace)
	String outputRep = if defined(outputPath) then "~{outputPath}/~{baseName}" else "~{baseName}"

	command <<<

		if [[ ! -d ~{outputRep} ]]; then
			mkdir -p ~{outputRep}
		fi

		~{path_exe} --runMode alignReads \
			--runThreadN ~{threads}
			--genomeDir ~{genomeDir} \
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
			--clip3pNbases ~{clip3pNbases} \
			~{true="--clip3pAdapterSeq" false="" clip3pAdapterSeq} \
			--clip3pAdapterMMp ~{clip3pAdapterMMp} \
			--clip3pAfterAdapterNbases ~{clip3pAfterAdapterNbases} \
			--clip5pNbases ~{clip5pNbases} \
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
			--winReadCoverageRelativeMin ~{winReadCoverageRelativeMin} \
			--winReadCoverageBasesMin ~{winReadCoverageBasesMin} \
			--quantMode ~{true="TranscriptomeSAM" false="GeneCounts" quantModeSAM} \
			--quantTranscriptomeBAMcompression ~{quantTranscriptomeBAMcompression} \
			--quantTranscriptomeBan ~{true="Singleend" false="IndelSoftclipSingleend" quantTranscriptomeBanSingleend} \
			--twopassMode ~{true="Basic" false="None" twopassMode} \
			--twopass1readsN ~{twopass1readsN} \
			~{true="--waspOutputMode SAMtag" false="" waspOutputMode} \
			--outFileNamePrefix ~{outPrefix} \
			--readFilesCommand ~{readFilesCommand}

	>>>

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
