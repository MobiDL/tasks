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
		String path_exe = "vep"

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
		~{path_exe} --help | grep -A5 "Versions"
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
			description: 'Path used as executable [default: "vep"]',
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

task vep_cache {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-09-29"
	}

	input {
		String path_exe = "vep"

		File in
		String? outputPath
		String? name
		String suffix = ".vep"
		String subString = ".vcf(.gz)?"
		String subStringReplace = ""

		## Basic options
		String species = "homo_sapiens"
		String assembly = "GRCh37"
		Boolean overwrite = false
		Boolean statsOutput = true

		## Cache options
		String dirBaseCache = "$HOME/.vep/"
		String dirCache = "$HOME/.vep/"
		String dirPlugins = "$HOME/.vep/"
		Boolean offlineMode = true
		File refFasta
		File refFai
		Boolean mergedCache = true
		String? cacheVersion
		Int bufferSize = 5000

		## Other annotations sources
		Array[String]? plugins
		Array[String]? customs
		Array[String]? gff
		Array[String]? gtf

		## Output format
		String outputFormat = "vcf"
		Boolean compressGzip = true
		Array[String]? fields

		## Output options
		Boolean variantClass = true
		Boolean sift = true
		Boolean polyphen = true
		Boolean humDiv = false
		Boolean nearestTranscript = true
		Int upstreamDist = 5000
		Int downstreamDist = upstreamDist
		Boolean overlaps = true
		Boolean genePhenotype = true
		Boolean regulatory = true
		Boolean showRefAllele = true
		Boolean totalLength = true
		Boolean numbers = true
		Boolean escapeHGVS = true
		Boolean overwriteCSQ = true
		String vcfInfoField = "CSQ"
		String terms = "SO"
		Boolean header = true
		Boolean shift3prime = false
		Boolean shiftGenomic = false
		Boolean shiftLength = false

		## Identifiers
		Boolean hgvs = true
		Boolean hgvsg = true
		Boolean hgvsShift = true
		Boolean transcriptEnsemblIDVersion = true
		Boolean proteinEnsemblID = true
		Boolean geneSymbol = true
		Boolean ccds = true
		Boolean uniprot = true
		Boolean canonical = true
		Boolean biotype = true
		Boolean domain = true
		Boolean xRefSeq = true
		Boolean tsl = true
		Boolean appris = true
		Boolean mane = true
		File? chrSynonyms

		## Co-located
		Boolean checkExisting = true
		Boolean clinsigAllele = true
		Boolean keepNullAllele = true
		Boolean checkAlleles = true
		Boolean AF = true
		Boolean maxAF = true
		Boolean AF1000g = true
		Boolean AFEsp = true
		Boolean AFGnomad = true
		Boolean AFExac = false
		Boolean pubmed = true
		Boolean varSynonyms = true
		Boolean includeFailed = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(in),subString,subStringReplace)
	String ext = if compressGzip then ".vcf.gz" else ".vcf"
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}~{ext}" else "~{baseName}~{suffix}~{ext}"
	String outputFileStats = if defined(outputPath) then "~{outputPath}/~{baseName}~{suffix}.stats.html" else "~{baseName}~{suffix}.stats.html"

	Boolean pluginsDefined = defined(plugins)
	Boolean customsDefined = defined(customs)
	Boolean gffDefined = defined(gff)
	Boolean gtfDefined = defined(gtf)
	Boolean fieldsDefined = defined(fields)

	command <<<

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} \
			--species ~{species} \
			--assembly ~{assembly} \
			--input_file ~{in} \
			--output_file ~{outputFile} \
			~{true="--force_overwrite" false="" overwrite} \
			--stats_file ~{outputFileStats} \
			~{true="" false="--no_stats" statsOutput} \
			--fork ~{threads} \
			--cache \
			--dir ~{dirBaseCache} \
			--dir_cache ~{dirCache} \
			--dir_plugins ~{dirPlugins} \
			~{true="--offline" false="" offlineMode} \
			--fasta ~{refFasta} \
			~{true="--merged" false="--refseq" mergedCache} \
			~{default="" "--cache_version " + cacheVersion} \
			--buffer_size ~{bufferSize} \
			~{true="--plugin " false="" pluginsDefined}~{sep="--plugin " plugins} \
			~{true="--custom " false="" customsDefined}~{sep="--custom " customs} \
			~{true="--gff " false="" gffDefined}~{sep="--gff " gff} \
			~{true="--gtf " false="" gtfDefined}~{sep="--gtf " gtf} \
			--~{outputFormat} \
			~{true="--compress gzip" false="" compressGzip} \
			~{true="--fields " false="" fieldsDefined}~{sep="," fields} \
			~{true="--variant_class" false="" variantClass} \
			~{true="--sift b" false="" sift} \
			~{true="--polyphen b" false="" polyphen} \
			~{true="--humdiv" false="" humDiv} \
			~{true="--nearest transcript" false="" nearestTranscript} \
			--distance ~{upstreamDist},~{downstreamDist} \
			~{true="--overlaps" false="" overlaps} \
			~{true="--gene_phenotype" false="" genePhenotype} \
			~{true="--regulatory" false="" regulatory} \
			~{true="--show_ref_allele" false="" showRefAllele} \
			~{true="--total_length" false="" totalLength} \
			~{true="--numbers" false="" numbers} \
			~{true="" false="--no_escape" escapeHGVS} \
			~{true="" false="--keep_csq" overwriteCSQ} \
			--vcf_info_field ~{vcfInfoField} \
			--terms ~{terms} \
			~{true="" false="--no_headers" header} \
			--shift_3prime ~{true="1" false="1" shift3prime} \
			--shift_genomic ~{true="1" false="1" shiftGenomic} \
			~{true="--shift_length" false="" shiftLength} \
			~{true="--hgvs" false="" hgvs} \
			~{true="--hgvsg" false="" hgvsg} \
			--shift_hgvs ~{true="1" false="0" hgvsShift} \
			~{true="--transcript_version" false="" transcriptEnsemblIDVersion} \
			~{true="--protein" false="" proteinEnsemblID} \
			~{true="--symbol" false="" geneSymbol} \
			~{true="--ccds" false="" ccds} \
			~{true="--uniprot" false="" uniprot} \
			~{true="--tsl" false="" tsl} \
			~{true="--appris" false="" appris} \
			~{true="--canonical" false="" canonical} \
			~{true="--mane" false="" mane} \
			~{true="--biotype" false="" biotype} \
			~{true="--domains" false="" domain} \
			~{true="--xref_refseq" false="" xRefSeq} \
			~{default="" "--synonyms " + chrSynonyms} \
			~{true="--check_existing" false="" checkExisting} \
			--clin_sig_allele ~{true="1" false="0" clinsigAllele} \
			~{true="" false="--exclude_null_alleles" keepNullAllele} \
			~{true="" false="--no_check_alleles" checkAlleles} \
			~{true="--af" false="" AF} \
			~{true="--max_af" false="" maxAF} \
			~{true="--af_1kg" false="" AF1000g} \
			~{true="--af_esp" false="" AFEsp} \
			~{true="--af_gnomad" false="" AFGnomad} \
			~{true="--af_exac" false="" AFExac} \
			~{true="--pubmed" false="" pubmed} \
			~{true="--var_synonyms" false="" varSynonyms} \
			--failed ~{true="1" false="0" includeFailed}

	>>>

	output {
		File outputFile = outputFile
		File? outputFileStats = outputFileStats
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "cutadapt"]',
			category: 'System'
		}
		outputPath: {
			description: 'Output path where files were generated. [default: pwd()]',
			category: 'Output path/name option'
		}
		name: {
			description: 'Name to use for output file name [default: sub(basename(in),subString,subStringReplace)]',
			category: 'Output path/name option'
		}
		in: {
			description: 'Input reads (format: fastq, fastq.gz)',
			category: 'Required'
		}
		suffix: {
			description: 'Suffix to add to the output [default: .adaptersTrim]',
			category: 'Output path/name option'
		}
		subString: {
			description: 'Extension to remove from the input file [default: "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		species: {
			description: 'Species for your data. This can be the latin name e.g. "homo_sapiens" or any Ensembl alias e.g. "mouse". [default: "homo_sapiens"]',
			category: 'Tool option (basics)'
		}
		assembly: {
			description: 'Select the assembly version to use if more than one available. [default: "GRCh37"]',
			category: 'Tool option (basics)'
		}
		overwrite: {
			description: 'Force the overwrite of the existing file. By default, VEP will fail with an error if the output file already exists. [default: false]',
			category: 'Tool option (basics)'
		}
		statsOutput: {
			description: 'Create an HTML file containing a summary of the VEP run. [default: true]',
			category: 'Tool option (basics)'
		}
		dirBaseCache: {
			description: 'Specify the base cache/plugin directory to use. [default: "$HOME/.vep/"]',
			category: 'Tool option (cache)'
		}
		dirCache: {
			description: 'Specify the cache directory to use. [default: "$HOME/.vep/"]',
			category: 'Tool option (cache)'
		}
		dirPlugins: {
			description: 'Specify the plugin directory to use. [default: "$HOME/.vep/"]',
			category: 'Tool option (cache)'
		}
		offlineMode: {
			description: 'Enable offline mode. [default: true]',
			category: 'Tool option (cache)'
		}
		refFasta: {
			description: 'Reference in fasta format',
			category: 'Required'
		}
		refFai: {
			description: 'Path to the reference file index (format: fai)',
			category: 'Required'
		}
		mergedCache: {
			description: 'Use the merged Ensembl and RefSeq cache, or only RefSeq cache if false. [default: true]',
			category: 'Tool option (cache)'
		}
		cacheVersion: {
			description: 'Use a different cache version than the assumed default (the VEP version).',
			category: 'Tool option (cache)'
		}
		bufferSize: {
			description: 'Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously. [default: 5000]',
			category: 'Tool option (cache)'
		}
		plugins: {
			description: 'Plugins to use (needs to be installed).',
			category: 'Tool option (others annotations)'
		}
		customs: {
			description: 'Custom annotations files (BigWig or tabix indexed)',
			category: 'Tool option (others annotations)'
		}
		gff: {
			description: 'Use GFF transcript annotations as an annotation source.',
			category: 'Tool option (others annotations)'
		}
		gtf: {
			description: 'Use GTF transcript annotations as an annotation source.',
			category: 'Tool option (others annotations)'
		}
		outputFormat: {
			description: 'Specify the output format (choices: vcf, tab or json) [default: "vcf"]',
			category: 'Tool option (output format)'
		}
		compressGzip: {
			description: 'Writes output compressed using either gzip. [default: true]',
			category: 'Tool option (output format)'
		}
		fields: {
			description: 'Configure the output format using a comma separated list of fields.',
			category: 'Tool option (output format)'
		}
		variantClass : {
			description: 'Output the Sequence Ontology variant class. [default: true]',
			category: 'Tool option (output options)'
		}
		sift : {
			description: 'Output SIFT prediction term and score. [default: true]',
			category: 'Tool option (output options)'
		}
		polyphen : {
			description: 'Output POLYPHEN prediction term and score [default: true]',
			category: 'Tool option (output options)'
		}
		humDiv : {
			description: 'Retrieve the humDiv PolyPhen prediction instead of the default humVar. [default: false]',
			category: 'Tool option (output options)'
		}
		nearestTranscript : {
			description: 'Retrieve the transcript with the nearest protein-coding transcription start site (TSS) to each input variant. [default: true]',
			category: 'Tool option (output options)'
		}
		upstreamDist : {
			description: 'Modify the distance upstream (and downstream by default) between a variant and a transcript for which VEP will assign the upstream_gene_variant (and downstream_gene_variant) consequences. [default: 5000]',
			category: 'Tool option (output options)'
		}
		downstreamDist : {
			description: 'Modify the distance downstream between a variant and a transcript for which VEP will assign the downstream_gene_variant consequences. [default: upstreamDist]',
			category: 'Tool option (output options)'
		}
		overlaps: {
			description: 'Report the proportion and length of a transcript overlapped by a structural variant in VCF format. [default: true]',
			category: 'Tool option (output options)'
		}
		genePhenotype : {
			description: 'Indicates if the overlapped gene is associated with a phenotype, disease or trait. [default: true]',
			category: 'Tool option (output options)'
		}
		regulatory : {
			description: 'Look for overlaps with regulatory regions. [default: true]',
			category: 'Tool option (output options)'
		}
		showRefAllele : {
			description: 'Adds the reference allele in the output. [default: true]',
			category: 'Tool option (output options)'
		}
		totalLength : {
			description: 'Give cDNA, CDS and protein positions as Position/Length. [default: true]',
			category: 'Tool option (output options)'
		}
		numbers : {
			description: 'Adds affected exon and intron numbering to to output. [default: true]',
			category: 'Tool option (output options)'
		}
		escapeHGVS : {
			description: 'URI escape HGVS strings. [default: true]',
			category: 'Tool option (output options)'
		}
		overwriteCSQ : {
			description: 'Overwrite existing CSQ entry in VCF INFO field. [default: true]',
			category: 'Tool option (output options)'
		}
		vcfInfoField : {
			description: 'Change the name of the INFO key that VEP write the consequences to in its VCF output. [default: "CSQ"]',
			category: 'Tool option (output options)'
		}
		terms : {
			description: 'The type of consequence terms to output (choices: SO, display, NCBI) [default: "SO"]',
			category: 'Tool option (output options)'
		}
		header : {
			description: 'Write header lines in output files. [default: true]',
			category: 'Tool option (output options)'
		}
		shift3prime : {
			description: 'Right aligns all variants relative to their associated transcripts prior to consequence calculation. [default: false]',
			category: 'Tool option (output options)'
		}
		shiftGenomic : {
			description: 'Right aligns all variants, including intergenic variants, before consequence calculation and updates the Location field. [default: false]',
			category: 'Tool option (output options)'
		}
		shiftLength : {
			description: 'Reports the distance each variant has been shifted when used in conjuction with "shift3prime" [default: false]',
			category: 'Tool option (output options)'
		}
		hgvs: {
			description: 'Add HGVS nomenclature based on Ensembl stable identifiers to the output. [default: true]',
			category: 'Tool option (identifiers)'
		}
		hgvsg: {
			description: 'Add genomic HGVS nomenclature based on the input chromosome name. [default: true]',
			category: 'Tool option (identifiers)'
		}
		hgvsShift: {
			description: 'Enable or disable 3-prime shifting of HGVS notations. [default: true]',
			category: 'Tool option (identifiers)'
		}
		transcriptEnsemblIDVersion: {
			description: 'Add version numbers to Ensembl transcript identifiers. [default: true]',
			category: 'Tool option (identifiers)'
		}
		proteinEnsemblID: {
			description: 'Add the Ensembl protein identifier to the output where appropriate. [default: true]',
			category: 'Tool option (identifiers)'
		}
		geneSymbol: {
			description: 'Adds the gene symbol (e.g. HGNC) (where available) to the output. [default: true]',
			category: 'Tool option (identifiers)'
		}
		ccds: {
			description: 'Adds the CCDS transcript identifer (where available) to the output. [default: true]',
			category: 'Tool option (identifiers)'
		}
		uniprot: {
			description: 'Adds best match accessions for translated protein products from three UniProt-related databases (SWISSPROT, TREMBL and UniParc) to the output. [default: true]',
			category: 'Tool option (identifiers)'
		}
		canonical: {
			description: 'Adds a flag indicating if the transcript is the canonical transcript for the gene. [default: true]',
			category: 'Tool option (identifiers)'
		}
		biotype: {
			description: 'Adds the biotype of the transcript or regulatory feature. [default: true]',
			category: 'Tool option (identifiers)'
		}
		domain: {
			description: 'Adds names of overlapping protein domains to output. [default: true]',
			category: 'Tool option (identifiers)'
		}
		xRefSeq: {
			description: 'Output aligned RefSeq mRNA identifier for transcript. [default: true]',
			category: 'Tool option (identifiers)'
		}
		tsl: {
			description: 'Adds the transcript support level for this transcript to the output (only GRCh38). [default: true]',
			category: 'Tool option (identifiers)'
		}
		appris: {
			description: 'Adds the APPRIS isoform annotation for this transcript to the output (only GRCh38). [default: true]',
			category: 'Tool option (identifiers)'
		}
		mane: {
			description: 'Adds a flag indicating if the transcript is the MANE Select transcript for the gene (only GRCh38). [default: true]',
			category: 'Tool option (identifiers)'
		}
		chrSynonyms		: {
			description: 'Load a file of chromosome synonyms.',
			category: 'Tool option (identifiers)'
		}
		checkExisting: {
			description: 'Checks for the existence of known variants that are co-located with your input. [default: true]',
			category: 'Tool option (co-located)'
		}
		keepNullAllele: {
			description: 'Include variants with unknown alleles when checking for co-located variants. [default: true]',
			category: 'Tool option (co-located)'
		}
		checkAlleles: {
			description: 'Compare using coordinates alone. [default: true]',
			category: 'Tool option (co-located)'
		}
		clinsigAllele: {
			description: 'Return allele specific clinical significance. (true: allele specific; false: locus specific) [default: true]',
			category: 'Tool option (co-located)'
		}
		AF: {
			description: 'Add the global allele frequency (AF) from 1000 Genomes Phase 3 data. [default: true]',
			category: 'Tool option (co-located)'
		}
		maxAF: {
			description: 'Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD. [default: true]',
			category: 'Tool option (co-located)'
		}
		AF1000g: {
			description: 'Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3. [default: true]',
			category: 'Tool option (co-located)'
		}
		AFEsp: {
			description: 'Include allele frequency from NHLBI-ESP populations. [default: true]',
			category: 'Tool option (co-located)'
		}
		AFGnomad: {
			description: 'Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. [default: true]',
			category: 'Tool option (co-located)'
		}
		AFExac: {
			description: 'Include allele frequency from ExAC project populations. [default: false]',
			category: 'Tool option (co-located)'
		}
		pubmed: {
			description: 'Report Pubmed IDs for publications that cite existing variant. [default: true]',
			category: 'Tool option (co-located)'
		}
		varSynonyms: {
			description: 'Report known synonyms for colocated variants. [default: true]',
			category: 'Tool option (co-located)'
		}
		includeFailed: {
			description: 'Include variants tags as failed [default: false]',
			category: 'Tool option (co-located)'
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

task install {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2020-12-02"
	}

	input {
		String path_exe = "vep_install"

		String destDir = "./"
		String cacheVersion = "101"
		String? cacheDir
		String? pluginsDir

		Boolean refseq = true
		Boolean ensembl = true

		Boolean api = true
		Boolean cache = true
		Boolean fasta = true
		Boolean plugins = false

		Boolean quiet = false

		String species = "homo_sapiens"
		String assembly = "GRCh37"
		Array[String]? pluginsList

		Boolean prefer_bin = false
		Boolean checkHTSLIB = false
		Boolean checkBioPerl = false

		Boolean convert = false

		String? cacheURL
		String? fastaURL

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	Boolean merged = if refseq && ensembl then true else false
	String species_complete = if merged then "~{species}_merged" else if refseq then "~{species}_refseq" else species

	String cacheDirFinal = if defined(cacheDir) then "~{cacheDir}" else "~{destDir}/.vep/"
	String pluginsDirFinal = if defined(pluginsDir) then "~{pluginsDir}" else "~{destDir}/.vep/Plugins/"

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	command <<<
		~{path_exe} \
			--NO_UPDATE \
			--DESTDIR ~{destDir} \
			--CACHE_VERSION ~{cacheVersion} \
			--CACHEDIR ~{cacheDirFinal} \
			--PLUGINSDIR ~{pluginsDirFinal} \
			--AUTO ~{true="a" false="l" api}~{true="c" false="" cache}~{true="f" false="" fasta}~{true="p" false="" plugins} \
			--SPECIES ~{species_complete} \
			--ASSEMBLY ~{assembly} \
			~{true="--QUIET" false="" quiet} \
			~{true="--PREFER_BIN" false="" prefer_bin} \
			~{true="" false="--NO_HTSLIB" checkHTSLIB} \
			~{true="" false="--NO_BIOPERL" checkBioPerl} \
			~{true="--CONVERT" false="" convert} \
			~{default="" "--CACHEURL " + cacheURL} \
			~{default="" "--FASTAURL " + fastaURL}
	>>>

	output {
		File info = "~{cacheDirFinal}/~{species_complete}/~{cacheVersion}_~{assembly}/info.txt"
	}

 	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
 	}

 	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "vep_install"]',
			category: 'System'
		}
		destDir: {
			description: "Set destination directory for API install [default: './']",
			category: "Output path"
		}
		cacheVersion: {
			description: "Set data (cache, FASTA) version to install [default: 101]",
			category: "Set cache"
		}
		cacheDir: {
			description: "Set destination directory for cache files",
			category: "Output path"
		}
		pluginsDir: {
			description: "Set destination directory for VEP plugins files",
			category:  "Output path"
		}
		refseq: {
			description: "Defined if install refseq cache [default: true]",
			category: "Set cache"
		}
		ensembl: {
			description: "Defined if install ensembl cache [default: true]",
			category: "Set cache"
		}
		api: {
			description: "Defined if install api cache [default: true]",
			category: "Set cache"
		}
		cache: {
			description: "Defined if install cache [default: true]",
			category: "Set cache"
		}
		fasta: {
			description: "Defined if install fasta in cache [default: true]",
			category: "Set fasta"
		}
		plugins: {
			description: "Defined if install plugins in auto [default: false]",
			category: "Set plugins"
		}
		quiet: {
			description: "Don't write any status output when using --AUTO [default: false]",
			category: "Mode"
		}
		species: {
			description: "Defined species to install [default: 'homo_sapiens']",
			category: "Set cache"
		}
		assembly: {
			description: "Defined assembly to install [default: 'GRCh37']",
			category: "Set cache"
		}
		pluginsList: {
			description: "List plugins to install",
			category: "Set plugins"
		}
		prefer_bin: {
			description: "Active mode 'prefer_bin' [default: false]",
			category: "Mode"
		}
		checkHTSLIB: {
			description: "Deactive mode 'NO_HTSLIB' [default: false]",
			category: "Mode"
		}
		checkBioPerl: {
			description: "Deactive mode 'NO_BIOPERL' [default: false]",
			category: "Mode"
		}
		convert: {
			description: "Convert downloaded caches to use tabix for retrieving co-located variants [default: false]",
			category: "Mode"
		}
		cacheURL: {
			description: "Override default cache URL",
			category: "Set cache"
		}
		fastaURL: {
			description: "Override default fasta URL",
			category: "Set Fasta"
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
