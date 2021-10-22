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



task annovarForMpa {
  meta {
    author: "Thomas Guignard; Olivier Ardouin (modifications)"
    email: "t-guignard(at)chu-montpellier.fr;o-ardouin(at)chu-montpellier.fr"
    version : "0.0.2"
    date: "2021-10-22"
  }

  input {
    File CustomXref
    File vcf
    File RefAnnotateVariation
    File RefCodingChange
    File RefConvert2Annovar
    File RefRetrieveSeqFromFasta
    File RefVariantsReduction
    File TableAnnovarExe
    String HumanDb
    String outputName
    String OutDir = "/."
    String perlExe = "perl"
    #databases
    String Genome
    String Clinvar
    String Dbnsfp
    #String Spidex
    String Dbscsnv
    String GnomadExome
    String GnomadGenome
    String PopFreqMax
    String Intervar
    String SpliceAI
    ## run time
    Int threads = 1
		Int memoryByThreads = 768
		String? memory
  }

  String Dollar = "$"

  String outputName = if defined(name) then name else sub(basename(in),".vcf", "")
	String outputFile = if defined(outputPath) then outputPath + "/" + outputName else outputName

  String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
  Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
  Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
  Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
  Int memoryByThreadsMb = floor(totalMemMb/threads)

  command <<<
    OPERATION_SUFFIX=',f'
    COMMA=','
    POPFREQMAX=',${PopFreqMax}'
    #REFGENE='refGeneWithVer'
    if [ ${Genome} == 'hg38' ];then
    OPERATION_SUFFIX=''
      COMMA=''
      POPFREQMAX=''
      #REFGENE='refGene'
      fi
      "~{perlExe}" "${TableAnnovarExe}" \
        "~{vcf}" \
        "${HumanDb}" \
        -thread "~{threads}" \
        -buildver "${Genome}" \
        -out "~{outputFile}" \
        -remove \
        -intronhgvs 80 \
        -protocol refGeneWithVer,refGeneWithVer,"${Clinvar}","${Dbnsfp}","${Dbscsnv}","${GnomadExome}","${GnomadGenome}","${Intervar}",regsnpintron,"${SpliceAI}""${Dollar}{POPFREQMAX}" \
        -operation gx,g,f,f,f,f,f,f,f,f"${Dollar}{OPERATION_SUFFIX}" \
        -nastring . \
        -vcfinput \
        -otherinfo \
        -arg '-splicing 5','-hgvs',,,,,,,,"${Dollar}{COMMA}" \
        -xref "${CustomXref}"
 >>>

 output {
  File outAnnotationVcf = "~{outputFile}.${Genome}_multianno.vcf"
  File outAnnotationAvinput = "~{outputFile}.avinput"
  File outAnnotationTxt = "~{outputFile}.${Genome}_multianno.txt"
 }
 runtime {
   cpu: "~{threads}"
   requested_memory_mb_per_core: "${memoryByThreadsMb}"
 }
}
