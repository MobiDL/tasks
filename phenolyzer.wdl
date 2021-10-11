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


task phenolyzer {
	Boolean IsPrepared
	String PhenolyzerExe
	String DiseaseFile
	String WorkflowType
	String SampleID
	String OutDir
	String PerlPath
	#runtime attributes
	Int Cpu
	Int Memory
	command <<<
		cd ${PhenolyzerExe}
		${PerlPath} disease_annotation.pl ${DiseaseFile} -f -p -ph -logistic -out ../..${OutDir}${SampleID}/${WorkflowType}/disease/${SampleID}
	>>>
	output {
		String? outPhenolyzer = "${OutDir}${SampleID}/${WorkflowType}/disease/${SampleID}.predicted_gene_scores"
	}
	runtime {
		cpu: "${Cpu}"
		requested_memory_mb_per_core: "${Memory}"
	}
}
