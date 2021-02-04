## Configuration file
import os
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['annotation', 'organism', 'build', 'release', 'txome', 'genome', 'gtf', 'chrnames', 'salmonindex', 'STARindex', 'readlength', 'fldMean', 'fldSD', 'metatxt', 'design', 'contrast', 'genesets', 'ncores', 'FASTQ', 'fqsuffix', 'output', 'useCondaR', 'Rbin', 'run_trimming', 'run_STAR', 'run_camera']
for k in configvars:
	if k not in config:
		config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
	if str is None:
		str = ''
	return str

config['txome'] = sanitizefile(config['txome'])
config['gtf'] = sanitizefile(config['gtf'])
config['genome'] = sanitizefile(config['genome'])
config['STARindex'] = sanitizefile(config['STARindex'])
config['salmonindex'] = sanitizefile(config['salmonindex'])
config['chrnames'] = sanitizefile(config['chrnames'])
config['metatxt'] = sanitizefile(config['metatxt'])

## Read metadata
if not os.path.isfile(config["metatxt"]):
  sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

import pandas as pd
samples = pd.read_csv(config["metatxt"], sep='\t')

if not set(['names','type']).issubset(samples.columns):
  sys.exit("Make sure 'names' and 'type' are columns in " + config["metatxt"])


## Sanitize provided input and output directories
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

outputdir = getpath(config["output"])
FASTQdir = getpath(config["FASTQ"])

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

## Define the R binary
Rbin = config["Rbin"]

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		outputdir + "MultiQC/multiqc_report.html",
		outputdir + "outputR/edgeR_dge.rds",
		outputdir + "outputR/qc_summary.rds"

rule setup:
	input:
		outputdir + "Rout/pkginstall_state.txt",
		outputdir + "Rout/softwareversions.done"

## Install R packages
rule pkginstall:
	input:
		script = "scripts/install_pkgs.R"
	output:
		outputdir + "Rout/pkginstall_state.txt"
	params:
		flag = config["annotation"],
		ncores = config["ncores"],
		organism = config["organism"]
	priority:
		50
	conda:
		Renv
	log:
		outputdir + "Rout/install_pkgs.Rout"
	benchmark:
		outputdir + "benchmarks/install_pkgs.txt"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args outtxt='{output}' ncores='{params.ncores}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Salmon quantification
rule runsalmonquant:
	input:
		expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.names.values.tolist())

## STAR alignment
rule runstar:
	input:
		expand(outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist()),
		expand(outputdir + "STARbigwig/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.names.values.tolist())

## List all the packages that were used by the R analyses
rule listpackages:
	log:
		outputdir + "Rout/list_packages.Rout"
	params:
		Routdir = outputdir + "Rout",
		outtxt = outputdir + "R_package_versions.txt",
		script = "scripts/list_packages.R"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {params.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	output:
		touch(outputdir + "Rout/softwareversions.done")
	conda:
		"envs/environment.yaml"
	shell:
		"salmon --version; trim_galore --version; "
		"echo -n 'cutadapt ' && cutadapt --version; "
		"fastqc --version; STAR --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##
## Extract chromosome names
rule chrnames:
	input:
		genome = config["genome"]
	output:
		config["chrnames"]
	log:
		outputdir + "logs/get_chrnames.log"
	benchmark:
		outputdir + "benchmarks/get_chrnames.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"""
		grep ">" {input.genome} | cut -d ">" -f 2 | cut -d " " -f 1 > {output}
		"""

## Generate Salmon index from transcriptome
rule salmonindex:
	input:
		txome = config["txome"],
		genome = config["genome"],
		chrnames = config["chrnames"]
	output:
		config["salmonindex"] + "/versionInfo.json"
	log:
		outputdir + "logs/salmon_index.log"
	benchmark:
		outputdir + "benchmarks/salmon_index.txt"
	params:
		salmonoutdir = config["salmonindex"],
		anno = config["annotation"],
		salmonextraparams = config["additional_salmon_index"]
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"""
		echo 'Salmon version:\n' > {log}; salmon --version >> {log};
		salmon index -t <(cat {input.txome} {input.genome}) -i {params.salmonoutdir} \
		--gencode {params.salmonextraparams} -d {input.chrnames} -p {threads}
		"""

## Generate STAR index
rule starindex:
	input:
		genome = config["genome"],
		gtf = config["gtf"]
	output:
		config["STARindex"] + "/SA",
		config["STARindex"] + "/chrNameLength.txt"
	log:
		outputdir + "logs/STAR_index.log"
	benchmark:
		outputdir + "benchmarks/STAR_index.txt"
	params:
		STARindex = config["STARindex"],
		readlength = config["readlength"],
		starextraparams = config["additional_star_index"]
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.STARindex} "
		"--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength} "
		"{params.starextraparams}"

## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_{sample}.log"
	benchmark:
		outputdir + "benchmarks/fastqc_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## FastQC, trimmed reads
rule fastqctrimmed:
	input:
		fastq = outputdir + "FASTQtrimmed/{sample}.fq.gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_trimmed_{sample}.log"
	benchmark:
		outputdir + "benchmarks/fastqc_trimmed_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"


# The config.yaml files determines which steps should be performed
def multiqc_input(wildcards):
	input = []
	input.extend(expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()))
	input.extend(expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.names.values.tolist()))
	if config["run_trimming"]:
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()))
	if config["run_STAR"]:
		input.extend(expand(outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist()))
	return input

## Determine the input directories for MultiQC depending on the config file
def multiqc_params(wildcards):
	param = [outputdir + "FastQC",
	outputdir + "salmon"]
	if config["run_trimming"]:
		param.append(outputdir + "FASTQtrimmed")
	if config["run_STAR"]:
		param.append(outputdir + "STAR")
	return param

## MultiQC
rule multiqc:
	input:
		multiqc_input
	output:
		outputdir + "MultiQC/multiqc_report.html"
	params:
		inputdirs = multiqc_params,
		MultiQCdir = outputdir + "MultiQC"
	log:
		outputdir + "logs/multiqc.log"
	benchmark:
		outputdir + "benchmarks/multiqc.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"

# Summary QC plots
rule qcsummary:
	input:
		outputdir + "MultiQC/multiqc_report.html",
		script = "scripts/summarize_qc.R",
		metatxt = config["metatxt"]
	output:
		outputdir + "outputR/qc_summary.rds"
	log:
		outputdir + "Rout/summarize_qc.Rout"
	benchmark:
		outputdir + "benchmarks/summarize_qc.txt"
	params:
		multiqcdir = outputdir + "MultiQC",
		fastqcdir = outputdir + "FastQC"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args metafile='{input.metatxt}' multiqcdir='{params.multiqcdir}' fastqcdir='{params.fastqcdir}' outrds='{output}'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	params:
		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
	log:
		outputdir + "logs/trimgalore_{sample}.log"
	benchmark:
		outputdir + "benchmarks/trimgalore_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"

## ------------------------------------------------------------------------------------ ##
## Salmon abundance estimation
## ------------------------------------------------------------------------------------ ##
# Estimate abundances with Salmon
rule salmonSE:
	input:
		index = config["salmonindex"] + "/versionInfo.json",
		fastq = outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "salmon/{sample}/quant.sf"
	log:
		outputdir + "logs/salmon_{sample}.log"
	benchmark:
		outputdir + "benchmarks/salmon_{sample}.txt"
	threads:
		config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		salmondir = outputdir + "salmon",
		salmonextraparams = config["additional_salmon_quant"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -r {input.fastq} "
		"-o {params.salmondir}/{wildcards.sample} -p {threads} {params.salmonextraparams}"

## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
rule starSE:
	input:
		index = config["STARindex"] + "/SA",
		fastq = outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/STAR_{sample}.log"
	benchmark:
		outputdir + "benchmarks/STAR_{sample}.txt"
	params:
		STARindex = config["STARindex"],
		STARdir = outputdir + "STAR",
		starextraparams = config["additional_star_align"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"

## Index bam files
rule bamindex:
	input:
		bam = outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_index_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
		chrl = config["STARindex"] + "/chrNameLength.txt"
	output:
		outputdir + "STARbigwig/{sample}_Aligned.sortedByCoord.out.bw"
	params:
		STARbigwigdir = outputdir + "STARbigwig"
	log:
		outputdir + "logs/bigwig_{sample}.log"
	benchmark:
		outputdir + "benchmarks/bigwig_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | LC_COLLATE=C sort -k1,1 -k2,2n > "
		"{params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## ------------------------------------------------------------------------------------ ##
## Transcript quantification
## ------------------------------------------------------------------------------------ ##
## tximeta
rule tximeta:
	input:
	    outputdir + "Rout/pkginstall_state.txt",
		expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.names.values.tolist()),
		metatxt = config["metatxt"],
		salmonidx = config["salmonindex"] + "/versionInfo.json",
		script = "scripts/run_tximeta.R"
	output:
		outputdir + "outputR/tximeta_se.rds"
	log:
		outputdir + "Rout/tximeta_se.Rout"
	benchmark:
		outputdir + "benchmarks/tximeta_se.txt"
	params:
		salmondir = outputdir + "salmon",
		flag = config["annotation"],
		organism = config["organism"]
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args salmondir='{params.salmondir}' metafile='{input.metatxt}' outrds='{output}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Input variable check
## ------------------------------------------------------------------------------------ ##
def geneset_param(wildcards):
	if config["run_camera"]:
                gs = config["genesets"].replace(" ", "") if config["genesets"] is not None else "NOTDEFINED"
		return "genesets='" + gs + "'"
	else:
		return ""


## check design matrix and contrasts
rule checkinputs:
    input:
        "config.yaml",
        script = "scripts/check_input.R"
    output:
        outputdir + "Rout/check_input.txt"
    log:
        outputdir + "Rout/check_input.Rout"
    benchmark:
    	outputdir + "benchmarks/check_input.txt"
    params:
        gtf = config["gtf"],
        genome = config["genome"],
        txome = config["txome"],
        fastqdir = config["FASTQ"],
        metatxt = config["metatxt"],
        design = config["design"].replace(" ", "") if config["design"] is not None else "NOTDEFINED",
        contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "NOTDEFINED",
        annotation = config["annotation"].replace(" ", "") if config["annotation"] is not None else "NOTDEFINED",
        genesets = geneset_param,
        fqsuffix = str(config["fqsuffix"]),
        run_camera = str(config["run_camera"]),
        organism = config["organism"]    
    conda:
	    Renv
    shell:
        '''{Rbin} CMD BATCH --no-restore --no-save "--args metafile='{params.metatxt}' design='{params.design}' contrast='{params.contrast}' outFile='{output}' gtf='{params.gtf}' genome='{params.genome}' fastqdir='{params.fastqdir}' fqsuffix='{params.fqsuffix}' txome='{params.txome}' run_camera='{params.run_camera}' organism='{params.organism}' {params.genesets} annotation='{params.annotation}'" {input.script} {log};
        cat {output}
        '''
       

## ------------------------------------------------------------------------------------ ##
## Differential expression
## ------------------------------------------------------------------------------------ ##
rule edgeR:
	input:
		outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "outputR/tximeta_se.rds",
		script = "scripts/run_render.R",
		template = "scripts/edgeR_dge.Rmd"
	output:
		html = outputdir + "outputR/edgeR_dge.html",
		rds = outputdir + "outputR/edgeR_dge.rds"
	params:
		directory = outputdir + "outputR",
		organism = config["organism"],        
		design = config["design"].replace(" ", "") if config["design"] is not None else "",
		contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "",
		genesets = geneset_param
	log:
		outputdir + "Rout/run_dge_edgeR.Rout"
	benchmark:
		outputdir + "benchmarks/run_dge_edgeR.txt"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' organism='{params.organism}' design='{params.design}' contrast='{params.contrast}' {params.genesets} rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='edgeR_dge.html'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
