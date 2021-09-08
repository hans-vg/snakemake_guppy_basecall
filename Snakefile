import glob

configfile: "config.yaml"

inputdirectory=config["directory"]
SAMPLES, = glob_wildcards(inputdirectory+"/{sample}.fast5", followlinks=True)
print(SAMPLES)

wildcard_constraints:
    sample="\w+\d+_\w+_\w+\d+_.+_\d"


##### target rules #####
rule all:
    input: 
       expand("basecall/{sample}/sequencing_summary.txt", sample=SAMPLES),
       "qc/multiqc.html"


rule make_indvidual_samplefiles:
    input:
        inputdirectory+"/{sample}.fast5",
    output:
        "lists/{sample}.txt",
    shell:
        "basename {input}  > {output}"


rule guppy_basecall_persample:
    input:
        directory=directory(inputdirectory),
        samplelist="lists/{sample}.txt",
    output:
        summary="basecall/{sample}/sequencing_summary.txt",
        directory=directory("basecall/{sample}/"),
    params: 
        config["basealgo"]
    shell:
        "guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params} --trim_barcodes --compress_fastq -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"

#def aggregate_input(wildcards):
#    checkpoint_output = checkpoints.guppy_basecall_persample.get(**wildcards).output[1]
#    print(checkpoint_output)
#    exparr = expand("basecall/{sample}/pass/{runid}.fastq.gz", sample=wildcards.sample, runid=glob_wildcards(os.path.join(checkpoint_output, "pass/", "{runid}.fastq.gz")).runid) 
#    print(exparr)
#    return exparr
#
##SAMPLES, RUNIDS, = glob_wildcards("basecall/{sample}/pass/{runid}.fastq.gz", followlinks=True)
##print(RUNIDS)
##print(SAMPLES)
#
#
#rule fastqc_pretrim:
#    input:
#        aggregate_input
#    output:
#        html="qc/fastqc_pretrim/{sample}.html",
#        zip="qc/fastqc_pretrim/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#    params: ""
#    log:
#        "logs/fastqc_pretrim/{sample}.log"
#    threads: 1
#    wrapper:
#        "0.77.0/bio/fastqc"
#
#rule multiqc:
#    input:
#        #expand("basecall/{sample}.fastq.gz", sample=SAMPLES)
#        #expand("qc/fastqc_pretrim/{sample}_fastqc.zip", sample=SAMPLES)
#        expand(rules.fastqc_pretrim.output.zip, sample=SAMPLES)
#    output:
#        "qc/multiqc.html"
#    params:
#        ""  # Optional: extra parameters for multiqc.
#    log:
#        "logs/multiqc.log"
#    wrapper:
#        "0.77.0/bio/multiqc"

##rule fastqc_pretrim:
##    input:
##        "basecall/{sample}/{failpass}/{runid}.fastq.gz",
##    output:
##        html="qc/fastqc_pretrim/{sample}_{failpass}_{runid}.html",
##        zip="qc/fastqc_pretrim/{sample}_{failpass}_{runid}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
##    params: ""
##    log:
##        "logs/fastqc_pretrim/{sample}_{failpass}_{runid}.log"
##    #resources: time_min=320, mem_mb=8000, cpus=1
##    threads: 1
##    wrapper:
##        "v0.75.0/bio/fastqc"
