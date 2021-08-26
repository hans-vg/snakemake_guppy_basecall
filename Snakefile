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
       expand('basecall/{sample}/sequencing_summary.txt', sample=SAMPLES),
#       "qc/multiqc.html"


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
        "guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params} --compress_fastq -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"


#rule guppy_linkfastq:
#    input:
#        #glob_wildcards("basecall/{sample}/*/*.fastq.gz"),
#        "basecall/{sample}/pass/*.fastq.gz",
#    output:
#        "basecall/{sample}.fastq.gz",
#    shell:
#        "ln -s {input} {output}"
#
#rule fastqc_pretrim:
#    input:
#        #"basecall/{sample}/{failpass}/{runid}.fastq.gz",
#        "basecall/{sample}.fastq.gz"
#    output:
#        html="qc/fastqc_pretrim/{sample}.html",
#        zip="qc/fastqc_pretrim/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#    params: ""
#    log:
#        "logs/fastqc_pretrim/{sample}.log"
#    threads: 1
#    wrapper:
#        "v0.75.0/bio/fastqc"
#
#rule multiqc:
#    input:
#        #expand("basecall/{sample}.fastq.gz", sample=SAMPLES)
#        expand("qc/fastqc_pretrim/{sample}_fastqc.zip", sample=SAMPLES)
#    output:
#        "qc/multiqc.html"
#    params:
#        ""  # Optional: extra parameters for multiqc.
#    log:
#        "logs/multiqc.log"
#    wrapper:
#        "0.77.0/bio/multiqc"

#rule fastqc_pretrim:
#    input:
#        "basecall/{sample}/{failpass}/{runid}.fastq.gz",
#    output:
#        html="qc/fastqc_pretrim/{sample}_{failpass}_{runid}.html",
#        zip="qc/fastqc_pretrim/{sample}_{failpass}_{runid}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#    params: ""
#    log:
#        "logs/fastqc_pretrim/{sample}_{failpass}_{runid}.log"
#    #resources: time_min=320, mem_mb=8000, cpus=1
#    threads: 1
#    wrapper:
#        "v0.75.0/bio/fastqc"
