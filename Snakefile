configfile: "config.yaml"


rule all:
    input:
        "report.html"


rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    threads: 4
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    log:
        "logs/samtools/{sample}.sort.log"
    shell:
        "(samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}) 2> {log}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    log:
        "logs/samtools/{sample}.index.log"
    shell:
        "(samtools index {input}) 2> {log}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    params:
        pmr="0.001"
    log:
        "logs/bcftools/all.vcf.log"
    shell:
        "(samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv -P '{params.pmr}' - > {output}) 2> {log}"


rule report:
    input:
        T1="calls/all.vcf",
        T2=expand("benchmarks/{sample}.bwa.benchmark.txt", sample=config["samples"])
    output:
        "report.html"
    script:
        "scripts/report.py"
