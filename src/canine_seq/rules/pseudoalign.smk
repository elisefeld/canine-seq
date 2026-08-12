rule pseudoalign:
    input:
        fastq="data/{sample}.fastq",
        idx=config["reference"]["fasta"] + ".idx"
    output:
        abundance="results/{sample}.txt"
    conda:
        "../envs/kallisto.yaml"
    threads: config["threads"]
    shell:
        "kallisto quant --index {input.idx} \
        --output-dir results/{wildcards.sample} \
        --threads {threads} {input.fastq} && mv results/{wildcards.sample}/abundance.tsv {output.abundance}"