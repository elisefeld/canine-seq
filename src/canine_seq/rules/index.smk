rule index:
    input:
        fasta=config["reference"]["fasta"]
    output:
        idx=config["reference"]["fasta"] + ".idx"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index --index {output.idx} {input.fasta}"