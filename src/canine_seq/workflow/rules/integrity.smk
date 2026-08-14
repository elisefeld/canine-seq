rule integrity:
    input:
        ".../resources/data/raw/ORBIT/{sample}.cram" #eventually direct this to the raw_data_dir in config.yaml for more flexibility

    output:
        

    log:
        "logs/integrity/{sample}.log"
    conda:
     
    shell:
        samtools quickcheck {input} > {output}