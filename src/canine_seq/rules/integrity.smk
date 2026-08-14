rule integrity:
    input:
        ".../resources/data/raw/ORBIT/{sample}.cram" #eventually direct this to the raw_data_dir in config.yaml for more flexibility

    output:
        
    conda:
        "../envs/kallisto.yaml"
    shell:
        samtools quickcheck {input} > {output}