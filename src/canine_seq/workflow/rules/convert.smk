rule convert:
    input:
        "{sample}.cram",
    output:
        bam="{sample}.bam",
        idx="{sample}.bai",
    log:
        "{sample}.log",
    params:
        extra="",  # optional params string
        region="",  # optional region string
    threads: 2
    wrapper:
        "v9.15.0/bio/samtools/view"