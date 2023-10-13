
rule get_sketch_dist:
    """
    currently dashing v2.1.19
    emits similarity by default - if changing args to a different metric, make sure to change script process_sketch_dist.R accordingly
    """
    input:
        fa =  config.get("TE_FA"),
    output:
        sketch_dist = "results/sketch/te_sketch_dist.txt",
    params:
        args = config.get("DASHING_ARGS"),
        path = config.get("DASHING_PATH"),
    threads:
        10
    shell:
        """
        {params.path} dist {params.args} --parse-by-seq --square {input.fa} > {output.sketch_dist}
        """
        
