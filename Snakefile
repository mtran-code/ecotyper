configfile: "config/config.yml"

include: "workflow/Snakefile"


rule all:
    default_target: True
    input:
        ALL_TARGETS
