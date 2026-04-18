shell.executable("/bin/bash")

""" Snakemake pipeline for 10X Genomics 3' single-cell RNA-seq and ONT long read joint isoform analysis, Part 2 """

configfile: "config.yaml"

WORKDIR = config["WORKDIR"]
RUN_FILES = config["RUN_FILES"]
SAMPLES = config["SAMPLES"]
JOB_NAME = config["JOB_NAME"]
OUTPUTDIR = config["OUTPUTDIR"]
METADATADIR = config["METADATADIR"]
LIBRARY_TYPE = config["LIBRARY_TYPE"]
LIBRARY_END = config["LIBRARY_END"]

rule all:
  input:
    #expand("{workdir}/{outputdir}/leafcutter_outputs/{samples}_exons.tsv.gz", workdir = WORKDIR, outputdir = OUTPUTDIR, samples = SAMPLES),
    expand("{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/exon.meta/{samples}_all.exons.info.wRIannotation.csv", workdir = WORKDIR, outputdir = OUTPUTDIR, samples = SAMPLES),
    expand("{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_strand_adjusted_metadata/strand_adjusted_metadata.csv", workdir = WORKDIR, outputdir = OUTPUTDIR, samples = SAMPLES),
    expand("{workdir}/{outputdir}/DT_input/1.JuncCounts_matrix/{samples}_perind_numbers.counts.txt", workdir = WORKDIR, outputdir = OUTPUTDIR, samples = SAMPLES),
    expand("{workdir}/{outputdir}/DT_input/3.Junction_metadata/{samples}_strand_adjusted_metadata.csv", workdir = WORKDIR, outputdir = OUTPUTDIR, samples = SAMPLES),
    expand("{workdir}/{outputdir}/DT_input/1.JuncCounts_matrix_exon/{samples}_per.exon_numbers.counts.txt", workdir = WORKDIR, outputdir = OUTPUTDIR, samples = SAMPLES),
    expand("{workdir}/{outputdir}/DT_input/3.Junction_metadata_exon/{samples}_all.exons.info.wRIannotation.csv", workdir = WORKDIR, outputdir = OUTPUTDIR, samples = SAMPLES)


#""" Generate Intron/3p/5p databases from GTF file (only need to do this once) """
#rule generate_db:

""" Clean intron.bed file """
rule clean_bed:
  input:
    "{workdir}/{outputdir}/{samples}/04b.matrices/molecules.tags.GE.bam"
    #"{workdir}/{outputdir}/{samples}/sicelore_outputs_{samples}/consensus_outputs/consensus.sorted.tags.GE.bam"
  output:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/leafviz_all_introns_cleaned.bed"
  params:
    job_name = "clean_bed",
    memory = "select[mem>100] rusage[mem=100]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    run_files = RUN_FILES,
    samples = "{samples}"
  log:
    "{workdir}/{outputdir}/logs/{samples}_clean_bed.txt"
  threads:
    10
  shell:
    """
    cd {params.workdir}
    cd {params.outputdir}
    mkdir -p {params.samples}/leafcutter_outputs_{params.samples}
    cd {params.samples}/leafcutter_outputs_{params.samples}
    mkdir -p intron.meta
    cp {params.run_files}/annotation_reference/leafviz_all_introns.bed.gz intron.meta/
    cd intron.meta
    gunzip -f leafviz_all_introns.bed.gz
    
    python {params.run_files}/intron_bed_cleaner.py {params.workdir}/{params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta/leafviz_all_introns.bed {params.workdir}/{params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta/leafviz_all_introns_cleaned.bed
  """

""" Run the python junction calling script """
rule junc_call:
  input:
    "{workdir}/{outputdir}/{samples}/04b.matrices/molecules.tags.GE.bam"
    #"{workdir}/{outputdir}/{samples}/sicelore_outputs_{samples}/consensus_outputs/consensus.sorted.tags.GE.bam"
  output:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/{samples}_exons.tsv.gz",
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/{samples}_introns.tsv.gz"
  params:
    job_name = "junc_call",
    memory = "select[mem>100] rusage[mem=100]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    run_files = RUN_FILES,
    samples = "{samples}",
    library_type = LIBRARY_TYPE
    library_end = LIBRARY_END

  log:
    "{workdir}/{outputdir}/logs/{samples}_junc_call.txt"
  threads:
    10
  shell:
    """
    cd {params.workdir}
    mkdir -p {params.outputdir}
    cd {params.outputdir}
    mkdir -p {params.samples}/leafcutter_outputs_{params.samples}
    cd {params.samples}/leafcutter_outputs_{params.samples}

    # Run the python junction calling script
    python {params.run_files}/count_introns_exons_ONT_and_PB.py {input} {params.samples} {params.library_type} {params.library_end}
    echo "junction calling done"
    """

""" Run the pre-processing scripts """
rule preproc:
  input:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/{samples}_exons.tsv.gz"
  output:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/exon.meta/{samples}_all.exons.info.txt",
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/exon.meta/{samples}_per.exon_numbers.counts.txt"
  params:
    job_name = "preproc",
    memory = "select[mem>100] rusage[mem=100]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    run_files = RUN_FILES,
    samples = "{samples}",
  log:
    "{workdir}/{outputdir}/logs/{samples}_preproc.txt"
  threads:
    10
  shell:
    """
    cd {params.workdir}
    cd {params.outputdir}
    cd {params.samples}/leafcutter_outputs_{params.samples}
    mkdir -p exon.meta
    
    Rscript {params.run_files}/exon.annotation.final.V4.R {params.samples}_exons.tsv.gz exon.meta/ {params.samples} {params.run_files}/annotation_reference/leafviz
    echo "Done"
    """

""" RI annotation script """
rule ri_anno:
  input:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/exon.meta/{samples}_all.exons.info.txt",
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/leafviz_all_introns_cleaned.bed"
  output:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/exon.meta/{samples}_all.exons.info.wRIannotation.csv"
  params:
    job_name = "ri_anno",
    memory = "select[mem>100] rusage[mem=100]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    run_files = RUN_FILES,
    samples = "{samples}",
  log:
    "{workdir}/{outputdir}/logs/{samples}_ri_anno.txt"
  threads:
    10
  shell:
    """
    cd {params.workdir}
    cd {params.outputdir}
    cd {params.samples}/leafcutter_outputs_{params.samples}

    
    Rscript {params.run_files}/retained.intron.annotation.V2.R exon.meta/{params.samples}_all.exons.info.txt {params.workdir}/{params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/exon.meta/ {params.samples} {params.workdir}/{params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta/leafviz_all_introns_cleaned.bed
    echo "Done RI annotation script"
    """

""" Running intron annotation script """
rule intron_anno:
  input:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/{samples}_introns.tsv.gz"
  output:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_perind_numbers.counts.txt",
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_all.introns.info.txt"
  params:
    job_name = "intron_anno",
    #memory = "select[mem>100] rusage[mem=100]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    run_files = RUN_FILES,
    samples = "{samples}"
  log:
    "{workdir}/{outputdir}/logs/{samples}_intron_anno.txt"
  threads:
    10
  shell:
    """
    cd {params.workdir}
    cd {params.outputdir}
    cd {params.samples}/leafcutter_outputs_{params.samples}
    mkdir -p intron.meta
        
    Rscript {params.run_files}/junc_calling_script.R {params.samples}_introns.tsv.gz intron.meta/ {params.samples} {params.run_files}/annotation_reference/leafviz
    echo "Done intron annotation"
    """

""" Exon-skipping annotation script """
rule es_anno:
  input:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_all.introns.info.txt",
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/leafviz_all_introns_cleaned.bed"
  output:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_all.introns.info.w.primary.annotations.exonSkipping.txt"
  params:
    job_name = "es_anno",
    memory = "select[mem>100] rusage[mem=100]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    run_files = RUN_FILES,
    samples = "{samples}"
  log:
    "{workdir}/{outputdir}/logs/{samples}_es_anno.txt"
  threads:
    10
  shell:
    """
    cd {params.workdir}
    cd {params.outputdir}
    cd {params.samples}/leafcutter_outputs_{params.samples}
    #cp {params.run_files}/annotation_reference/leafviz_all_introns.bed.gz intron.meta/
    cd intron.meta
    #gunzip -f leafviz_all_introns.bed.gz

    #clean the intron.bed file for exon skipping
    #python {params.run_files}/intron_bed_cleaner.py {params.workdir}/{params.outputdir}/leafcutter_outputs/intron.meta/leafviz_all_introns.bed {params.workdir}/{params.outputdir}/leafcutter_outputs/intron.meta/leafviz_all_introns_cleaned.bed
    python {params.run_files}/new_annotator_with_skipping.py {params.workdir}/{params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta {params.samples}_all.introns.info.txt {params.workdir}/{params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta/leafviz_all_introns.bed {params.workdir}/{params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta/leafviz_all_introns_cleaned.bed ONT {params.samples}_all.introns.info.w.primary.annotations.exonSkipping.txt
    """

""" Run Strand Adjustment of Metadata """
rule strand_adjust:
  input:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_all.introns.info.w.primary.annotations.exonSkipping.txt"
  output:
    "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_strand_adjusted_metadata/strand_adjusted_metadata.csv"
  params:
    job_name = "strand_adjust",
    memory = "select[mem>100] rusage[mem=100]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    run_files = RUN_FILES,
    samples = "{samples}"
  log:
    "{workdir}/{outputdir}/logs/{samples}_strand_adjust.txt"
  #threads:
  #  10
  shell:
    """
    cd {params.workdir}
    mkdir -p {params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta/{params.samples}_strand_adjusted_metadata
    
    cd {params.outputdir}/{params.samples}/leafcutter_outputs_{params.samples}/intron.meta

    Rscript {params.run_files}/strand_adjustment.R {params.samples}_all.introns.info.w.primary.annotations.exonSkipping.txt {params.samples}_strand_adjusted_metadata
    """

""" Copy files as input for next step """
rule copy_files:
  input:
    junc_mat = "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_perind_numbers.counts.txt",
    cell_meta = "{workdir}" + config["METADATADIR"] + "/{samples}_metadata.txt",
    junc_meta = "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/intron.meta/{samples}_strand_adjusted_metadata/strand_adjusted_metadata.csv",
    junc_mat_exon = "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/exon.meta/{samples}_per.exon_numbers.counts.txt",
    junc_meta_exon = "{workdir}/{outputdir}/{samples}/leafcutter_outputs_{samples}/exon.meta/{samples}_all.exons.info.wRIannotation.csv"
  output:
    junc_mat = "{workdir}/{outputdir}/DT_input/1.JuncCounts_matrix/{samples}_perind_numbers.counts.txt",
    cell_meta = "{workdir}/{outputdir}/DT_input/2.Cell_metadata/{samples}_metadata.txt",
    junc_meta = "{workdir}/{outputdir}/DT_input/3.Junction_metadata/{samples}_strand_adjusted_metadata.csv",
    junc_mat_exon = "{workdir}/{outputdir}/DT_input/1.JuncCounts_matrix_exon/{samples}_per.exon_numbers.counts.txt",
    junc_meta_exon = "{workdir}/{outputdir}/DT_input/3.Junction_metadata_exon/{samples}_all.exons.info.wRIannotation.csv"
  params:
    job_name = "copy_files",
    memory = "select[mem>10] rusage[mem=10]",
    workdir = "{workdir}",
    outputdir = "{outputdir}",
    samples = "{samples}"
  log:
    "{workdir}/{outputdir}/logs/{samples}_copy_files.txt"
  threads:
    5
  shell:
    """
    cd {params.workdir}
    cd {params.outputdir}
    mkdir -p DT_input
    cd DT_input
    mkdir -p 1.JuncCounts_matrix
    mkdir -p 2.Cell_metadata
    mkdir -p 3.Junction_metadata
    mkdir -p 1.JuncCounts_matrix_exon
    mkdir -p 3.Junction_metadata_exon
    cp {input.junc_mat} {output.junc_mat}
    cp {input.cell_meta} {output.cell_meta}
    cp {input.junc_meta} {output.junc_meta}
    cp {input.junc_mat_exon} {output.junc_mat_exon}
    cp {input.junc_meta_exon} {output.junc_meta_exon}
    """
