# ONT-sc-Splice: Generation of splicing junction matrices 

<!-- TOC -->

- [ONT-sc-Splice: Generation of splicing junction matrices](#ont-sc-splice-generation-of-splicing-junction-matrices)
    - [Preparation](#preparation)
        - [Config File](#config-file)
    - [Running workflow](#runing-workflow)
    - [Adapting the resources](#adapting-the-resources)
    - [Output files](#output-files)

<!-- /TOC -->


This part of the pipeline will generate the junction and exon counts as well as annotate the splicing events recovered in every sample. 

## Preparation

In the folder were we want to output the results we need to generate a folder with the name of the sample. Inside this folder we will add the bam result of the sequencing result preprocessing. 

If using ONT data processed by SiCeLoRe 2.1, the file that we need will be found in a folder with the following structure:

```bash
/path/to/sicelore/output/SAMPLEID/04b.matrices/molecules.tags.GE.bam
```

If using PacBio, processed by our custom pipeline, the file that we need will be found in a folder with the following structure:

```bash
/path/to/pacbio/output/SAMPLEID/mapped_corrected/SAMPLEID_me_corrected_modified.bam
```

We can copy these bams or instead generate soft links (`ln -s ${original_file} ${new_file_name}`) that will allow us to have the following folder structure in our target output folder: 

```bash
/path/to/your/workdir/splicing_runs_folder
├── SAMPLEID_ONT
│   └─── 04b.matrices
│       ├── molecules.tags.GE.bam -> /path/to/sicelore/output/SAMPLEID_ONT/04b.matrices/molecules.tags.GE.bam
│       └── molecules.tags.GE.bam.bai -> /path/to/sicelore/output/SAMPLEID_ONT/04b.matrices/molecules.tags.GE.bam.bai
└── SAMPLEID_PACBIO
    └── 04b.matrices
        ├── molecules.tags.GE.bam -> /path/to/pacbio/output/SAMPLEID_PACBIO/mapped_corrected/SAMPLEID_PACBIO_me_corrected_modified.bam
        └── molecules.tags.GE.bam.bai -> /path/to/pacbio/output/SAMPLEID_PACBIO/mapped_corrected/SAMPLEID_PACBIO_me_corrected_modified.bam.bai
```

One level above the location of out output folder, we need to add our metadata files. The metadata files have to contain the following specifications: 

1. The header first two characters are spaces (no TAB), all the other column separations can be tabs. 
2. The barcodes are rows of our metadata and have the format: "BARCODELETTERS_1"
3. A column that contains the categories that we want to compare (i.e. WT vs MUT or CELLTYPE1 vs CELLTYPE2)
4. A column that contains an indicator string of the cells that will be used in a comparison. For instance if we are comparing MUT and WT only on HSC cells, this will be the "CELLTYPE" column, and our filter will be "HSC"


Make sure of adding this metadata tables and giving them the name in a format `SAMPLEID_metadata.txt`. 

The structure of the target folder will then be: 

```bash
/path/to/your/workdir
├── sample_metadata_ont
└── splicing_runs_folder
    ├── SAMPLEID_ONT
    │   └─── 04b.matrices
    │       ├── molecules.tags.GE.bam -> /path/to/sicelore/output/SAMPLEID_ONT/04b.matrices/molecules.tags.GE.bam
    │       └── molecules.tags.GE.bam.bai -> /path/to/sicelore/output/SAMPLEID_ONT/04b.matrices/molecules.tags.GE.bam.bai
    └── SAMPLEID_PACBIO
        └── 04b.matrices
            ├── molecules.tags.GE.bam -> /path/to/pacbio/output/SAMPLEID_PACBIO/mapped_corrected/SAMPLEID_PACBIO_me_corrected_modified.bam
            └── molecules.tags.GE.bam.bai -> /path/to/pacbio/output/SAMPLEID_PACBIO/mapped_corrected/SAMPLEID_PACBIO_me_corrected_modified.bam.bai
```


### Config File

This file contains the details about the samples that we have placed in our target folder according with the description above. This file has to be placed in the folder where you launch the pipeline 
You can process multiple samples at the time 

```yaml
WORKDIR:
  "/path/to/your/workdir"
RUN_FILES:
  "path/to/snakemake/pipeline/JunctionCalling.Pipeline.Scripts"
SICERLORE_VERSION:
  "2.1"
SAMPLES:
  - "SAMPLEID1"
  #....
  - "SAMPLEIDN"
JOB_NAME:
  "snakemake_part2"
OUTPUTDIR:
  "splicing_runs_folder"
METADATADIR:
  "/sample_metadata_ont/"
LIBRARY_TYPE:
  "ONT" # Run only similar sample types together
```

## Running workflow

To run this part of the workflow, you can use the `launch.sh` script. Make sure to load the corresponding conda environment. You can find a recipe for the ont-sc-splice environment [here](https://github.com/landau-lab/ONT-sc-splice/blob/main/ont_sc_splice.yml)

## Adapting the resources 

In the `ne1` cluster, the current Snakemake version is `8.28.0` therefore an alternative that has been adapted to administrate resources for this pipeline is contained in the `new_slurm_profile folder`. The corresponding `config.yalm` file declares resource requirements for the jobs, as well as SLURM specific configurations. If your sample requires less resources, that file is the ideal place to adapt. More on how to work with profiles in slurm can be found [here](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) and in the [official Snakemake documentation](https://snakemake.readthedocs.io/en/v8.28.0/snakefiles/configuration.html)

You can monitor the running of each process in the `logs` folder. 

Sometimes, if you need to re-run individual steps, you can follow a configuration similar to: 

```bash
snakemake --profile new_slurm_profile \
          --configfile config.yaml \
          --rerun-triggers mtime \
          --rerun-incomplete \
          --force /location/of/your/output/file
```

## Output files 

When successful, the pipeline will produce the following files per sample: 

```bash
/path/to/your/workdir/splicing_runs_folder/SAMPLEID/leafcutter_outputs_SAMPLEID
├── exon.meta
│   ├── SAMPLEID_all.exons.info.Rdata
│   ├── SAMPLEID_all.exons.info.txt
│   ├── SAMPLEID_all.exons.info.wRIannotation.csv
│   ├── SAMPLEID_counts_sc.rda
│   ├── SAMPLEID_junc.annotation.intermediates.Rdata
│   ├── SAMPLEID_junc.meta.rda
│   └── SAMPLEID_per.exon_numbers.counts.txt
├── intron.meta
│   ├── leafviz_all_introns.bed
│   ├── leafviz_all_introns_cleaned.bed
│   ├── SAMPLEID_all.introns.info.Rdata
│   ├── SAMPLEID_all.introns.info.txt
│   ├── SAMPLEID_all.introns.info.w.primary.annotations.exonSkipping.txt
│   ├── SAMPLEID_perind_numbers.counts.txt
│   └── SAMPLEID_strand_adjusted_metadata
│       └── strand_adjusted_metadata.csv
├── SAMPLEID_exons.tsv.gz
└── SAMPLEID_introns.tsv.gz
```

At the same time, in the working directory you expect to find a folder with all the files that we will need for step 3: 

```bash
/path/to/your/workdir/splicing_runs_folder/DT_input
├── 1.JuncCounts_matrix
│   ├── SAMPLEID1_perind_numbers.counts.txt
│   └── SAMPLEIDN_perind_numbers.counts.txt
├── 1.JuncCounts_matrix_exon
│   ├── SAMPLEID1_per.exon_numbers.counts.txt
│   └── SAMPLEIDn_per.exon_numbers.counts.txt
├── 2.Cell_metadata
│   ├── SAMPLEID1_metadata.txt
│   └── SAMPLEIDn_metadata.txt
├── 3.Junction_metadata
│   ├── SAMPLEID1_strand_adjusted_metadata.csv
│   └── SAMPLEIDn_strand_adjusted_metadata.csv
└── 3.Junction_metadata_exon
    ├── SAMPLEID1_all.exons.info.wRIannotation.csv
    └── SAMPLEIDN_all.exons.info.wRIannotation.csv
```

