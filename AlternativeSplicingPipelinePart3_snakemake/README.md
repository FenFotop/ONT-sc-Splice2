# ONT-sc-Splice: Differential Splicing Analysis
<!-- TOC -->
<!-- TOC -->

- [ONT-sc-Splice: Differential Splicing Analysis](#ont-sc-splice-differential-splicing-analysis)
    - [Preparation](#preparation)
        - [Config File](#config-file)
    - [Running workflow](#running-workflow)
    - [Adapting the resources](#adapting-the-resources)
    - [Output files](#output-files)

<!-- /TOC -->
This part of the pipeline will perform the local differential splicing analysis. To do this, we perform a permutation test, comparing the changes in PSI (dPSI) between one condition and another against a null distribution.
The analysis is based in the conceptual framework of [LeafCutter](https://github.com/davidaknowles/leafcutter).

Given that the analysis relies in junctions and exon-junctions (intron retention analysis), our algorithm can produce information of the following local splicing events: 

- Alternative 3' splice sites
- Alternative 5' splice sites
- Exon skipping * 
- Intron retention *

\* *Limited* 

Thus additional events such as mixed exons or alternative 5' TSS and 3' alternative polyadenylation sites are limited. 

## Preparation

We start from the `DT_input` folder produced by the Part2 of this pipeline. 

### Config File

This file contains the details about the comparisons that we will perform. This file has to be placed in the folder where you launch the pipeline. 
You can process all samples or only a subset by specifying the samples parameter. 

```yaml
WORKDIR:
  "/path/to/your/workdir/"
RUN_FILES:
  "path/to/snakemake_folder/MultiSample.MultiGroupComp.AdjVars.Scripts"
RUN_FILESIR:
  "path/to/snakemake_folder/MultiSample.MultiGroupComp.AdjVars.Scripts.IR"
INTRONBED_FILE:
  "path/to/snakemake_folder/JunctionCalling.Pipeline.Scripts/annotation_reference/leafviz_all_introns_cleaned.bed"
SAMPLES:
  - "SAMPLEID"
NPERM:
  100000 # Can be adapted to the limit of recovery - Look at the p.value distribution
JOB_NAME:
  "snakemake_part3"
COMP_GROUPS_COLUMN:
  "GENOTYPECOLUMN" # Replace with the column name where the comparison groups are defined
COMP1:
  "WT"
COMP2:
  "MUT"
CELL_GROUPS_COLUMN:
  "CELLTYPE" # Replace with the column name that allows to filter the cells that will be used for the analysis
GROUPLIST:
  "ALL"
OUTPUTDIR:
  "splicing_runs_folder"
JUNC_COV_THRESH:
  5
N_CHUNK:
  5
```

## Running workflow

To run this part of the workflow, you can use the `launch.sh` script. Make sure to load the correspoding conda environment. You can find a recipe for the ont-sc-splice environment [here](https://github.com/landau-lab/ONT-sc-splice/blob/main/ont_sc_splice.yml)

## Adapting the resources 

In the `ne1` cluster, the current Snakemake version is `8.28.0` therefore an alternative that has been adapted to administrate resources for this pipeline is contained in the `new_slurm_profile folder`. The corresponding `config.yalm` file declares resource requiriments for the jobs, as well as SLURM specific configurations. If your sample requires less resources, that file is the ideal place to adapt. More on how to work with profiles in SLURM can be found [here](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) and in the [official Snakemake documentation](https://snakemake.readthedocs.io/en/v8.28.0/snakefiles/configuration.html)

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

When successful, the pipeline will produce the following files per sample in a folder with the same name as the parameter used as filter (cell type or any other).  

```bash
/path/to/your/workdir/splicing_runs_folder/ALL
├── logs
├── merge_annotate.log
├── merge_final_output
│   ├── logOR_within_cell_type_ALT_3P_Junctions.txt
│   ├── logOR_within_cell_type_ALT_3P_Junctions_with_threshold_info.txt
│   ├── logOR_within_cell_type_ALT_5P_Junctions.txt
│   ├── logOR_within_cell_type_ALT_5P_Junctions_with_threshold_info.txt
│   ├── logOR_within_cell_type_ONLY_CRYPTIC_3P_Junctions_with_threshold_info.txt
│   ├── logOR_within_cell_type_ONLY_CRYPTIC_5P_Junctions_with_threshold_info.txt
│   └── splicing_threshold_distribution.pdf
├── merge_final_output_IR
│   ├── logOR_within_cell_type_ALT_3P_Junctions.txt
│   └── logOR_within_cell_type_IR_Junctions_with_threshold_info.txt
├── split_cluster_files
|   ...
├── split_cluster_files_IR
|   ...
├── split_cluster_output
|   ...
└── split_cluster_output_IR
```

The most important files for downstream analysis are all the ones that end with `Junctions_with_threshold_info.txt`. They contain the following columns: 

1. **pvalue** - Calculated from the permutation test
2. **three.obs.logOR.ratio** - Log odds ratio 
3. **intron_junction** - Location of the junction
4. **chr** - Chromosome 
5. **start** - Junction start
6. **end** - Junction end 
7. **strand_1** - strand 
8. **five_prime** - Location of the closest 5' annotated in reference 
9. **three_prime** - Location of the closest 3' annotated in reference 
10. **five_prime_ID** - ID of the 5' junction cluster
11. **three_prime_ID** - ID of the 3' junction cluster
12. **unique_ID** - Unique ID for the cluster
13. **obs.group1** - Total counts for group 1 
14. **obs.group2** - Total counts for group 2
15. **total.reads.per.junction** - Total reads per junction (group 1 + group 2)
16. **startClass** - Classification of the splice site based on the start site 
17. **endClass** - Classification of the splice site based on the end site 
18. **alt_(three/five)_prime_intron_junction** - Class of junction (alternative or main, according to reference)
19. **(tree/five).group1.cluster.cov** - Total counts for group 1 cluster
20. **(tree/five).group2.cluster.cov** - Total counts for group 2 cluster
21. **strand** - Strand
22. **gene** - Gene symbol 
23. **verdict** - Classification level 1 
24. **fivep_class** - Classification of the junction based on the 5' end site
25. **threep_class** - Classification of the junction based on the 3' end site 
26. **fivep_distance** - Distance (in nts) to the closest 5' annotated in reference 
27. **threep_distance** -  Distance (in nts) to the closest 3' annotated in reference 
28. **relStartExon_skipping** - Start of the closest exon to account for skipping 
29. **relEndExon_skipping** - End of the closest exon to account for skipping 
30. **group2.psi** - Percent Splice In (PSI) in group 1
31. **group1.psi** - Percent Splice In (PSI) in group 2
32. **dPSI** - difference in PSI group 2 vs group 1 
33. **Final_Verdict** - Aggregated classification from the two splice sites 
34. **num.skipped.exons** - Total skipped exons
35. **exon.skip** - Classification based on the presence of skipped exons
36. **Final_Verdict.w.skipping** - Aggregated classification including exon skipping
37. **dPSI_threshold_0** - Is the absolute dPSI > than 0? 
38. **dPSI_threshold_5** - Is the absolute dPSI > than 5? 
39. **pvalue_threshold** - Is the pvalue < 0.05?
40. **dPSI_0_pvalue_threshold** - Is the pvalue < 0.05 & absolute dPSI > than 0 ?
41. **dPSI_5_pvalue_threshold** - Is the pvalue < 0.05 & absolute dPSI > than 5 ?

Ideally a change of PSI > than 5 represents a change in splicing, between two conditions, of at least 5%. The p. value can be adjusted. 

Note: Multiple junctions can be in the same cluster, some of them can be significant but with opposite dPSI values. Sometimes is good to check for duplicates that can come from these types of events.  

An example of possible filters and downstream visualizations can be found in the [Splicing_report_updated.Rmd](./Splicing_report_updated.Rmd) script. 