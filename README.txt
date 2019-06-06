Necessary files:
* 'lowef_med_meta_bn.txt'
* 'lowef_med_heals_meta_bn.txt'
* 'lowef_med_best_meta_bn.txt'
* 'lowef_med_heals_epic_bn.txt'
  ** These files already have binary and weighted scores coded

-------------
Ensure that you are in the MR_arsenic directory when executing
each R script in the pipeline.
-------------

Pipeline for replication of analysis:
1) Run arsenic_corr.R to generate Figure S1
2) Run descriptive_stats.R to generate source Excel file ('describe.csv')
 for Table 2
3) Run general_regression.R to generate results of regression analyses.
4) Run MR.R to generate results of MR analyses.
5) Run consistency_test.R to make the consistency files and pie charts.
6) Run MR_plots.R to make MR plots and forest plots for top CpGs.
7) Run forest_plots.R to make forest plots for all analyses.