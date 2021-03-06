# GI_ReportingTool
Genetic Interaction Reporting Tool to collect information of mutually exclusive and co-occurring altered gene pairs in childhood cancers. For more details about the candidates we refer to the accompanying manuscript **_A comprehensive map of genetic interactions in childhood cancer reveals multiple underlying biological mechanisms_** (https://doi.org/10.1101/2020.11.17.385120)


## Running the app
To aid in reproducibility and ease of deployment a _Dockerfile_ is provided. This file is not needed when starting the Shiny server from R. To do that, the _app.R_ file can be used. The _setup.R_ and _functions.R_ will be read in via the main _app.R_ file.  


## Data information

**_www/GI_map.png_**

Simplified version of Figure 3 of the manuscript.

**_data/cand_target_dkfz.txt_**

Same table as table S2 of the manuscript.

**_data/muts_cand_target_dkfz.txt_**

This file contains all mutations from the TARGET and DKFZ data set, but only in genes that were part of a candidate genetic interaction in one of the two data sets. Only samples that were included in the test are considered, so hypermutators etc. are filtered out.

**_data/vep_out_muts_target_dkfz.txt_**

This file is the output file of VEP after submitting a filtered version of _muts_cand_target_dkfz.txt_. The VEP input file contains only columns needed for VEP and the rows where _candidate=FALSE_ and _vartype=SYN_ were filtered out.


## Subset of the data
Unfortunately, we cannot share the complete tables, because of sensitive information. However, we will share the HGG-K27M, HGG-other, and MB-SHH information of DKFZ to get familiar with the tool. Please note, the files in the data folder only contains the information of these three cancer types and therefore only the associated candidates are shown in the reporting tool when running it from here. All candidates can be analyzed via http://gi-analysis.bioinf.prinsesmaximacentrum.nl/. 
