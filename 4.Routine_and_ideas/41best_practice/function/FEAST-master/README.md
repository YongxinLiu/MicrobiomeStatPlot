FEAST - a scalable algorithm for quantifying the origins of complex microbial communities
-----------------------

A major challenge of analyzing the compositional structure of microbiome data is identifying its potential origins. Here, we introduce Fast Expectation-mAximization microbial Source Tracking (FEAST), a ready-to-use scalable framework that can simultaneously estimate the contribution of thousands of potential source environments in a timely manner, thereby helping unravel the origins of complex microbial communities. The information gained from FEAST may provide insight into quantifying contamination, tracking the formation of developing microbial communities, as well as distinguishing and characterizing bacteria-related health conditions. 
For more details see Shenhav et al. 2019, Nature Methods (https://www.nature.com/articles/s41592-019-0431-x).

Support
-----------------------

For support using FEAST, please email: liashenhav@gmail.com


Software Requirements and dependencies
-----------------------

FEAST is written R. In addition to R 3.4.4 (and higher), it has the following dependencies::

"doParallel", "foreach",  "dplyr", "vegan", "mgcv", "reshape2", "ggplot2", "philentropy", "MCMCpack", "lsei", "Rcpp", "RcppArmadillo" and "cowplot".


Input format
-----------------------
The input to FEAST is composed of two tab-separated ASCII text files :

count table  - A matrix of samples by taxa with the sources and sink. The first row contains the sample ids ('SampleID'). The first column contains taxa ids. Then every consecutive column contains read counts for each sample. Note that this order must be respected (see example below).

metadata -  The first row contains the headers ('SampleID', 'Env', 'SourceSink', 'id'). The first column contains the sample ids. The second column is a description of the sampled environment (e.g., human gut), the third column indicates if this sample is a source or a sink (can take the value 'Source' or 'Sink'). The forth column is the Sink-Source id. 
When using multiple sinks, each tested with the same group of sources, only the rows with 'SourceSink' =  Sink will get an id (between 1 -  number of sinks in the data). In this scenatio, the sources ids are blank. 
When using multiple sinks, each tested with a distinct group of sources, each combination of sink and its corresponding sources should get the same id (between 1 -  number of sinks in the data). 
Note that these names must be respected  (see examples below).



Output format
-----------------------

The output is a vector of  contributions of the known and unknown sources (with the pre-defined source environments as headers). 

Usage instructions
---------------------------

FEAST will be available on Qiime II in July 2019. Until then you can easily run it on your computer in just a few easy steps which I will walk you through in the following lines. 

	1. Clone this repository ('FEAST') and save it on your computer.
	2. Save your input data (metadata and count table) in the directory 'Data_files'.
	3. Run the file 'FEAST_main' from 'FEAST_src' after inserting the following arguments as input:


| ARGUMENT | DEFAULT |DESCRIPTION |
| ------------- | ------------- |------------- |
| path  |   |The path in which you saved the repository 'FEAST' (e.g., "~/Dropbox/Microbial_source_Tracking") |
| metadata_file  |   |The full name of you metadata file, including file type (e.g., "my_metadata.txt) |
| count_matrix   |   |The full name of your taxa count matrix, including file type (e.g., "my_count_matrix.txt)  |
| different_sources_flag  |   |Relevant only when using multiple sinks. If you use different sources for each sink, different_sources_flag = 1, otherwise = 0 |
| EM_iterations  | 1000  |Number of EM iterations. We recommend using this default value.   |




Example
---------------------------

To run FEAST on example data (using multiple sinks) do:

	
	1. Clone this repository ('FEAST') and save it on your computer.
	2. Run the file 'FEAST_example_Multiple_sinks' which takes the following arguments as input:
	path = The path in which you saved the repository 'FEAST' (e.g., "~/Dropbox/Microbial_source_Tracking") 
	

Input - 

metadata

*using multiple sinks, each tested with the same group of sources:

| SampleID | Env |SourceSink | id |
| ------------- | ------------- |------------- |-------------|
| ERR525698  |  infant gut 1 | Sink | 1
| ERR525693  |  infant gut 2 | Sink | 2 |
| ERR525688   |  Adult gut 1 | Source| NA |
| ERR525699  |  Adult gut 2 | Source | NA |
| ERR525697  |  Adult gut 3 | Source | NA |


*using multiple sinks, each tested with a different group of sources:

| SampleID | Env |SourceSink | id |
| ------------- | ------------- |------------- |-------------|
| ERR525698  |  infant gut 1 | Sink | 1
| ERR525688   |  Adult gut 1 | Source| 1 |
| ERR525691  |  Adult gut 2 | Source | 1 |
| ERR525699  |  infant gut 2 | Sink | 2 |
| ERR525697  |  Adult gut 3 | Source | 2 |
| ERR525696  |  Adult gut 4 | Source | 2 |


count matrix (first 4 rows and columns):

| | ERR525698 |ERR525693 | ERR525688| ERR525699|
| ------------- | ------------- |------------- |------------- |------------- |
| taxa_1  |  0 | 5 | 0|20 |
| taxa_2  |  15 | 5 | 0|0 |
| taxa_3  |  0 | 13 | 200|0 |
| taxa_4  |  4 | 5 | 0|0 |

 

Output - 

| infant gut 2  |Adult gut 1 | Adult gut 2| Adult gut 3| Adult skin 1 |  Adult skin 2|  Adult skin 3| Soil 1 | Soil 2 | unknown|
| ------------- | ------------- |------------- |------------- |------------- |------------- |------------- |------------- |------------- |------------- |
|  5.108461e-01  |  9.584116e-23 | 4.980321e-12 | 2.623358e-02|5.043635e-13 | 8.213667e-59| 1.773058e-10 |  2.704118e-14 |  3.460067e-02 |  4.283196e-01 |



This is an exmaple illustrating the use of FEAST with multiple sinks. To use FEAST with only one sink, please see 'FEAST_example.R'

© 2018 Big Data and Genomics Lab at UCLA All Rights Reserved
