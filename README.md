# CORNAS #
CORNAS (Coverage-dependent RNA-Seq) is a fast method that is able to reliably make DEG calls in non-replicated conditions. CORNAS works on gene lists with number of mapped reads (counts). CORNAS is written in [R](https://www.r-project.org/), and can either be run in the commandline with Rscript, or be run in the R console. The package contains two files:

1. `CORNAS.R`
2. `cornas.config`
	
The program had been tested on R version 3.0.0


## HOWTO ##
There are two ways to use this package:

1. USAGE 1 (Rscript):

		Rscript CORNAS.R <config> <datatable>
	
2. USAGE 2 (R console):

		source("/path/to/CORNAS.R")
		cornas("/path/to/config" , "/path/to/datatable")

The `datatable` input file must contain at least three tab-separated columns; one for the Gene/Transcript/Fragment ID and two for sample counts. The file must not contain any column headers.

Crucial to the run is the preparation of the configuration file. The template is provided in `cornas.config`. The next section explains the configuration file, followed by a section with a prepared example (`example_run` directory).



## Configuration file ##
The configuration file stores the parameters needed for the run. The value for each parameter is written after the parameter name with ":", with one parameter per line. Lines that begin with "#" are ignored by the program. The following are the parameters:

### Compulsory: 
1. **Gene Name** = The column which contains the IDs.
2. **Sample A column** = The column with observed counts for the first sample
3. **Sample B column** =  The column with observed counts for the second sample
	
### Compulsory coverage option choice of either:
1. Option A: If coverage is known.

	+ **Sample A Coverage**
	+ **Sample B Coverage**

2. Option B: With the total read counts for each sample. 
	
The sequencing coverage is the number of total reads (observed counts) divided by the actual amount of fragments present in the PCR mix. While we have calculated it to be close to 300,000,000 (3M), you should calculate the sample coverages in your experiments and add in the the coverages (Option A). Otherwise, the program will calculate the coverage with the detected total reads over 3M per sample from the `datatable` file (Option B).


### Optional:
1. **Alpha** = If you wish to change the alpha from it's default of 99%. Reducing this percentage will increase sensitivity of making DEG calls.
2. **Fold threshold** = Sets the fold cut-off for signifance. Default is 1.5. Reducing this percentage will increase sensitivity of making DEG calls.

For further information on the parameters, please refer to the CORNAS paper.


## Example run ##
In our example, we will conduct a run using data from Marioni *et. al's* study of human kidney and liver tissues in 2008 ([RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays](http://genome.cshlp.org/content/18/9/1509.long)). The data consists of 32000 genes with their counts from two samples of a Illumina Genome Analyzer run compiled in a tab-delimited file called `test4_kidneyliver_example.tab`. Note that there are no headers in this file. Just for your reference, the columns refer to the Gene ID, gene length, Sample R2L2Kidney counts and Sample R2L3Liver counts respectively. For more information on the coverage and data processing, please refer to the CORNAS paper.

We will need to prepare a configuration file,`cornas.config.test4`, to tell CORNAS which two samples to compare, as well as to set the run parameters. We refer to columns by numbers, starting with 1 from the left-most column. In our example, the sample counts are in column 3 and column 4, while the gene ID is in column 1. We had estimated the coverage for each sample as described in our paper as 0.0631 and 0.0681 respectively. The rest of the parameters were left as default. To run with USAGE 1, the full command is:

		Rscript CORNAS.R cornas.config.test4 test4_kidneyliver_example.tab >cornas_test4_example2.out

You can then compare the output above with the one provided in `cornas_test4_example1.out`.

If you wish to run on the R console instead, load CORNAS.R as a source in R first:

		source("/path/to/CORNAS.R")

Given your working directory is in `example_run` and you wish to save the output as `cornasExample1`, run the following command:

		cornasExample1 <- cornas("cornas.config.test4" , "test4_kidneyliver_example.tab")

You may than print `cornasExample1` to file, or do further processing in R.


## References ##
Ihaka, Ross, and Robert Gentleman. "R: a language for data analysis and graphics." Journal of computational and graphical statistics 5.3 (1996): 299-314.

Marioni, John C., et al. "RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays." Genome research 18.9 (2008): 1509-1517.


## AUTHORS & COPYRIGHT ##
CORNAS was developed by Joel Low Zi-Bin, Khang Tsung Fei & Martti Tapani Tammi.

Copyright is under MIT license (see LICENSE.txt).
