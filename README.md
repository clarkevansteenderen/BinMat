## **BINMAT: FOR FRAGMENT ANALYSIS DATA** 
[![DOI](https://zenodo.org/badge/238669382.svg)](https://zenodo.org/badge/latestdoi/238669382)   

### **USER GUIDE**

<img src="https://github.com/clarkevansteenderen/BinMat/blob/master/www/clevercow.png" height = 120>

---

*Created by:*   
*Clarke van Steenderen*   
*Department of Zoology and Entomology*   
[*The Centre for Biological Control*](https://www.ru.ac.za/centreforbiologicalcontrol/)   
*Rhodes University, Grahamstown, Eastern Cape, South Africa*   
*e-mail:* vsteenderen@gmail.com   

---

**OVERVIEW**

This program was written to rapidly consolidate binary matrices derived from dominant marker genetic analyses (such as for 
inter-simple sequence repeats (ISSRs) and amplified fragment length polymorphisms (AFLPs)). 
This can be applied to data consisting of replicate pairs, or it can be used to consolidate the entire dataset (i.e. treating all the samples uploaded as replicates for one
sample). 
The user is also able to generate interactive non-metric MDS plots, and filter their data according to a specified minimum peak number threshold.		
Please go through the [vignette](https://cran.r-project.org/web/packages/BinMat/vignettes/BinMat.html) for more information and example output. 

---

**INSTALLING THE R PACKAGE**

BinMat is available on [CRAN](https://cran.r-project.org/web/packages/BinMat/index.html), and can be installed by using the command: 

```{r}
install.packages("BinMat")
```

Or it can be run via GitHub as an R shiny application:   

```{r}
install.packages("shiny")      
library(shiny)      
shiny::runGitHub("BinMat", "clarkevansteenderen")    

```

Or accessed via the R Shiny Apps platform:   
https://clarkevansteenderen.shinyapps.io/BINMAT/   

To cite BinMat, use:

```{r}
citation("BinMat")
```

---

**UPLOADING DATA**

Binary data needs to be in the following example format, saved as a .csv file:

---

|                      		| 	Locus 1     	  | Locus 2   			| Locus 3   		| Locus 4 	    	| Locus 5 	    	|
|----------------------		|:------------:	      | :------------:	    |:------------:	    |:------------:	    |:------------:	    |
| Sample A replicate 1 		| 	0       	      |  	0       	    | 	1       	    | 	1       	    | 	1       	    |
| Sample A replicate 2 		| 	0       	      |  	0       	    | 	1       	    | 	1       	    | 	1       	    |
| Sample B replicate 1 		| 	1       	      |  	1       	    | 	0       	    | 	0       	    | 	0       	    |
| Sample B replicate 2 		| 	0       	      |  	1       	    | 	0       	    | 	0       	    | 	1       	    |

---

Note that replicate pairs must be directly underneath each other, and that each sample needs to have a unique name.

The following conditions are applied to binary data replicates:
 1. A 0 and a 1 produces a "?"
 2. A 0 and a 0 produces a "0"
 3. A 1 and a 1 produces a "1"    
 
The consolidated output for the above example would thus be:

Sample A: 0 0 1 1 1

Sample B: ? 1 0 0 ?

---

In the case of all the rows being replicates for one sample, an average is calculated for each locus column. If this average is > 0.5, a '1' is recorded. 
If it = 0.5, a '?', and < 0.5, a '0' is recorded. For example:

|                      		| 	Locus 1     	  | Locus 2   			| Locus 3   		| Locus 4 	    	| Locus 5 	    	|
|----------------------		|:------------:	      | :------------:	    |:------------:	    |:------------:	    |:------------:	    |
| Sample A replicate 1 		| 	0       	      |  	0       	    | 	1       	    | 	1       	    | 	1       	    |
| Sample A replicate 2 		| 	0       	      |  	0       	    | 	1       	    | 	1       	    | 	1       	    |
| Sample A replicate 3 		| 	1       	      |  	1       	    | 	0       	    | 	0       	    | 	0       	    |
| Sample A replicate 4 		| 	0       	      |  	1       	    | 	0       	    | 	0       	    | 	1       	    |
																	
---

The output would be:

Sample A: 0 ? ? ? 1
																	
---

**FUNCTIONALITY: R SHINY APPLICATION**

Once the .csv file is uploaded, the 'Preview' button can be clicked to view the data. 

The 'Check my data for unwanted values' button can be clicked to ensure that the data contains only 1's, 0's and ?'s.

The 'Consolidate Matrix' button consolidates replicates, and displays the output on the screen. 

The 'Download Matrix' button enables the user to download the consolidated matrix to their PC. It is important to name the file with a .csv extension. E.g. MyISSRData.csv

**Summary Tab:**

The 'Summary info' button provides information regarding the average number of peaks and its standard deviation, the 
maximum and minimum number of peaks, and the total number of loci.
The output table can be saved as a .csv file.

**Error Rates Tab:**

The 'Check Error rates' button calculates the average Euclidean error rate and its
standard deviation, as well as the average Jaccard error rate and its standard deviation.
This information can be saved as a .csv file.

The 'Remove' button removes samples with a Jaccard error equal to and above the specified value (between 0 and 1),
and recalculates the error rates. It also records the percentage of samples removed.

**UPGMA Tree Tab:**

The 'Create Clustering Tree' creates a UPGMA hierarchical clustering tree based upon a distance matrix calculated using the Jaccard
index.  The user can specify the number of bootstrap replicates to run, and the resulting plot can be downloaded as a .svg file.
When prompted for a name to save the file as, ensure that the .svg extension follows the name given (i.e. tree01.svg).

**nMDS Plot Tab:**

This is for the creation of a non-metric multi-dimensional scaling (nMDS) plot. 
Upload a binary matrix that has already been consolidated. Ensure that the first column contains the sample name, and that the second column contains grouping information.

For example:

---

|           |    Group        	| 	Locus 1     	  | Locus 2   			| Locus 3   		| Locus 4 	    	| Locus 5 	    	|
|-----------|:----------:		|:------------:	      | :------------:	    |:------------:	    |:------------:	    |:------------:	    |
| Sample A  |    Africa	    	| 	?       	      |  	0       	    | 	1       	    | 	1       	    | 	1       	    |
| Sample B  |    Asia	    	| 	0       	      |  	0       	    | 	1       	    | 	1       	    | 	?       	    |
| Sample C  |    Europe	    	| 	1       	      |  	?       	    | 	0       	    | 	0       	    | 	0       	    |
| Sample D  |    USA	        | 	?       	      |  	1       	    | 	0       	    | 	?       	    | 	1       	    |
                            

---
							
An editable table will appear, displaying your group names as well as default colours and character shapes. Change these if desired.

Select the desired number of dimensions for the nMDS plot. Only 2 and 3 are available for selection.

Click the "Plot nMDS" button to display the plot. Point size can be altered by adjusting the slide bar, and sample labels can be displayed or hidden.
If any changes are made, you need to click the "Plot nMDS" button again to incorporate them. 
						
							
**nMDS Validation Tab:**

The Scree plot shows a red line at a stress value of 15%. K-values below this are acceptable to implement.
If a K-value of both 2 and 3 yield stress values above 15%, a different ordination method should be used.
The Shepard plot gives an indication of how well the ordination distances capture the original distances in the data. 
A high R-squared value is desirable with points that cluster around the line of best fit.

**Filter data Tab:**

You can remove samples from your dataset with a peak count less than a specified value. The 'CHECK' button tells you how many samples were removed,
and the 'download filtered' and 'download removed' samples buttons allows you to download the relevant data as .csv files.
The filtered data can then be re-uploaded in the 'nMDS PLOT' tab, and a new plot created.

---

**USING THE R PACKAGE**

	library(BinMat)
	data = read.csv("file_name.csv") # read in the binary data for all replicate sample pairs
	check.data(data) # check for any unwanted characters apart from zeros and ones. Here, 'none found' is printed
	peaks.original(data) # access peak summary information
	cons = consolidate(data) # consolidate all replicate pairs into consensus reads
	errors(cons) # get the mean Jaccard and Euclidean error rates for all replicate samples
	peaks.consolidated(cons) # access information about the peak numbers in the consolidated matrix
	filtered_data = peak.remove(cons, thresh = 10) # if desired, remove samples with peak counts less than 10
	write.csv(filtered_data, file = "file_name_2.csv") # write the consolidated matrix (either the original matrix, or the filtered one) to the working directory, open in Microsoft Excel, and add grouping information as a second column in the matrix
	data_consolidated = read.csv("file_name_2.csv") # read in this consolidated file, with the added grouping information
	scree(data_consolidated) # choose a k-value below the 15% threshold
	shepard(data_consolidated, k_val = 2) # create a shepard plot for k = 2 or k = 3 dimensions
	group.names(data_consolidated) # access the group names
	colrs = c("red", "blue", "green", "gold") # assign appropriate colors to each group. This is an example for 4 groups.
	shps = c(16, 16, 16, 16) # assign point shapes to each group
	nmds(data_consolidated, colours = colrs, shapes = shps, labs=F, pt_size = 2) # create the nMDS plot with, or without labels
	upgma(data_consolidated, fromFile = T, size = 1, bts = 500) # plot a UPGMA hierarchical clustering tree. fromFile = T, because this data was read in again from Excel after adding the grouping column

--- 

**METHODS**

Mismatch/Euclidean error = (f01 + f10) / (f01 + f10 + f00 + f11)

The Euclidean error rate includes the shared absence of a band (f00).

Jaccard error = (f01 + f10) / (f01 + f10 + f11)

The Jaccard error rate does not take shared absences of bands into account. The error rate will thus be inflated compared to the Euclidean.

In the formulae above, 'f' refers to the frequency of each combination (i.e. 'f01 + f10' means the sum of all the occurences of a zero and a one, and a one and a zero).
An error value is calculated for each replicate pair, and an average obtained representing the whole dataset. 

The error rates for the example samples below would be:

**Sample X rep 1: 0 1 1 0**

**Sample X rep 2: 1 0 1 0**

Euclidean error = (1+1) / (1+1+1+1) = 2/4 = 0.5

Jaccard error = (1+1) / (1+1+1) = 2/3 = 0.67

**Sample Y rep 1: 1 1 1 0**

**Sample Y rep 2: 1 1 0 0**

Euclidean error = (1) / (1+1+2) = 1/4 = 0.25

Jaccard error = (1) / (1+2) = 1/3 = 0.33

---

Average Euclidean error = (0.5+0.25) / 2 = 0.38

Average Jaccard error = (0.67+0.33) / 2 = 0.5

---

**OTHER R PACKAGES USED BY THIS PROGRAM**

[**pvclust**](https://cran.r-project.org/web/packages/pvclust/index.html)

[**shinyhelper**] (https://cran.r-project.org/web/packages/shinyhelper/index.html)

[**readr**](https://cran.r-project.org/web/packages/readr/index.html)

[**rhandsontable**](https://cran.r-project.org/web/packages/rhandsontable/rhandsontable.pdf)

[**vegan**](https://cran.r-project.org/web/packages/vegan/vegan.pdf)

[**MASS**](https://cran.r-project.org/web/packages/MASS/MASS.pdf)

[**magrittr**](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html)
