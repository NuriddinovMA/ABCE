# ABCE
The ABCE scripts are designed to calculate A/B-compartments in complicated cases such as the Rabl-configuration of chromosomes, 
chromosomes with large inversions, heterchromatin blocks and etc.

ABCE consists of 2 steps, each step has 2 options: 

1) Data preparation: 
    * easy option
    * complicated option
2) A/B-calculation
    * local (framed) approach 
    * whole-chromosome approach 

<details>
  <summary> Is my case easy or complicated? </summary>
  
* easy: you think the problem with compartments calculation is caused by large 
inversion, heterochromatin blocks or mis-assemblies

* complicated: you have Rabl-orientation of chromosomes, or solution 
 under "Easy CASE" was not efficient
</details>

You can use any combionation of the options from these two steps. The best strategy depends on your input data and should be 
determined empirically. 

## Requerments:
* java 8 with juicer_tools (for contacts extraction; https://github.com/aidenlab/juicer/wiki/Data-Extraction)
* R 3.4.4
* python 2.7
* numpy
* matplotlib for picture output

## Quick usage guide

## STEP I - Data preparation

### Option 1. EASY CASE
Manually define region of interest: single contiguous locus (free of inversions/misassemblies/heterochromatin blocks) which should be split into A/B compartments  
```bash
# Dump observed/expected contact matrix in dense format (for more details see: https://github.com/aidenlab/juicer/wiki/Data-Extraction).
java -jar juicertools.jar dump oe KR -d you_hic_map_name.hic chr_name chr_name BP resolution path_to_your_oe_matrix

# Crop region of interest   
python cropping_enhancing.py -i your_oe_matrix -l locus_start locus_end -r resolution
```

Parameters:
* your_oe_matrix - the data obtained from `juicertools.jar dump` command

* resolution - resolution used for `juicertools.jar dump` command

* locus_start and locus_end - genomic coordinates of region of interest

Use `python cropping_enhancing.py -h` for more details of usage and 
https://github.com/aidenlab/juicer/wiki/Data-Extraction for details of `juicertools.jar dump` parameters.

The script will produce *your_matrix.cropped.prs* file
 
This file should be used in the next step **A/B-CALCULATION**


### Option 2. COMPLICATED CASE
   
```bash
# Dump observed/expected contact matrix in dense format (for more details see: https://github.com/aidenlab/juicer/wiki/Data-Extraction).
java -jar juicertools.jar dump oe KR -d you_hic_map_name.hic chr_name chr_name BP resolution your_oe_matrix

# Enhance oe-matrix
python contrast_enhancing.py -i your_oe_matrix -l locus_start locus_end -r resolution 
-d distance -c radius -e ce_radius
```
Parameters:

* your_oe_matrix - the data obtained from `juicertools.jar dump` command

* resolution - resolution used for `juicertools.jar dump` command

* locus_start locus_end - coordinates of the locus of interest (i.e. 1000000 99000000).

* distance, radius, ce_radius - see _example_command.txt_ for default values and [parameters optimization](#params) section to find how to choose 
distance, radius and ce_radius 
 
The script will produce *your_matrix.[ce_radius].ce.prc.prs* and *your_matrix.[ce_radius].ce.range.prs* files.


Each of these files could be used in the next step (A/B-CALCULATION). 
We advice to proceed with both and manually check which will provide better results.

## STEP II - A/B-CALCULATION

### Option 1. Compute whole-matrix PCs
```bash
# Compute Principle components
Rscript eig_whole.r --args your_enhanced_matrix PC_num resolution chrm_name locus_start locus_end
```

Parameters:
* your_enhanced_matrix - should be one of the files produced in the previous steps 
(your_matrix.cropped.prs/your_matrix.[ce_radius].ce.prc.prs/your_matrix.[ce_radius].ce.range.prs)

* resolution - resolution used for `juicertools.jar dump` command

* chrm_name - chromosome name to write in output file in a .bedGraph format

* locus_start locus_end - coordinates of the locus of interest (i.e. 1000000 99000000).

* PC_num - number of PC to compute. For A/B-compartments you need 1st PC, so this should be set to 1


The resulting file will contain desired PC1 vector: *your_enhanced_matrix.whole.pc1.eig.bedGraph* 

### Option 2. Compute PCs within local frames

`Rscript eig_framed.r your_enhanced_matrix track_for_correlation.bedGraph resolution frame_length chrm_name locus_start locus_end`

Parameters:
* your_enhanced_matrix - should be one of the files produced in the previous steps 
(your_matrix.cropped.prs/your_matrix.[ce_radius].ce.prc.prs/your_matrix.[ce_radius].ce.range)

* resolution - resolution used for `juicertools.jar dump` command

* locus_start locus_end - coordinates of the locus of interest (i.e. 1000000 99000000).

* frame_length - *eig_CE.r* generates PC1 within local frames. frame_length corresponds to the lengths of 
these frames in bp  

* track_for_correlation.bedGraph - Values of PC1 calculated within individual frames should 
be correlated with external standard. This standard (.badGraph format) should reflect chromatin state, i.e. represent
GC-content, gene density or transcriptional activity. See [parameters optimization](#params) and [usage notes](#use_notes) for details.

The resulting file will contain desired PC1 vector: *your_enhanced_matrix.framed.pc1.eig.bedGraph*

## <a name="params"></a> Parameters optimization

The ABCE pipline requeres oprimization of several parameters.

**1.) paramateres of `contrast_enhancing.py` script: distance and radius**

This script perfroms averaging of contacts within submatrices of Hi-C matrix (matrix smoothing).
Depending on data quality (sequencing deapth and noise level) and distance from diagonal the optimal
submatrix size could be different. 

Submatrix size is defined by radius parameter, which could be different at different distances from
diagonal. For example, near main diagonal, where the data is reach, we use radius=0, which means
that no averaging will be performed. For larger distances, we increase radius. Paired distance-radius 
values could be passed to `contrast_enhancing.py` as follows:
`python contrast_enhancing.py -d distance1 distance2 distance3 -c radius1 radius2 radius3 ...other arguments...` 

We provide a script `error_estimation.py` which estimates the radius for different distances.
Use this script as follows:

```bash
# first, get NOT NORMALIZED (raw) data:
java -jar juicertools.jar dump observed NONE -d you_hic_map_name.hic chr_name chr_name BP resolution your_raw_matrix
# then, run error_estimation.py on the obtained raw matrix
python error_estimation.py -i your_raw_matrix -o output -t threshold
```

The _threshold_ determines radius _c_ so that for each point _i,j_ of the raw Hi-C matrix, 
the submatix _i-c..i+c ; j-c..j+c_ contains >=_threshold_ reads. 

Averaging/smoothing over larger submatrices (larger _threshold_ values) reduces noise, but will not allow inferring fine compartments structure. 
Using smaller _threshold_ would allow capturing smaller compartments, but is more sensitive to noise.
We recommend starting with _threshold_ value equal to 25 at 25 KB resolution.
  
The output of `error_estimation.py` looks like this:

```
example/example.25000.none
distance radius
0 0
46 1
510 2
511 3
512 4
```

This means that desirable radius for contacts averaging for loci lying between 0 and 46 bins is 0, 
between 46 and 510 bins is 1, etc.
 
This information is a advise, not a requirement.
 
After certain distance the radius often increases only slightly, so at this point we fix radius value to the 
maximal number observed in the file. In the example above, we would use radius 0 for distances 
below 46 bins, 1 for 46..510 and 4 for distances >510 bins   
   
**2.) paramater of `contrast_enhancing.py` script: ce_radius**

Within any submatrix of size 2*ce_radius+1, all correlation values will be normalized to be 
in the -1..+1 range. Thus, ce_radius should be high enough to include several compartment transistion.
On the other hand, large ce_radius will allow long-range effects, such as telomere/centromere 
clustering dominate local compartments.

We typically use ce_radius values (measured in bins) 3,5 or 7, but users may test other values as well.
For convenience, multiple space-separated values could be used as an input of `contrast_enhancing.py`:

`contrast_enhancing.py -e 3 5 7 ...other arguments...`

This will produce multiple output files, each corresponding to one of the provided ce_radius values.

**3.) paramater of `eig_framed.r`: track_for_correlation**

In our experience, the highest correlation is obtained when using RNA-seq data (fpkm) .bedGraph. 
However, gene density (calculated per Hi-C bin) also works well. GC-contents gives very low correlation
values and should not be used if other options available.

## <a name="use_notes"></a> Usage notes

0) ABCE works only with intrachromosomal contacts!
1) The track for correlation must contain ONLY a chromosome of locus of interest.
2) The track must start from 0.
3) The end of track must be equal the locus end or more.
4) There are not missing bins.
