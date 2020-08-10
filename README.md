# ABCE
The ABCE script pack is designed to calculate A/B-compartment in complicated cases such as the ruble-configuration, large invesrions, hard heterchromatin blocks etc.

At first, you should define to be your case very complicated or not to be.

NO COMPLICATED CASE

If you think the problem with compartment calculation is caused only by large inversion, or some hard heterohromatin blocks, or some local problem with genome assembly, we can enhanced the contact matrix by cutting off all problems.
Than:
1) Dump observed/expected contact matrix in density format (for more details see: https://github.com/aidenlab/juicer/wiki/Data-Extraction).

java -jar juicertools.jar dump oe KR -d you_hic_map_name.hic chr_name chr_name BP resolution path_to_your_oe_matrix
ATTENTION! ABCE works only with intrachromosomal contacts!

2) Use cropping_enhancing.py on dumped observed/expected matrix (use -h to more details).

python cropping_enhancing.py -i path_to_your_oe_matrix -o path_to_output_directory -l locus_start locus_end -r matrix_resolution_in_bp

The script output named as "your_matrix.cropped.prs" contains a pearson correlation of contact matrix within locus of interest. This matrix are given as parameter to eig_CE.r or eig_CR.r

THE HARD CASE

If you think the yor case is complicated, for example, been caused by rubl-orientaion, or easy way was unefficient, we should use the scripts error_estimation.py and contrast_enhancing.

In many case? the problem with calculating of A/B-compartments are resulted from a weak contrast between A- and B-compartments. This can be caused by rubl-orientation of chromosomes or a noise. To avoid this problem, our script combines and smoothes contacts between distant loci in relation of distance.

1) Dump raw contact matrix in density format (for more details see: https://github.com/aidenlab/juicer/wiki/Data-Extraction).

java -jar juicertools.jar dump observed NONE -d you_hic_map_name.hic chr_name chr_name BP resolution path_to_your_raw_matrix
ATTENTION! ABCE works only with intrachromosomal contacts!

2) Estimating of minimal radius for combining and smoothing contacts (use -h and see ExampleCommand.txt to more details).
python error_estimation.py -i path_to_your_raw_matrix -o output -t threshold

The "threshold" is a number of minimal contacts in the area (mean for given distance).
The output looks like this:

example/example.25000.none
46 1
510 2
511 3
512 4

That means the desirable radius of contact combining for loci lying between 0 and 46 bins is 0, between 46 and 510 bins is 1, etc. This information is a advise, not a requirement. 

3) Enhancing matrix using information about distance and combine radius (use -h and see ExampleCommand.txt to more details).

python contrast_enhancing.py -i path_to_your_oe_matrix -o path_to_output_directory -l locus_start locus_end -r matrix_resolution_in_bp -d distance -c combine_radius -e contrast_enhancing_radius

The script outputs are named as "your_matrix.[combine_radius].ce.prc.prs" and "your_matrix.[combine_radius].ce.range.prs".

If you ignored step 2 you can ignore -d and -c flags.

A/B-CALCULATION

Second stage is calculating PC from scripts output. 

1) eig_CR.r generates PC1/2/... from full matrix.

Using (see ExampleCommand.txt to more details):
eig_CR.r --args your_enhanced_matrix number_of_PC resolution_in_bp chrm_name locus_start_in_bp locus_end_in_bp 
if you wish calculate A\B-compartment, than number_of_PC is 1.

2) eig_CE.r - generates PC1 from local frames.
This script splits locus of interest on several frames and calculate independly PC1 for each frame. Than values of PC1 are correlated with given track (.badGraph) and the data from different frame are combined. Than operation help us to avoid problems with distant interaction.

Using  (see ExampleCommand.txt to more details):
eig_CE.r --args your_enhanced_matrix track_for_correlation.bedGraph resolution_in_bp frame_length_in_bp chrm_name locus_start_in_bp locus_end_in_bp 

The output is named as your_enhanced_matrix.framed.eig.bedGraph

ATTENTION! 
1) The track for correlation must contain ONLY a chromosome of locus of interest.
2) The track must start from 0.
3) The end of track must be equal the locus end or more.
4) There are not missing bins.
