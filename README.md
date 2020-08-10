# ABCE
The ABCE script pack is designed to calculate A/B-compartment in complicated cases such as the ruble-configuration, large invesrions, hard heterchromatin blocks etc.

At first, you should define to be your case very complicated or not to be.

NO COMPLICATED CASE

If you think the problem with compartment calculation is caused only by large inversion, or some hard heterohromatin blocks, or some local problem with genome assembly, we can enhanced the contact matrix by cutting off all problems.
Than:
1) Dump observed/expected contact matrix in density format (for more details see: https://github.com/aidenlab/juicer/wiki/Data-Extraction).

java -jar juicertools.jar dump oe KR -d you_hic_map_name.hic chr_name chr_name BP resolution path_output
ATTENTION! ABCE works only with intrachromosomal contacts!

2) Use cropping_enhancing.py on dumped observed/expected matrix (use -h to more details).

python cropping_enhancing.py -i path_to_your_matrix -o path_to_output_directory -l start_locus end_locus -r matrix_resolution_in_bp

The script output named as "your_matrix.cropped.prs" contains a pearson correlation of contact matrix within locus of interest.

THE HARD CASE

A/B-calculation


Second stage is calculating PC from scripts output. 

eig_CE.r - generates PC1 from local frames.

eig_CR.r -  generates PC1/2/... from full matrixes.
