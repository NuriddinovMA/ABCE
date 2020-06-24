# ABCE

First stage is the enhancing of contact matrixes. This may be doed by two ways: 

1) cropped_pc.py - processing of OE matrix with cropping heterochomatine
2) contrast_enhancing.py - processing of OE matrix with cropping heterochomatine and contrast enhancing

This scripts require a observed/expected matrix in dencity format. They analyze only interchromosome contacts. 

If you use juicertools to generate hi-c matrix, we can dump required matrix:
java -jar juicertools.jar dump oe KR -d you_hic_map_name.hic chr_name chr_name BP resolution path_output

Second stage is calculating PC1 from scripts output. 

to_r.py - r-script manager

eig_CE.r - generates contrast enhanced PC1

eig_CR.r -  generates PC1 from cropped matrixes

eig_CR2.r - generates PC2 from cropped matrixes

eig_FR.r - generates PC1 from cropped matrixes with framing

