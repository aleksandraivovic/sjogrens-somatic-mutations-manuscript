README Mapscape

Example code provided for patient PD42769
Input files in PD42769 folder

Mapscape_files may be needed to run code locally

Inputs needed:
1. Clonal prevalence file:
    sample_id   clone_id     		clonal_prev
    lo0057      PD42769_24     		0.07865169\
    lo0057      PD42769_23 		0.02247191
    lo0057      PD42769a_CD20_8  	0.01123596
    lo0057      PD42769a_CD20_2  	0.00000000
    lo0057      PD42769_20  		0.02247191
    lo0057      PD42769_19  		0.00000000

2. Tree edge file:
 source          	target
 PD42769_16      	PD42769_17
 PD42769_17 		PD42769a_CD20_1
 PD42769_17 		PD42769a_CD20_2
 PD42769_17 		PD42769a_CD20_3
 PD42769_17 		PD42769a_CD20_8
 PD42769_17 		PD42769a_lo0073

3. Image file
Histology image in png format

4. Sample locations file;
X and y pixel coordinates for sample locations in histology image. Can be generated with ImageJ, Photoshop, etc.

5. Clone color table (optional):
Specific color palettes corresponding to each clone and its descendants
  clone_id        	colour 
 PD42769_17      	#66A61E
 PD42769a_CD20_1 	#7FB443
 PD42769a_CD20_2 	#99C369
 PD42769a_CD20_3 	#B2D28E
 PD42769a_CD20_8 	#CCE1B3
 PD42769a_lo0073 	#E5F0D9

For full Mapscape documentation see https://bioconductor.posit.co/packages/3.24/bioc/vignettes/mapscape/inst/doc/mapscape_vignette.html
