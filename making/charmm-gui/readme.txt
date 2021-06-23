Step-1.1:

-remove CYZ and rename glutamate "chains" to HETA-D


Step 1.2:

-Rename engineered residues: GLU -> GLU (yeah, don't know why, probs the CYZ confusing it)



- Add glycosylation
	HNK-1 GRS below, add to Asn413 of each chain:
	glyc papers:
	https://doi.org/10.1371/journal.pone.0135644
	https://doi.org/10.1016/j.bbagen.2017.06.025


1 BGLCNA
2 - 14B: BGLCNA
3 - - 14B: BMAN
4 - - - 13B: BMAN
5 - - - - 12B: BGLCNA
6 - - - - - 14B: BGAL
7 - - - - - - 13B: BGLCA_3SUF
8 - - - 14B: BGLCNA
9 - - - 16B: BMAN
10- - - - 12B: BGLCNA


( no palmitoylation
	Lipid tails:
	https://doi.org/10.1016/j.neuron.2005.06.035
	CYSP to resid 836 of each chain
)

Step 1.3:

-Use PPM server and send PROA-D only
-flip along z-axis (ppm server guesses intra- and extracellular sides wrong)

Step 2:

- z-length: water thickness 35
- xy: ratios, with size 200
- POPC ratios 1 and 1

Step 3:

- Include 0.15 M NaCl

Step 4.x:

-

Step 5:

- ff: CHARMM36m
- both cation pi thingy and H-mass partitioning
- More charmm minimization
- gromacs
- temperature 310.15
