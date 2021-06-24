Steps to remaking the systems:

TL;DR: modeller->charmm-gui->equilibr->prod

Starting from pdb structure 5weo, use modeller to fit the human sequence.
Following the details in the modeller/ folder, this gives us 5weo_human.pdb.
Use the modeller output for CHARMM-GUI, with details and output under charmm-gui/
Finally copy the gromacs part of the output to equilibr/ and follor the README in there.
