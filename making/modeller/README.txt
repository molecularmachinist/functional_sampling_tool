1. removed STZ from 5weo and saved to 5weo_no_STZ.pdb (copy_struct.py).
2. using the first function of model-default.py output the sequence in the pdb into 5weo_no_STZ.seq
3. using uniprot alignment to align to GRIA2_HUMAN (P442262(-1))
4. Manually made alignment.ali to have all 4 chains (note only chains A and B have the last Ala in the PDB and "/../../../.." to include ligands)
5. run second function of model-python.py to make the new protein.
6. run fix_resids.py to fix the residue numbering (the previous step outputs a pdb with incremental numbering instead of each chain starting from 1), which also removes ATD and CTD
7. ???
8. profit
