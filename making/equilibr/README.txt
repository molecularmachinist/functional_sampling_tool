0. Copy charmm structures and change h-bonds to all-bonds in all mdp files, except minimization

1. remove water with gmx trjconv
2. change z-size of box vector manually
3. solvate and add ions
4. make an index file for gromping


Commands (in solv/):
gmx make_ndx -f ../charmm_model/step5_input.gro <<___EOF___
18 | 19 | 20
name 21 Water_and_ions
!21

q
___EOF___
echo 22 | gmx trjconv -f ../charmm_model/step5_input.gro -s ../charmm_model/step5_input.gro -n index.ndx  -o unsolv.gro
#*** Manually remove som length from z box vector (~2.5nm-> 20.00000) ***
echo Protein System | gmx trjconv -f unsolv.gro -s unsolv.gro -center -o new_box.gro
gmx solvate -cp unsolv.gro -p topol.top -o solv.gro
./sol_remover.py -f input.dat
#*** Manually change SOL to TIP3 in topol.top ***
gmx grompp -f ions.mdp -c solv_removed.gro  -p topol.top -o ions.tpr -maxwarn 1
echo TIP3 | gmx genion -s ions.tpr -o system_solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral -conc 0.15

# We make an index file that has the correct groups, make sure manually this chooses right stuff
gmx make_ndx -f system_solv_ions.gro -o index_grompp.ndx <<___EOF___
18 | 19 | 20
name 21 Water_and_ions
!21
22 & !17
name 23 SOLU
17
name 24 MEMB
21
name 25 SOLV
23 | 24
0
name 27 SYSTEM
del 0-22

q
___EOF___
