#!/usr/bin/env python3
# -*- coding: utf-8 -*-

first_res = 1
last_res  = 826
add_res   = (1301, 1302)

def choose_res(ri):
    """ Helper func for choosing resids
    """
    return (first_res<=ri<=last_res) or (ri in add_res)

def copy_struct(fname,foutname):
    with open(fname) as f:
        with open(foutname, "w") as fout:
            ind=1
            hind=1
            seqr = {}
            # Dict to store indices of nonfiltered residues in the old and new pdbs
            accept = {}
            fout.write(f.readline())
            fout.write(f.readline())
            fout.write(f.readline())
            for line in f:
                if(line[:6] in ("HET   ","HETNAM","HETSYN","FORMUL")):
                    fout.write(line)
                elif(line.startswith("SSBOND")):
                    r1 = int(line[17:21])
                    r2 = int(line[31:35])
                    if(choose_res(r1) and choose_res(r2)):
                        fout.write(line)
                elif(line.startswith("HELIX  ")):
                    r1 = int(line[21:25])
                    r2 = int(line[33:37])
                    if(choose_res(r1) and choose_res(r2)):
                        fout.write("HELIX  %3d%s"%(hind%1000,line[10:]))
                        hind+=1
                elif(line.startswith("SHEET  ")):
                    r1 = int(line[22:26])
                    r2 = int(line[33:37])
                    if(choose_res(r1) and choose_res(r2)):
                        fout.write(line)

                elif(line.startswith("ATOM  ")):
                    ch = line[21]
                    anum = int(line[6:11])
                    resid = int(line[22:26])
                    res = line[17:20].strip()

                    if(choose_res(resid)):
                        fout.write("ATOM  %5d%s"%(ind%100000, line[11:]))
                        accept[anum] = ind
                        ind += 1
                elif(line.startswith("HETATM")):
                    ch = line[21]
                    anum = int(line[6:11])
                    resid = int(line[22:26])
                    res = line[17:20].strip()
                    fout.write("HETATM%5d%s"%(ind%100000, line[11:]))
                    accept[anum] = ind
                    ind += 1
                elif(line.startswith("CONECT")):
                    # Hoping these all come after the ATOM and HETATOM records, so accept is populated
                    ai1 = line[6:11]
                    ai2 = line[11:16]
                    ai3 = line[16:21]
                    ai4 = line[21:26]
                    ai5 = line[27:31]
                    out="CONECT"
                    skip = False

                    for ais in (ai1,ai2,ai3,ai4,ai5):
                        if(ais.strip()):
                            ai = int(ais)
                            if (ai in accept):
                                out += "%5d"%(accept[ai])
                            else:
                                # if any of the atoms is not in the chosen residues, we skip the record
                                skip=True
                                break

                    if(skip):
                        continue
                    else:
                        fout.write("%-80s\n"%out)

copy_struct("5weo.pdb", "modeller/5weo_no_STZ.pdb")
