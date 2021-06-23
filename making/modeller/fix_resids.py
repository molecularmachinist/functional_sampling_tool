#!/usr/bin/env python3
# -*- coding: utf-8 -*-


first_res = 400
last_res  = 841
add_res   = list(range(3533, 3540+1))
# ligand chains, mapped back to match prot chains
chain_map = {"E":"A","F":"B","G":"C","H":"D"}
# And map resdidue numbers by resname to match original
res_map = {"GLU":1301,"CYZ":1302}

def choose_res(ri):
    """ Helper func for choosing resids
    """
    return (first_res<=ri<=last_res) or (ri in add_res)

def copy_struct(fname, hetfname,foutname):
    chosen_res = {}
    resid=1
    prev_resid=-1
    prev_ch=""
    # First we'll read and filter residues
    with open(fname) as f:
        for line in f:
            if line.startswith("ATOM  "):
                ch = line[21]
                anum = int(line[6:11])
                curr_resid = int(line[22:26])
                if(ch!=prev_ch):
                    resid=1
                    prev_ch=ch
                    prev_resid=curr_resid
                elif(prev_resid!=curr_resid):
                    resid+=1
                    prev_resid=curr_resid
                res = line[17:20].strip()

                if(choose_res(resid)):
                    chosen_res[curr_resid] = resid

    # and only on second go we start writing output
    with open(fname) as f:
        with open(foutname, "w") as fout:
            ind=1
            resid=1
            # Dict to store indices of nonfiltered residues in the old and new pdbs
            accept = {}
            pos = f.tell()
            line=f.readline()
            while(line[:6] in ("EXPDTA", "REMARK")):
                fout.write(line)
                pos = f.tell()
                line=f.readline()
            f.seek(pos)

            with open(hetfname) as hetf:
                for line in hetf:
                    if(line[:6] in ("HET   ","HETNAM","HETSYN","FORMUL")):
                        fout.write(line)

            for line in f:
                if(line.startswith("SSBOND")):
                    r1 = int(line[17:21])
                    r2 = int(line[31:35])
                    if(r1 in chosen_res and r2 in chosen_res):
                        fout.write("%s%4d%s%4d%s"%(line[:17],chosen_res[r1],line[21:31],chosen_res[r2],line[35:]))
                elif(line.startswith("ATOM  ")):
                    resid = int(line[22:26])
                    anum = int(line[6:11])

                    if(resid in chosen_res):
                        fout.write("ATOM  %5d%s%4d%s"%(ind%100000, line[11:22],chosen_res[resid]%10000, line[26:]))
                        accept[anum] = ind
                        ind += 1
                elif(line.startswith("HETATM")):
                    ch = line[21]
                    anum = int(line[6:11])
                    resid = int(line[22:26])
                    res = line[17:20].strip()
                    fout.write("HETATM%5d%s%s%4d%s"%
                        (ind%100000, line[11:21], chain_map[ch], res_map[res], line[26:]))
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

if __name__ == '__main__':
    copy_struct("5weo_fill.B99990001.pdb", "5weo_no_STZ.pdb", "5weo_fill.pdb")
