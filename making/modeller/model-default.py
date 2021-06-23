from modeller import *
from modeller.automodel import *    # Load the AutoModel class
# Get the sequence of the 1qg8 PDB file, and write to an alignment file

def get_seq():
    code = '5weo_no_STZ'

    e = Environ()

    # Read in HETATM records from template PDBs
    e.io.hetatm = True
    m = Model(e, file=code)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')



def modelling():
    log.verbose()
    env = Environ()


    # Read in HETATM records from template PDBs
    env.io.hetatm = True

    a = AutoModel(env, alnfile = 'alignment.ali',
                  knowns = "5weo_no_STZ", sequence = '5weo_fill')
    a.starting_model= 1
    a.ending_model  = 1

    a.make()


if __name__ == '__main__':
    #get_seq()
    modelling()
