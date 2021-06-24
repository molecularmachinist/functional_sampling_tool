#!/usr/bin/env python3

from optparse import OptionParser
import subprocess as subp

def optP():
    """Parse command line options.
    """

    usage="[python3] %prog -f input.dat"
    description="Membrane solvent remover"
    version="\n%prog Version 1.0 \n\nRequires Python 3\n" \

    optParser = OptionParser(
                        usage=usage,
                        version=version,
                        description=description
    )

    optParser.set_defaults(
                    param_file=None,
                    input_file=None,
                    output_file=None,
                    res_name=None,
                    first_atom=None,
                    res_size=None,
                    lipid_name=None,
                    lim_atom=None,
                    adjust_low=None,
                    adjust_high=None
    )

    optParser.add_option(
                    '-f', type='str',
                    dest='param_file',
                    help="Parameter file (e.g. input.dat)"
                    " [default: %default]"
    )

    optParser.add_option(
                    '-i', type='str',
                     dest='input_file',
                     help="Input 3D structure (gro)"
                     " [default: %default]"
    )

    optParser.add_option(
                    '-o', type='str',
                     dest='output_file',
                     help="Output 3D structure (gro)"
                     " [default: %default] "
    )

    optParser.add_option(
                    '--rn', type='str',
                     dest='res_name',
                     help="Name of residues to be removed"
                     " [default: %default]"
    )

    optParser.add_option(
                    '--ra', type='str',
                     dest='first_atom',
                     help="First atom of residue rn"
                     " [default: %default] "
    )

    optParser.add_option(
                    '--rs', type='int',
                     dest='res_size',
                     help="Size of the residue rn"
                     " [default: %default]"
    )

    optParser.add_option(
                    '--ln', type='str',
                     dest='lipid_name',
                     help="Name of the lipid from which to calculate the membrane plane"
                     " [default: %default]"
    )

    optParser.add_option(
                    '--la', type='str',
                     dest='lim_atom',
                     help="Atom name of ln used for membrane plane calculation"
                     " [default: %default] "
    )

    optParser.add_option(
                    '--al', type='float',
                     dest='adjust_low',
                     help="Manual adjusting of lower leaflet"
                     " [default: %default]"
    )

    optParser.add_option(
                    '--ah', type='float',
                     dest='adjust_high',
                     help="Manual adjusting of upper leaflet leaflet"
                     " [default: %default]"
    )

    optParser.add_option(
                    '--top', type='str',
                     dest='topol',
                     help="Topology file"
                     " [default: %default]"
    )

    return optParser.parse_args()

def read_n_split(f):
    return f.readline().split()

def skip_lines(f, i):
    for i in range(i):
        next(f)

def tail(f):
    return subp.getoutput('tail -1 ' + f) + "\n"

def read_gro_line(f_in):
    l = f_in.readline()
    if len(l)<45:
        return ''

    # res_num = int(l[0:5])
    # a_number = int(l[15:20])
    # x = float(l[20:28])
    # y = float(l[28:36])

    res_name = l[5:10].strip()
    a_name = l[10:15].strip()
    z = float(l[36:44])

    return res_name, a_name, z, l

def read_params(fname):
    f = open(fname, 'r')
    params = []

    skip_lines(f,1)
    split = read_n_split(f)
    params.append(split[0])
    params.append(split[1])

    skip_lines(f,2)
    split = read_n_split(f)
    params.append(split[0])
    params.append(split[1])
    params.append(int(split[2]))

    skip_lines(f,2)
    split = read_n_split(f)
    params.append(split[0])
    params.append(split[1])

    skip_lines(f,2)
    split = read_n_split(f)
    params.append(float(split[0]))
    params.append(float(split[1]))

    skip_lines(f,2)
    split = read_n_split(f)
    params.append(split[0])

    f.close()
    return params

def change_mol_number(output_file, removed_lines):
    f_in = open("temp_out","r")
    f_out = open(output_file, "w")

    # Handling the first two lines of the .gro file ------------------
    f_out.write(f_in.readline())
    f_out.write("%d\n"%(int(f_in.readline().split()[0])-removed_lines))
    #  ---------------------------------------------------------------

    for line in f_in:
        f_out.write(line)

    f_in.close()
    f_out.close()

    subp.run(["rm", "temp_out"])

def change_top_number(output_file, removed_mols, resname):
    print("Changing number of molecules in %s"%output_file)
    f_out = open("temp_out","w")
    f_in = open(output_file, "r")

    line = ""

    while not "[ molecules ]" in line:
        line = f_in.readline()
        f_out.write(line)

    for line in f_in:
        if(line.split()[0]==resname):
            parts = line.split()
            parts[1] = str(int(parts[1])-removed_mols)
            f_out.write("{:<15s} {:7s}\n".format(*parts))
            break
        else:
            f_out.write(line)

    for line in f_in:
        f_out.write(line)

    f_in.close()
    f_out.close()

    subp.run(["mv", output_file, "#%s.sr#"%output_file])
    subp.run(["mv", "temp_out", output_file])


def find_min_n_max(input_file, lipid_name):
    f_in = open(input_file, 'r')
    skip_lines(f_in,2)
    first_meet = True

    while True:
        try:
            resname, atom, z, line = read_gro_line(f_in)
        except ValueError:
            break

        if resname == lipid_name:
            if first_meet:
                minim = z
                maxim = z
                first_meet = False
            else:
                if z < minim:
                    minim = z
                elif z > maxim:
                    maxim = z

    return minim, maxim

def find_center(input_file, res_name, atom_name):
    f = open(input_file, 'r')
    skip_lines(f,2)

    z_sum = 0
    z_numb = 0

    while True:
        try:
            res, atom, z, line = read_gro_line(f)
        except ValueError:
            break

        if res == res_name and atom == atom_name:
            z_sum += z
            z_numb += 1

    f.close()

    try:
        return z_sum/z_numb
    except ZeroDivisionError:
        print("\nNo matches found. Check lipid_name and lim_atom.\n")
        raise

def leaflets(input_file, lipid_name, lim_atom):
    f = open(input_file, "r")
    skip_lines(f,2)

    upper_zsum = 0
    upper_num = 0
    lower_zsum = 0
    lower_num = 0

    center = find_center(input_file, lipid_name, lim_atom)

    while True:
        try:
            res, atom, z, line = read_gro_line(f)
        except ValueError:
            break

        if res == lipid_name and atom == lim_atom:
            if z > center:
                upper_zsum += z
                upper_num += 1
            else:
                lower_zsum += z
                lower_num += 1

    f.close()

    try:
        return lower_zsum/lower_num, upper_zsum/upper_num
    except ZeroDivisionError:
        print("\nNo matches found. Check lipid_name and lim_atom.\n")
        raise

def remove_lims(
    input_file, output_file, zmin, zmax, res_name, first_atom, res_size, topol
):

    f = open(input_file, "r")
    out = open("temp_out", "w")

    out.write(f.readline())
    out.write(f.readline())

    removed_lines = 0

    while True:
        try:
            resname, atom, z, line = read_gro_line(f)
        except ValueError:
            break

        if resname == res_name:
            if z < zmax and z > zmin and atom == first_atom:
                skip_lines(f,res_size-1)
                removed_lines += res_size
                continue

        out.write(line)

    out.write(tail(input_file))

    f.close()
    out.close()

    change_mol_number(output_file, removed_lines)
    if(topol):
        change_top_number(topol, removed_lines//res_size, res_name)

    print(
    "Removed: %d molecules = %d atoms"%(removed_lines//res_size, removed_lines)
    )

def remove(
    input_file, output_file, res_name, first_atom, res_size, lipid_name,
    lim_atom, adjust_low, adjust_high, topol
):

    lower, upper = leaflets(input_file, lipid_name, lim_atom)
    lower += adjust_low
    upper -= adjust_high

    remove_lims(
        input_file, output_file, lower, upper, res_name, first_atom, res_size, topol
    )

def read_n_remove():
    opt, args = optP()

    params = [
        "input_file",
        "output_file",
        "res_name",
        "first_atom",
        "res_size",
        "lipid_name",
        "lim_atom",
        "adjust_low",
        "adjust_high",
        "topol"
    ]

    params_read = [None]*10

    if opt.param_file != None:
        params_read = read_params(opt.param_file)

    for i in range(len(params)):
        opt_dict = opt.__dict__

        if opt_dict[params[i]] != None:
            params_read[i] = opt_dict[key]

    if None in params_read:
        print("Some params missing.\n")
        sys.exit(1)

    remove(*params_read)

if __name__ == "__main__":
    read_n_remove()
