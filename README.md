# Functional Sampling Tool

## Requirements

1. Python 3.7 or higher (3.10 recommended)
1. NumPy and matplotlib (with python<3.10, also importlib_resources)
1. [MDAnalysis](https://docs.mdanalysis.org/stable/index.html)
1. [scikit-learn](https://scikit-learn.org/stable/)
1. GROMACS



Either install them manually, or follow the instructions below to let pip take care of that for you. GROMACS does of course need to be manualy installed.

The tool has only been tested on linux, so it might or might not work with other OS.

## Installation


### Quick and easy

Just run
```sh
pip install .
```
in the project root. This will install all dependencies and compile and install the package.

### Recommended way

In the recommended way we will use conda to make a new envirnonment. This way works independent of which Python version you use by default. If you do not already have it, [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).


First, we will make a new conda environment dedicated just for the tool and activate it (if you have [mamba](https://mamba.readthedocs.io/en/latest/) installed, use it for the first command to speed up the process):

```sh
conda create -c conda-forge -n fst_env python=3.10 numpy matplotlib mdanalysis scikit-learn
conda activate fst_env
```

This will also install the needed dependencies, so that they are handled by conda rather than pip.

Second, we run the command to install the package. 
Here I assume that you have downloaded the project and are within its root folder. From there run

```sh
pip install .
```

When the command finishes, it should be all done and the tool usable as `fst`. If you ever want to uninstall it just run (with the `fst_env` environment activated).

```sh
pip uninstall functional_sampling_tool
```

**Remember** that in this way you need to run `conda activate fst_env` once in every new terminal before using the tool. If you want to use this environment by default, add `conda activate fst_env` at the end of your `.bashrc`


## How it works

### Vocabulary

 1. **Epoch**: the set of N runs that are run at the same time.
 1. **Repetition** (rep): The individual runs within each epoch. i.e. there are N reps in each epoch.

### Epoch choosing

After each epoch, the user defined function is calculated for each frame, and a histogram is calculated, to get the distribution. The tool then smooths the distribution and finds the local minima within the boundaries defined by user. The starting points for each rep in the new epoch will chosen from these minima.

### Clustering

After a set amount of epochs, the tool will start to use a clustering method on the trajectories. This way it will try to identify different states of the structure, that have the same function value. A part of the repetitions will be started from the clusters with the lowest amount of data. The clustering is done by first making a histogram with much larger bins than before, and clustering the data in each bin separately.


## Setup

At minimum you will need to copy the templates for `config.py` and  `sbatch_launch.sh` to your project folder. This can be done with
```sh
fst make_templates
```
after which you should make sure the sbatch launch script matches your workflow and simulation methods. The parts in curly brackets will be populated with the info from `config.py` for email and account, and  the epoch number in the job name.

You will also need to have a starting structure, which should be put into `initil/start.gro`. You will also need a topology in `topol.top`, mdrun options in `mdrun.mdp` and an index file for grompping in `index_grompp.ndx` (for now this is required, in future version might only be needed if special groups are used in the `mdp`-file). Here is a tree view of how you project folder should look like to begin:

```
├── config.py *
├── index_grompp.ndx
├── start.gro
├── mdrun.mdp
├── topol.top
└── sbatch_launch.sh *
```

The files marked with an asterisk can be copied and modified from templates with `fst make_templates`, while the rest you should provide yourself. All the file names can differ from these and are defined in the config (except the config name, of course), these are just the defaults. The config name can be given as a command line argument **before the command** (e.g. `fst -c myconfig.py choose`).

## Usage

Run `fst -h` for help or `fst <cmd> -h` for help on specific command.

Almost all configuration should be done in the `config.py` file. The tool does have a few command line tools, but in principle these should only affect *what* is being done, not *how* it is being done (e.g. whether data should be pulled from teh remote is a command line argument, but the remote name, excludes and such are in config).

If you have the configuration file in the same directory, named as `config.py`, just run as below. Otherwise add `fst -c <path/to/config>.py <cmd>`.

### Basic commands

To initialize the first epoch and rsync it to the remote

```sh
fst init --push
```

After an epoch has finished, sync down from remote and make the choices for next epoch

```sh
fst choose --pull
```

Check the figure and if everything is fine, sync up to remote

```sh
fst push
```

Then go to the remote machine and run the sbatch script under the epoch folder. Repeat the last two commands and run simulations until satisfied with the results.

### Configuration

At the top of the config file (well, really anywhere in the file, but at the top in the examples) you can setup the remote options, eg. the directory on the remote machine to be synced, and the remote name. The syncing will be done using `rsync`, to `<remote_name>:<remote_path>/`, so remote name should either be the full url (with `<username>@` prepended if username is not the same as on local machine), or just the hostname if you've defined it in `~/.ssh/config`.


The next most important settings to change are the selections, where you should define the selections to be used firstly to calculate the function, and secondly to calculate the clustering. Both can be either an MDAnalysis selection string, or if the `index_file` variable points to a gromacs index file, optionally an index group in the file. Note though, when using atom names in the selection, that in order to get the chain information properly, the initial gro file will be made into a PDB file, which will be used as the structure file when reading trajectory info, and this process sometimes may change the atom naming scheme. The example case README has a bit more info on this.

Lastly you should define the function itself in `function_val`. The only parameter it will get the positions of one repetition at a time as a `ndarray(n,m,3)`, for n frames and a selection of m atoms. It will also be used to calculate the intial function value, so it should work when `n=1` (But the shape will still be (1,m,3), so in most cases any code that works for the general case works for the initial value). As per MDAnalysis default units, the positions will be in ångströms.


### Configuration options

Here is a listing of all the variables you should need in normal usage. The default values are defined in `src/default_config.py`. If there is no default value, this means that it must be defined in the config file. For now there is no check that everything is defined, and the program will simply crash at some point.

#### Remote options

| Variable | Description | Default value |
| --- | - | - |
| `remote_dir` | Remote dir to sync the project to |  - |
| `remote_name` | Hostname of the remote machine | - |
| `rsync_excludes` | A list of patterns to exclude from the syncing process. | *See below* |

By default the rsync excludes are

```python
rsync_excludes = ["fval_data.npz", ".mdrun.xtc_offsets.lock", ".mdrun.xtc_offsets.npz"]
```


#### Sampling variables


| Variable | Description | Default value |
| --- | - | - |
| `N` | The number of repetitions per epoch. Changing it after a few epochs, will only affect new epochs being run. | - |



#### Simulations grompping


| Variable | Description | Default value |
| --- | - | - |
| `gmx` | The name (including the path if it is not in your `PATH`) of the gromacs binary. | `"gmx"` |
| `mdp` | The name of the `.mdp` file. | `"mdrun.mdp"` |
| `topol` | The name of the topology file. | `"topol.top"` |
| `ndx` | The name of the index file used for grompping. | `"index_grompp.ndx"` |
| `restraint_file` | The file to give to the `-r` option of grompp, to use for restraints. `"start"` will use the starting structure of the repetitions, `"initial"` will use the initial structure. Any other (nonempty) string will be interpreted as a path to the file. If it evaluates to False, teh option will not be passed to gromacs. | `False` |
| `maxwarn` | How many warnings to ignore. Can be useful, e.g. when generating velocities while using Nosé-Hoover thermostat. | 0 |



#### Function calculations


| Variable | Description | Default value |
| --- | - | - |
| `select_str` | The selection string (or index group) used for the function calculation. | - |
| `select_str_clust` | The selection string (or index group) used for the clustering. | - |
| `initial_struct` | The initial structure file. Can be any file that MDAnalysis can read as a structure file without needing a trajectory. E.g. TPR files are (currently, as of MDAnalysis v2.2.0) not enough. Gro and pdb are guaranteed to work. | `"start.gro"` |
| `index_file` | `None` or the name of the index file, where the index groups of the selections are. | `None` |
| `minval`/`maxval` | The boundaries for the functional sampling. One of them can be the string "start", to use the initial value, and either can be `None` to use no boundary (only useful at the beginning). | - |
| `function_val` | The function to be sampled. | - |


#### Trajectory transformations


##### General

These options change the on-the-fly transformations that are done to the trajectory as it is being read.

| Variable | Description | Default value |
| --- | - | - |
| `unwrap_mols` | Whether to unwrap molecules (make them whole if broken over the PBC). | `False` |
| `unwrap_sel` | The selection string (or index group) used for unwrapping. Only atoms within this selection will be considered, so if parts of a chain are missing, only those parts that have unbroken bonded graphs will be made whole, each of them separately. In general you should make a selection with at least the backbone connecting whichever parts you want whole. This should be fast enough even with large proteins, but including the water can effect performance badly. The default selection of "protein" should work in most cases.| `"protein"` |
| `unwrap_starters` | `None` or the selection string (or index group) used as "starters". These atoms will be the starting points of making the molecules whole. In effect, they are guaranteed to be inside the box after the process. If `None`, or if a molecule does not have any atoms in the group, the atom with the smallest index will be used. | `None` |
| `mols_in_box` | Whether to put centre of mass of molecules back in box after unwrapping. Only the coordinates in `unwrap_sel` are moved and considered for the centre of mass, **however**, unlike for unwrapping, molecules are moved as a whole, even if they are missing parts in between. This option is ignored if `unwrap_mols=False`.| `False` |

##### Clustering

These options only affect the coordinates used for clustering. The transformations are added after the fval coordinates have been read, but before the clustering coordinates.

| Variable | Description | Default value |
| --- | - | - |
| `clust_centre` | Whether to move the centre of mass of the clustering coordinates to match the inital structure. Ignored if `clust_superpos=True` | `True` |
| `clust_superpos` | Whether to move the centre of mass of the clustering coordinates to match the inital structure **and** rotate for optimal fit (minimum mass weighted RMSD). | `False` |

#### Advanced options

In most cases you should not need these, but may be helpful in others.

| Variable | Description | Default value |
| --- | - | - |
| `ignore_reps` | List of repetitions (as tuples of epoch number and repetiton number, as ints) to ignore. E.g. `[(1,3),(4,5)]` to ignore repetition 3 of epoch 1 and repetition 5 of epoch 4. | `[]` |
| `ignore_epcs` | List of epoch numbers (as ints) to ignore. E.g. `[1,3]` to ignore epochs 1 and 3. | `[]` |
| `ignore_from_start` | Number of frames to ignore from the start. | 0 |
| `stride` | Only consider every N'th frame | 1 |
| `data_per_bin` | On average, at least this much data should be in the histogram bins within the boundaries. | 100 |
| `maxbins` | Maximum number of histogram bins within the boundaries. | 100 |
| `minchoice` | Minimum number of frames, from which a choice will be made. If the chosen bin has less frames, this many of the closest frames in function value are used as the pool of frames to choose from. | 100 |
| `peak_options` | A dictionary of keyword arguments for [scipy.signal.find_peaks](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html), which will be used for dfinding the local minima and maxima from the distributions. | `{"width":5, "distance":10}` |
| `choice_crit` | Maximum height to choose from, as fraction from minimum found height to maximum. The final maximum height to allow choosing from will be `(maxh-minh)*choice_crit+minh`, where maxh and minh are the highest and lowest histogram bin heights. Default is halfway between. | 0.5 |
| `allow_choice_duplicates` | Whether to allow choosing the same frame during the same epoch more than once. If this is False, minchoice **must** be greater than, or equal to N. This only affects duplicates from a single choice, so edge cases of multiple choices within a region of low sampling, the minchoice values might overlap and still produce duplicates. | `False` |
| `cumulative_histogram` | If True, the hstogram will be made as a cumulative one. Otherwise each epoch is plotted only including its own data | `True` |
| `histogram_max_epochs` | Maximum number of (latest) epochs to plot into histogram. For the cumulative case all epoch are plotted, but older ones with single color and no key in legend. | 15 |
| `clust_data_per_bin` | Same as `data_per_bin`, but for the histogram when clustering. | 1000 |
| `clust_maxbins` | Same as `maxbins`, but for the histogram when clustering. | 10 |
| `epochs_pre_clust` | Run this many epochs without clustering. | 3 |
| `maxclust` | Maximum number of clusters in bin. | 15 |
| `clust_choice_frac` | Fraction of the choices to be made from clustering. | 0.5 |
| `clust_tol` | Tolerance to choose the number of clusters. | 0.1 |



### Specific use cases and tips

#### Multiple starting structure

If you want to have multiple different starting points, or even different starting points for each repetition this is now possible. All pdb/gro/xtc files in the same dir as `initial_struct`, matching `<name>*.<suffix>` (where `<name>` is the name of the initial structure without file ending), excluding `initial_struct`, will be used as starting frames.

For example, you can set `initial_struct="initial/start.gro"`, make the folder "initial" and include gro files named `start.gro`,`start01.gro`,`start02.gro`,...,`start16.gro` to have 16 different frames. "start.gro" will be used for making all selections and as astructure file during and after initialization. The rest will only be used during initialization as the starting structures. For the same effect you can also have just `start.gro` and a single xtc file as `start.xtc` with 16 frames.

The frames will be read in alphabetic order, `pdb` files before `gro`, and `gro` files before `xtc`. If there are in total more frames than repetitions, a warning will be produced and the N first frames will be used (for N repetitions). If there are less frames than repetitions, they will be read in order and filled in "round robin" fashion (e.g.for three frames and 8 repetitions, the repetiton starting points would be [1,2,3,1,2,3,1,2]).

**Remember** that you do always need the `initial_struct` to be used as a structure file, and it will only be used as a starting structure, if it is the only one found. If you have three gro files under `initial`: `start.gro`, `start1.gro`, `start2.gro`, this will be counted as two starting frames, the first from `start1.gro` and second from `start2.gro`. If the `initial_struct` is the only file found, then it will be used as the only frame.


#### Using the MDAnalysis selection objects in the calculation of the function

To use the selection object, use the variable `sel` in the function (make sure that you do not overwrite it making it local, or use `global sel` at the beginning of the function). Similarily you can access the universe with `struct`.

As an example consider the external comand `command` that reads in a single frame from "tmp_analysis.pdb" and outputs the result to stdout. The function could then be


```python
import subprocess as subp

def function_val(positions):
    results = []
    # Iterate over the first axis (time)
    for p in positions:
        # Set the selection positions
        sel.positions = p
        # Write the selection
        sel.write("tmp_analysis.pdb")
        # Run external command
        compProc = subp.run("command", capture_output=True, text=True)
        # get output
        results.append(float(compProc.stdout))
    
    return results
```

For now the program only makes a  single selection. However, do remember that the config file is a python script, so you can write any arbitrary code. This means you can make your own selection objects by adding the following lines anywhere in the config:

```python
import MDAnalysis as mda
_univ  = mda.Universe("initial/start.pdb")
_mysel = univ.select_atoms("protein and resid 200-300 and backbone")
```

Here we used a leading underscore, so that we do not accidentally overwrite any variables the program adds to the module, like `sel` and `struct`. Also here we cannot yet use those variables, since they are of course not added yet when the config is imported. Another note about this is that the choose command has to be run, or a PDB produced manually before this works, since the automatic PDB maker only runs after the config is imported. You can also just use the gro file as structure, but then make sure the atom naming is correct and you do not need chain info.

In the future multiple selection might be supprted and this will become less 