# Functional Sampling Tool

## Requirements

1. Python 3 with NumPy and matplotlib
1. [MDAnalysis](https://docs.mdanalysis.org/stable/index.html)
1. [Numba](https://numba.pydata.org/)
1. [scikit-learn](https://scikit-learn.org/stable/)



If you have conda, you can make sure all dependencies are met by running

```
conda install -c conda-forge numpy matplotlib numba mdanalysis scikit-learn
```

The tool has been tested with python 3.7 and 3.9. Please let us know if you are using other versions, whether everything works (and especially if it doesn't).

## How it works

### Vocabulary

 1. **Epoch**: the set of N runs that are run at the same time.
 1. **Repetition** (rep): The individual runs within each epoch. i.e. there are N reps in each epoch.

### Epoch choosing

After each epoch, the user defined function is calculated for each frame, and a histogram is calculated, to get the distribution. The tool then smooths the distribution and finds the local minima within the boundaries defined by user. The starting points for each rep in the new epoch will chosen from these minima.

### Clustering

After a set amount of epochs, the tool will start to use a clustering method on the trajectories. This way it will try to identify different states of the structure, that have the same function value. A part of the repetitions will be started from the clusters with the lowest amount of data. The clustering is done by first making a histogram with much larger bins than before, and clustering the data in each bin separately.


## Setup

At minimum you will need to copy the templates folder from `rsc/`, and the `config.py` file to your project folder. You should then modify both the config file and `sbatch_launch.sh` template to fit your particular use case.

You will also need to have a starting structure, which should be put into `initil/start.gro`. You will also need a topology in `topol.top`, mdrun options in `mdrun.mdp` and an index file for grompping in `index_grompp.ndx` (for now this is required, in future version might only be needed if special groups are used in the `mdp`-file). Here is a tree view of how you project folder should look like to begin:

```
├── config.py *
├── index_grompp.ndx
├── initial
│   └── start.gro
├── mdrun.mdp
├── templates
│   └── sbatch_launch.sh *
└── topol.top
```

The files marked with an asterisk can be copied and modified, while the rest you should provide yourself. The names of the index, topology and mdp files can differ from these and are defined in the config. The config name can also differ and can be given as a command line argument.

## Usage

Assuming the tool is in your path, run `fst -h` for help or `fst <cmd> -h` for help on specific command. If it is not in your path run the rest of the commands as `/path/to/fst` instead of just `fst`.

Almost all configuration should be done in the `config.py` file. The tool does have a few command line tools, but in principle these should only affect what is being done, not how it is being done.

If you have the configuration file in the same directory, as `config.py`, just run as below. Otherwise add `fst -c <path/to/config>.py <cmd>`.

### Basic commands

To initialize the first epoch and rsync it to the remote

```
fst init --push
```

After an epoch has finished, sync down from remote and make the choices for next epoch

```
fst choose --pull
```

Check the figure and if everything is fine, sync up to remote

```
fst push
```

Then go to the remote machine and run the sbatch script under the epoch folder. Repeat the last two commands and run simulations until satisfied with the results.

### Configuration

At the top of the config file (well, really anywhere in the file, but at the top in the examples) you can setup the remote options, eg. the directory on the remote machine to be synced, and the remote name. The syncing will be done using `rsync`, to `<remote_name>:<remote_path>/`, so remote name should either be the full url (with `<username>@` prepended if username is not the same as on local machine), or just the hostname if you've defined it in `~/.ssh/config`.


The next most important settings to change are the selections, where you should define the selections to be used firstly to calculate the function, and secondly to calculate the clustering. Both can be either an MDAnalysis selection string, or if the `index_file` variable points to a gromacs index file, optionally an index group in the file. Note though, when using atom names in the selection, that in order to get the chain information properly, the initial gro file will be made into a PDB file, which will be used as the structure file when reading trajectory info, and this process sometimes may change the atom naming scheme.

Lastly you should define the function itself in `function_val`. As the only parameter it will get the positions of one repetition at a time as a `ndarray(n,m,3)`, for n frames and a selection of m atoms. As per MDAnalysis default units, the positions will be in ångströms.


### Configuration options

Here is a listing of all the variables you should need in normal usage. The default values are defined in `src/default_config.py`. If there is no default value, this means that it must be defined in the config file. For now there is no check that everything is defined, and the program will simply crash at some point.

#### Remote options

| Variable | Description | Default value |
| --- | - | - |
| `remote_dir` | Remote dir to sync the project to |  - |
| `remote_name` | Hostname of the remote machine | - |
| `rsync_excludes` | A list of patterns to exclude from the syncing process. | *See below* |

By default the rsync excludes are

```
rsync_excludes = ["config.py", "initial", "templates", "fval_data.npz", ".*.xtc_offsets.npz", "figs", "fst", "src"]
```

#### Sbatch template variables

These variables are copied into the `sbatch_launch.sh` script when it is copied for each epoch.

| Variable | Description | Default value |
| --- | - | - |
| `email` | Where to send an email after the job finishes on the remote | - |
| `account` | The account to be billed (project on CSC). | - |


#### Sampling variables


| Variable | Description | Default value |
| --- | - | - |
| `N` | The number of repetitions per epoch. Changing it asfter a few epochs, will only affect new epochs being run. | - |



#### Simulations grompping


| Variable | Description | Default value |
| --- | - | - |
| `gmx` | The name (including the path if it is not in your `PATH`) of the gromacs binary. | `"gmx"` |
| `mdp` | The name of the `.mdp` file. | `"mdrun.mdp"` |
| `topol` | The name of the topology file. | `"topol.top"` |
| `ndx` | The name of the index file used for grompping. | `"index_grompp.ndx"` |
| `maxwarn` | How many warnings to ignore. Can be useful, e.g. when generating velocities while using Nosé-Hoover thermostat. | 0 |



#### Function calculations


| Variable | Description | Default value |
| --- | - | - |
| `select_str` | The selection string (or index group) used for the function calculation. | - |
| `select_str_clust` | The selection string (or index group) used for the clustering. | - |
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
| `mols_in_box` | Whether to put centre of mas of molecules back in box after unwrapping. Only the coordinates in `unwrap_sel` are moved and considered for the centre of mass, **however**, unlike for unwrapping, molecules are moved as a whole, even if they are missing parts in between. This option is ignored if `unwrap_mols=False`.| |

##### Clustering

These options only affect the coordinates used for clustering. The transformations are added after the fval coordinates are read, but before the clustering coordinates.

| Variable | Description | Default value |
| --- | - | - |
| `clust_centre` | Whether to move the centre of mass of the clustering coordinates to match the inital structure. Ignored if `clust_superpos=True` | `True` |
| `clust_superpos` | Whether to move the centre of mass of the clustering coordinates to match the inital structure **and** rotate for optimal fit (minimum mass weighted RMSD). | `False` |

#### Advanced options

In most cases you should not need these, but may be helpful in others.

| Variable | Description | Default value |
| --- | - | - |
| `maxwarn_add` | In case the atom naming has changed between the topology and the structure (eg when making it into a PDB file), one more warning might be issued by grompp. In that case setting this to `True` will add one more to the number of warnings ignored, in all but the first epoch. | `False` |
| `data_per_bin` | On average, at least this much data should be in the histogram bins within the boundaries. | 100 |
| `maxbins` | Maximum number of histogram bins within the boundaries. | 100 |
| `minchoice` | Minimum number of frames, from which a choice will be made. If the chosen bin has less frames, this many of the closest frames in function value are used as the pool of frames to choose from. | 100 |
| `peak_options` | A dictionary of keyword arguments for [scipy.signal.find_peaks](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html), which will be used for dfinding the local minima and maxima from the distributions. | `{"width":5, "distance":10}` |
| `allow_choice_duplicates` | Whether to allow choosing the same frame during the same epoch more than once. If this is False, minchoice **must** be greater than, or equal to N. This only affects duplicates from a single choice, so edge cases of multiple choices within a region of low sampling, the minchoice values might overlap and still produce duplicates. | `False` |
| `clust_data_per_bin` | Same as `data_per_bin`, but for the histogram when clustering. | 1000 |
| `clust_maxbins` | Same as `maxbins`, but for the histogram when clustering. | 10 |
| `epochs_pre_clust` | Run this many epochs without clustering. | 3 |
| `maxclust` | Maximum number of clusters in bin. | 15 |
| `clust_choice_frac` | Fraction of the choices to be made from clustering. | 0.5 |
| `clust_tol` | Tolerance to choose the number of clusters. | 0.1 |
