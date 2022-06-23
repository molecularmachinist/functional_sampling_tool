# Functional Sampling Tool - Usage Example

## Setup

To begin, you should have downloaded this repo, and installed all the required python packages. To start off, `cd` into the `initialize` folder and take a look around. This is how you should set up your own project, when you start running one. You'll notice the `config.py` is already set up, a starting structure is in `initial/start.gro`, along with a corresponding topology in `topol.top`.

As in the root README, the commands here assume the tool to be in your PATH and can be run as `fst`. If this example folder is in the repo still, you also use `../../fst` instead.


## System initialization

The first step in the project is of course to make the first epoch. The command is simply

```
fst init
```

This makes a folder `epoch01`, copies the necessary files for each repetition into `epoch/repNN` and grompps the the systems. It will also copy the `sbatch_launch.sh` into `epoch01`.

At this point, you would usually push these to a supercomputer where you actually run the simulations. You can try to do the push command, and (if you can ssh onto localhost) it will copy these to a subfolder. The command would be

```
fst push
```

Now we won't actually run the simulations, and will simply move to the other folder in this example.

If, for whatever reason, you want to remake a specific repetition, just remove the folder and rerun the previous command. It will simply ignore those repetitions that already have an existing folder. To rerun everything, remove the epoch folder.


## Choosing frames for next epoch

Now `cd` into the other folder, `after_e1`. It is just as the other folder was after initialization, but also has simulated trajectories. To not make this repo any more bloated, the trajectories have very little frames and the `edr`, `log`, and `cpt` files have been removed. As these do not affect the behaviour and usage of the tool, no problem there. The resulting figures will look a bit horrific due to the low number of samples.

Normally to begin, we'd first pull the changes with

```
fst pull
```

but that is not necessary now. We then use the main command of the program.

```
fst choose
```

This will read the data from the xtc files, calculate the function values, make the choices and grompp the new repetitions. It will produce, just as before a folder `epoch02` with everything the folder `epoch01` had after initialization, as well as a file `origin.txt`, which records where the frame was started from. The command also produces a folder `figs`, with a subfolder named after the epoch the frames were made after (in this case `epoch01`). **Always check these figures**, to see what choices the tool made! If you are not satisfied, play with the option in `config.py` and rerun the command. You will have to remove the newly generated epoch folder first. To run the choosing, without grompping (so that it is faster and you do not need to remove the folder each time), use

```
fst choose --choose_only
```

After you are satisfied in these results (for the sake of this example, ignoring the minuscule amount of data the first epoch has generated), push the new epoch to the remote just as before.


## Cleaning up after ourselves

If you want to redo the example, or simply keep the folder tidy, `cd` back to the `example/` folder and run
```
./clean_example
```
Now any files the example run produced will have been deleted.
