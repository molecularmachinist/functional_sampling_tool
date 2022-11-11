# Functional Sampling Tool - Usage Example

Here we take a very small protein (α-Conotoxin MII, PDB id [1MII](https://www.rcsb.org/structure/1MII)) to have a lightweight usage example. As a function we will choose the distance between the GLY1 and ASN14 alpha carbons. To my knowledge this function does not bear any biological purpose, so this tutorial should only be thought about in terms of using the tool. Coming up with the function to study will be all up to you when you start running this on your own.

## Setup

To begin, you should have downloaded this repo, and installed all the required python packages, as well as have gromacs installed. To start off, `cd` into the `initialize` folder and take a look around. This is how you should set up your own project, when you start running one. You'll notice the `config.py` is already set up, a starting structure is in `initial/start.gro`, along with a corresponding topology in `topol.top`.

The `sbatch_launch.sh` is where you would normally add you SLURM commands to run the job on a supercomputer. In this example case it is just a bash script that runs the repetitions on the local machine. To produce an example of the `sbatch_launch.sh` and `config.py`, you can run `fst make-templates`. Don't run it now, or you would overwrite the current files.

`index_grompp.ndx` is, for now, needed even if you have no special groups in the mdp file. To make one with just the defaults, use
```
echo q | gmx make_ndx -f initial/start.gro -o index_grompp.ndx
```

As in the root README, the commands here assume the tool to be in your PATH and can be run as `fst`. If your gromacs installation is not in your path as `gmx`, you will need to edit the config with

```
gmx = "/path/to/gmx"
```


## System initialization

The first step in the project is of course to make the first epoch. The command is simply

```
fst init
```

This makes a folder `epoch01`, copies the necessary files for each repetition into `epoch/repNN` and grompps the systems. It will also copy the `sbatch_launch.sh` into `epoch01`, with the epoch number filled in on the job name.

In this case the command should actually fail, complaining about nonzero return code from grompp. Check the output file in the error message (e.g. `epoch01/rep01/output_grompp.txt`) to see what the error is caused by. It should be caused by a warning about generating velocities when using Parrinello-Rahman pressure coupling. In this case we know that we do want to do this. The system is already equilbrated, so the pressure coupling should work fine, and we do need to regenerate velocities for each repetition. As such this warning can be ignored. Do make sure that this was what caused the error, and never ignore a warning you are not sure should be ignored.

To ignore the warning we add the below line anywhere in the config

```
maxwarn = 1
```
Then, remove clean the failed command, run 
```
fst clean
```
It makes sure no simulation data is in the latest (highest number) epoch, and then removes it. Now you can rerun the init command and everything should work fine. If something is still off, rerun the clean command, fix your system and rerun the init command.


At this point, you would usually push these to a supercomputer where you actually run the simulations. You can try to do the push command, and (if you can ssh onto localhost) it will copy these to a local folder (`remote_folder` inside the `example` folder). The command to do it is

```
fst push
```


## Running the epoch

Normally, you would now ssh into the remote machine, cd into the `epoch01` folder and run `sbatch sbatch_launch.sh` (assuming the computer uss SLURM). Then you'd just wait for the simulations to finish.

If you want, you can run the simulations. `cd` into `epoch01` under the "remote" folder and and run the `sbatch_launch.sh` as a bash script. On my desktop this took about 2h 40min.

You can also skip this and change into the `after_e1` folder, where the simulated trajectories are already present.



## Choosing frames for next epoch

Now `cd` into the other folder, `after_e1`, or use `fst pull` in `initialize` to pull the simulated files from the "remote" directory.

We start by using the main command of the program:

```
fst choose
```

This will read the data from the xtc files, calculate the function values, make the choices and grompp the new repetitions. It will produce, just as before, a folder `epoch02` with everything the folder `epoch01` had after initialization, as well as a file `origin.txt`, which records where the frame was started from. The command also produces a folder `figs`, with a subfolder named after the epoch the frames were made after (in this case `epoch01`). **Always check these figures**, to see what choices the tool made! If you are not satisfied, use again the `clean` command, play with the options in `config.py` and rerun the `choose` command. In many cases simply changing the boundaries is enough.

To run the choosing, without grompping (so that it is faster and you do not need to remove the folder each time), use

```
fst choose --choose-only
```

After you are satisfied in these results, push the new epoch to the remote just as before.


## Cleaning up after ourselves

If you want to redo the example, or simply keep the folder tidy, `cd` back to the `example/` folder and run
```
./clean_example
```
Now any files the example run produced will have been deleted.
