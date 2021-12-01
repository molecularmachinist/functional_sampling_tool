# Functional Sampling Tool

## Requirements

1. Python 3 with NumPy and matplotlib
1. [MDAnalysis](https://docs.mdanalysis.org/stable/index.html)


If you have conda, you can install MDAnalysis with

```
conda install -c conda-forge mdanalysis
```

## How it works

## Setup

## Usage

Assuming the project is in your path, run `fst -h` for help or `fst <cmd> -h` for help on specific command.

If you have the config file in the same directory, as `config.py`, just run as below. Otherwise add `fst -c <path/to/config>.py <cmd>`.

### Basics

To initialize the first epoch and rsync it to the remote

```
fst init --push
```

After an epoch has finished, sync down from remote and make the choices for nect epoch

```
fst choose --pull
```

Check teh figure and if everything is fine, sync up to remote

```
fst push
```
