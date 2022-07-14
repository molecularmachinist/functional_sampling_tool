#!/usr/bin/env bash

rm -rf initialize/epoch01 initialize/epoch02 initialize/figs after_e1/epoch02 after_e1/figs remote_dir after_e1/analysis
rm -rf after_e1/epoch01/rep*/fval_data.npz  after_e1/epoch01/rep*/.mdrun.xtc_offsets.npz after_e1/epoch01/rep*/.mdrun.xtc_offsets.lock
rm -f after_e1/initial/start.pdb after_e1/initial/output_trjconv.txt
