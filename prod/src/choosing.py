# -*- coding: utf-8 -*-
import numpy as np

from . import config as cfg

"""
This module includes a class that holds the loaded data, histograms and some
methods for choosing frames.
"""

class FrameChooser():

    def __init__(self, reps, fval, epcs, frms):
        self.reps = reps
        self.fval = fval
        self.epcs = epcs
        self.frms = frms
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        for e in self.u_epcs:
            self.u_reps = np.unique(reps[self.epcs==e])

    def make_hist(self):
        print(f"Making histogram with total {len(self.fval)} data")
        # Get extent of data and boundaries
        largest_val = np.max(self.fval)
        lowest_val  = np.min(self.fval)
        self.binmax = min(largest_val, cfg.maxval)
        self.binmin = max(lowest_val,  cfg.minval)
        # Calculate binsize
        maxbins = np.sum((self.fval>cfg.minval)*(self.fval<cfg.maxval))//config.data_per_bin
        maxbins = min(maxbins,config.maxbins)
        self.binsize = (self.largest_val-self.lowest_val)/maxbins
        print(f"Calculated maxbins {maxbins}, final maxbins {maxbins}")
        print(f"data inside boundaries {np.sum((self.fval>cfg.minval)*(self.fval<cfg.maxval))}")
        # Make histogram
        self.bin_edges = np.arange(lowest_val,largest_val+self.binsize , self.binsize)
        self.hist, _ = np.histogram(self.fval, bins=self.bin_edges)
        self.bin_centers = (self.bin_edges[:-1]+self.bin_edges[1:])/2
        # Make a mask of which bins have higher edge larger than minval
        # AND lower edge lower than maxval
        self.hist_mask = (self.bin_edges[1:]>cfg.minval)*(self.bin_edges[:-1]<cfg.maxval)


    def plot_hist():

        for i in self.u_epcs:
            _ep_hist, _ = np.histogram(fval[epcs==i], bins=bin_edges)
            plt.plot(bin_centers, _ep_hist, "--", color="C%d"%(i+2), label="Epoch %d"%i)
        plt.plot(self.bin_centers, self.hist, color="C0", alpha=0.5)
        plt.plot(self.bin_centers[self.hist_mask], self.hist[self.hist_mask], color="C0", label="Total")
        plt.axvline(startval, linestyle="-.", color="C1", label="Start", alpha=0.5)
        plt.axvline(minval, linestyle="-.", color="C2", label="Boundaries")
        plt.axvline(maxval, linestyle="-.", color="C2")
        plt.legend()
        plt.show()
