# -*- coding: utf-8 -*-
import numpy as np
import os

from . import utils

"""
This module includes a class that holds the loaded data, histograms and some
methods for choosing frames.
"""

class FrameChooser():

    def __init__(self, cfg, reps, fval, epcs, frms):
        """
        Takes the config and data as 4 N-length arrays with the rep, fval, epoch
        and frame number (within teh specific epoch and rep) of each datapoint/frame.
        Saves the data and calculates histogram.
        """
        self.cfg  = cfg
        self.reps = reps
        self.fval = fval
        self.epcs = epcs
        self.frms = frms
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        for e in self.u_epcs:
            self.u_reps = np.unique(reps[self.epcs==e])

        self.make_hist()


    def update(self, reps, fval, epcs, frms):
        """
        Update the new data and recalculates histogram
        """
        self.reps = reps
        self.fval = fval
        self.epcs = epcs
        self.frms = frms
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        for e in self.u_epcs:
            self.u_reps = np.unique(reps[self.epcs==e])

        self.make_hist()


    def make_hist(self):
        print(f"Making histogram with total {len(self.fval)} data")
        # Get extent of data and boundaries
        largest_val = np.max(self.fval)
        lowest_val  = np.min(self.fval)
        self.binmax = min(largest_val, self.cfg.maxval)
        self.binmin = max(lowest_val,  self.cfg.minval)
        # Calculate binsize
        maxbins = np.sum((self.fval>self.cfg.minval)*(self.fval<self.cfg.maxval))//config.data_per_bin
        maxbins = min(maxbins,config.maxbins)
        self.binsize = (self.largest_val-self.lowest_val)/maxbins
        print(f"Calculated maxbins {maxbins}, final maxbins {maxbins}")
        print(f"data inside boundaries {np.sum((self.fval>self.cfg.minval)*(self.fval<self.cfg.maxval))}")
        # Make histogram
        self.bin_edges = np.arange(lowest_val,largest_val+self.binsize , self.binsize)
        self.hist, _ = np.histogram(self.fval, bins=self.bin_edges)
        self.bin_centers = (self.bin_edges[:-1]+self.bin_edges[1:])/2
        # Make a mask of which bins have higher edge larger than minval
        # AND lower edge lower than maxval
        self.hist_mask = (self.bin_edges[1:]>self.cfg.minval)*(self.bin_edges[:-1]<self.cfg.maxval)


    def make_choices(self):
        """
        Uses the histogram to choose bins and returns the bin indices of the choices.
        """
        smoothed = utils.rolling_mean(-hist)
        nanmask = np.isfinite(smoothed)*self.hist_mask
        maxims,max_crit = scipy.signal.find_peaks(smoothed[nanmask], width=5, distance=10)
        minims,min_crit = scipy.signal.find_peaks(-smoothed[nanmask], width=5, distance=10)


        maxh = np.max(self.hist[self.hist_mask])
        minh = np.min(self.hist[self.hist_mask])
        # If not ends have sampled been, make sure min height set to zero is
        if(largest_val<maxval or lowest_val>minval):
            minh=0

        # The criteria for choices is half way between max and min heights
        crith = (maxh-minh)/2+minh

        choices = []
        # convert from masked indices to unmasked
        indexes = np.arange(len(hist))[nanmask]
        for p in peaks:
            if(hist[nanmask][p]<crith):
                choices.append(indexes[p])

        zero_mask = mask*(hist!=0)
        indexes = np.arange(len(hist))[zero_mask]

        # TODO: refactor below
        #Handle not peaks-edge case
        if(maxims.size==0):
            if(self.hist[self.hist_mask][-1]<crith):
                choices.append(indexes[-1])
            if(self.hist[self.hist_mask][0]<crith):
                choices.append(indexes[0])
        else:
            # Check if first/last extrema is a maxima -> also include the far end(s)
            if(minims.size>0 and maxims[-1]<minims[-1] and self.hist[self.hist_mask][-1]<crith):
                choices.append(indexes[-1])
            if(minims.size>0 and maxims[0]>minims[0] and self.hist[self.hist_mask][0]<crith):
                choices.append(indexes[0])

        weights = [(crith-hist[c]) for c in choices]
        weights /= np.sum(weights)
        print(choices,(weights*(N-len(choices))))
        srt_ind = np.argsort(-weights)
        print([choices[i] for i in srt_ind],(weights[srt_ind]*(N-len(choices))))

        # Each point has been added once
        len_choice = len(choices)
        # Add more depending on weight. Multiplication is floored, so between 0 and len_choice-1 too few are added
        for i in srt_ind:
            for j in range(floor(weights[i]*(self.cfg.N-len_choice))):
                choices.append(choices[i])

        # Fill the rest in sorted order
        for i in range(N-len(choices)):
            choices.append(choices[srt_ind[i]])


        print(f"Chose {len(choices)} bins, {len(np.unique(choices))} unique")

        self.plot_choices(crith,smoothed,choices,nanmask,maxims,minims)
        return choices

    def choose_frames(self, chosen_bins):
        """ Input parameters:
                - chosen_bins : Length N list of the indices of the bins that have been chosen.
                                Duplicates (starting from the same bin) are simply many times in the list.
            Output:
                - v : Length N array of the fvals for each new rep
                - e : Length N array of the epochs each new rep comes from
                - r : Length N array of the rep within the epoch each new rep comes from
                - f : Length N array of the frame within the rep each new rep comes from

            Other "inputs":
                Saved already in the object:
                - bin_edges : Edges of the bins, so bin_edges[i] is the left, and bin_edges[i+1]
                              the right edge of the i'th bin. As per numpy.histogram, each bin is half open
                              [left_edge, right_edge), except for the last one [left_edge, right_edge].
                - fval      : An array containing all the fvals from all the simulations so far.
                - epcs      : An array of same shape as fval, with each element being the corresponding
                              epoch each value comes from.
                - reps      : Same as epcs, but with the corresponding rep.
                - frms      : Same as epcs and reps, but with the information of the frame.
        """


        rng = np.random.default_rng()
        # lists for value, rpoch, rep and frame
        v,e,r,f=[],[],[],[]
        for bi in chosen_bins:
            vals_in_bin = (self.fval >= self.bin_edges[bi])*(self.fval < self.bin_edges[bi+1])
            ndx = rng.choice(np.sum(vals_in_bin))

            v.append(self.fval[vals_in_bin][ndx])
            e.append(self.epcs[vals_in_bin][ndx])
            r.append(self.reps[vals_in_bin][ndx])
            f.append(self.frms[vals_in_bin][ndx])

        return np.array(v), np.array(e), np.array(r), np.array(f)

    def plot_hist(self):

        for i in self.u_epcs:
            _ep_hist, _ = np.histogram(fval[epcs==i], bins=bin_edges)
            plt.plot(bin_centers, _ep_hist, "--", color="C%d"%(i+2), label="Epoch %d"%i)
        plt.plot(self.bin_centers, self.hist, color="C0", alpha=0.5)
        plt.plot(self.bin_centers[self.hist_mask], self.hist[self.hist_mask], color="C0", label="Total")
        plt.axvline(startval, linestyle="-.", color="C1", label="Start", alpha=0.5)
        plt.axvline(minval, linestyle="-.", color="C2", label="Boundaries")
        plt.axvline(maxval, linestyle="-.", color="C2")
        plt.legend()
        os.makedirs("figs/epoch%02d"%(self.u_epcs[-1]),exist_ok=True)
        plt.savefig("figs/epoch%02d/hist.png"%(self.u_epcs[-1]))
        plt.clf()


    def plot_choices(self,crith,smoothed,choices,nanmask,maxims,minims):
        plt.axhline(crith, linestyle="--", c="r", alpha=0.5, label="Max height for choices")
        plt.plot(self.bin_centers, hist, color="C0", alpha=0.5)
        plt.plot(self.bin_centers[self.hist_mask], self.hist[self.hist_mask], color="C0", label="data")
        plt.plot(self.bin_centers, -smoothed, color="C1", label="Smoothed")
        plt.plot(self.bin_centers[nanmask][minims], -smoothed[nanmask][maxims], "rv", label="Minima")
        plt.plot(self.bin_centers[nanmask][maxims], -smoothed[nanmask][minims], "bv", label="Maxima")
        plt.plot(self.bin_centers[choices], self.hist[choices], "g^", label="Choices")
        plt.legend()
        os.makedirs("figs/epoch%02d"%(self.u_epcs[-1]),exist_ok=True)
        plt.savefig("figs/epoch%02d/choices.png"%(self.u_epcs[-1]))
        plt.clf()
