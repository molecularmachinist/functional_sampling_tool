# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os,math
from scipy.signal import find_peaks

from . import utils
from . import inout

# Type hints
from typing import Any, Tuple, List
from numpy.typing import NDArray

"""
This module includes a class that holds the loaded data, histograms and some
methods for choosing frames.
"""

class FrameChooser():

    def __init__(self, cfg: Any, fval: NDArray[np.float_], frms: NDArray[np.int_], reps: NDArray[np.int_], epcs: NDArray[np.int_]):
        """
        Takes the config and data as 4 N-length arrays with the rep, fval, epoch
        and frame number (within teh specific epoch and rep) of each datapoint/frame.
        Saves the data and calculates histogram.
        """
        self.cfg  = cfg
        self.fval = fval
        self.frms = frms
        self.reps = reps
        self.epcs = epcs
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        self.u_reps = []
        for e in self.u_epcs:
            self.u_reps.append(np.unique(reps[self.epcs==e]))

        self.nextepoch = self.u_epcs[-1]+1
        # make sure this is not an ignored epoch
        while (self.nextepoch in cfg.ignore_epcs):
            self.nextepoch+=1

        self._make_hist()

    @classmethod
    def fromReadData(cls,cfg: Any, load_fval: bool):
        """
        Factory method to easily load data and make the object
        """
        fval,_,frms,reps,epcs = inout.load_data(cfg,load_fval)
        return cls(cfg, fval, frms, reps, epcs)


    def update(self, fval: NDArray[np.float_], frms: NDArray[np.int_], reps: NDArray[np.int_], epcs: NDArray[np.int_]) -> None:
        """
        Update the new data and recalculates histogram
        """
        self.reps = reps
        self.fval = fval
        self.epcs = epcs
        self.frms = frms
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        self.u_reps = []
        for e in self.u_epcs:
            self.u_reps.append(np.unique(reps[self.epcs==e]))

        self.make_hist()

    def _make_hist_no_cfg(self, maxbins: int, data_per_bin: int, minval: float, maxval: float, hard_boundaries: bool = False) -> None:
        """
        The function that actually makes the histogram. Does not read values from cfg,
        but gets them as parameters. This way different implementations can use
        different values.
        hard_boundaries=True means that bins are calculated strictly between the boundaries only.
        Otherwise the bin edges do not have to got with the boundaries and bins exist outside boundries also.
        """
        print(f"Making histogram with total {len(self.fval)} data")
        # Get extent of data and boundaries
        largest_val = np.max(self.fval)
        lowest_val  = np.min(self.fval)
        self.binmax = min(largest_val, maxval)
        self.binmin = max(lowest_val,  minval)
        # Calculate binsize
        maxbins_calc = np.sum((self.fval>minval)*(self.fval<maxval))//data_per_bin
        maxbins = min(maxbins_calc,maxbins)
        self.binsize = (self.binmax-self.binmin)/maxbins
        print(f"Calculated maxbins {maxbins_calc}, final maxbins {maxbins}")
        print(f"{np.sum((self.fval>minval)*(self.fval<maxval))} data inside boundaries")
        # Make histogram
        if(hard_boundaries):
            self.bin_edges = np.linspace(self.binmin,self.binmax, maxbins)
        else:
            self.bin_edges = np.arange(lowest_val,largest_val+self.binsize , self.binsize)
        self.hist, _ = np.histogram(self.fval, bins=self.bin_edges)
        self.bin_centers = (self.bin_edges[:-1]+self.bin_edges[1:])/2
        # Make a mask of which bins have higher edge larger than minval
        # AND lower edge lower than maxval
        self.hist_mask = (self.bin_edges[1:]>minval)*(self.bin_edges[:-1]<maxval)


    def _make_hist(self) -> None:
        self._make_hist_no_cfg(
                self.cfg.maxbins,
                self.cfg.data_per_bin,
                self.cfg.minval,
                self.cfg.maxval
            )


    def make_choices(self, prechoices: int = 0, plot: bool = True) -> Tuple[NDArray[np.float_],NDArray[np.int_],NDArray[np.int_],NDArray[np.int_]]:
        """
        Uses the histogram to choose bins and returns the bin indices of the choices.
        prechoices is the number of choices already done.
        """
        smoothed = utils.rolling_mean(self.hist)
        nanmask = np.isfinite(smoothed)*self.hist_mask
        maxims,max_crit = find_peaks(smoothed[nanmask], **self.cfg.peak_options)
        minims,min_crit = find_peaks(-smoothed[nanmask], **self.cfg.peak_options)


        maxh = np.max(self.hist[self.hist_mask])
        minh = np.min(self.hist[self.hist_mask])
        # If not ends have sampled been, make sure min height set to zero is
        if(np.max(self.fval)<self.cfg.maxval or np.min(self.fval)>self.cfg.minval):
            minh=0

        # The criteria for choices is half way between max and min heights
        crith = (maxh-minh)/2+minh

        choices = []

        # convert from masked indices to unmasked
        indexes = np.arange(len(self.hist))[nanmask]
        for p in minims:
            if(self.hist[nanmask][p]<crith):
                choices.append(indexes[p])

        zero_mask = self.hist_mask*(self.hist!=0)
        indexes = np.arange(len(self.hist))[zero_mask]

        # TODO: refactor below
        #Handle not peaks-edge case
        if(maxims.size==0):
            if(self.hist[self.hist_mask][-1]<crith):
                choices.append(indexes[-1])
            if(self.hist[self.hist_mask][0]<crith):
                choices.append(indexes[0])
        else:
            # Check if first/last extrema is a maxima -> also include the far end(s)
            if((minims.size==0 or maxims[-1]>minims[-1]) and self.hist[self.hist_mask][-1]<crith):
                choices.append(indexes[-1])
            if((minims.size==0 or maxims[0]<minims[0]) and self.hist[self.hist_mask][0]<crith):
                choices.append(indexes[0])

        weights = [(crith-self.hist[c]) for c in choices]
        weights /= np.sum(weights)
        srt_ind = np.argsort(-weights)
        #sort choices
        choices = [choices[i] for i in srt_ind]

        # Each point has been added once
        len_choice = len(choices)
        # Add more depending on weight. Multiplication is floored, so between 0 and len_choice-1 too few are added
        for i in srt_ind:
            for j in range(math.floor(weights[i]*(self.cfg.N-len_choice-prechoices))):
                choices.append(choices[i])

        # Fill the rest in sorted order
        for i in range(self.cfg.N-len(choices)-prechoices):
            choices.append(choices[srt_ind[i]])

        # Make sure too many were not added
        # This would either be because more than needed minima found, or the one
        # bug, I have yet to find
        choices = choices[:(self.cfg.N-prechoices)]
        print(f"Chose {len(choices)} bins, {len(np.unique(choices))} unique")

        if(plot):
            self._plot_choices(crith,smoothed,choices,nanmask,maxims,minims)
        return self.choose_frames(choices)


    def choose_frames(self, chosen_bins: List[int]) -> Tuple[NDArray[np.float_],NDArray[np.int_],NDArray[np.int_],NDArray[np.int_]]:
        """ Input parameters:
                - chosen_bins : Length N list of the indices of the bins that have been chosen.
                                Duplicates (starting from the same bin) are simply many times in the list.
            Output:
                - v : Length N list of the fvals for each new rep
                - e : Length N list of the epochs each new rep comes from
                - r : Length N list of the rep within the epoch each new rep comes from
                - f : Length N list of the frame within the rep each new rep comes from

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
        for bi,ci in zip(*np.unique(chosen_bins, return_counts=True)):
            if(bi<len(self.bin_edges)-2):
                vals_in_bin = (self.fval >= self.bin_edges[bi])*(self.fval < self.bin_edges[bi+1])
            else:
                vals_in_bin = (self.fval >= self.bin_edges[bi])*(self.fval <= self.bin_edges[bi+1])
            
            # How many are in bin
            n_in_bin = np.sum(vals_in_bin)
            # Check if enough are to satisfy minchoice criteria
            if(n_in_bin<self.cfg.minchoice):
                # If not, we choose minchoice closest to bin center
                n_in_bin = self.cfg.minchoice
                cnt = (self.bin_edges[bi]+self.bin_edges[bi+1])/2
                # Sort by distance to cnt and choose the n_in_bin first ones
                vals_in_bin = np.argsort(np.abs(self.fval-cnt))[:n_in_bin]

            # Choose ci from bin
            ndx = rng.choice(n_in_bin,size=ci,replace=self.cfg.allow_choice_duplicates)

            v.extend(self.fval[vals_in_bin][ndx])
            e.extend(self.epcs[vals_in_bin][ndx])
            r.extend(self.reps[vals_in_bin][ndx])
            f.extend(self.frms[vals_in_bin][ndx])

        return v,e,r,f


    def print_choices(self,val: NDArray[np.float_], epc: NDArray[np.int_], rep: NDArray[np.int_], frm: NDArray[np.int_]) -> None:
        print("Final choices:")
        for v,e,r,f in zip(val, epc, rep, frm):
            print(f"frm {f}, rep {r}, epc {e}, fval={v}")


    def plot_hist(self) -> None:
        """
        Make a histogram of the current total data, as well as maxepochs latest
        OR cumulative histograms after each epoch, with labels in legend only for
        maxepoch latest.
        """
        maxepochs = -self.cfg.histogram_max_epochs
        if(maxepochs==1):
            maxepochs = 0
        #if doing cumulative, also plot the epochs previous to maxepochs, but without labels
        if(self.cfg.cumulative_histogram):
            for i in self.u_epcs[:maxepochs]:
                _ep_hist, _ = np.histogram(self.fval[self.epcs<=i], bins=self.bin_edges)
                plt.plot(self.bin_centers, _ep_hist, "--", color="C0", alpha=0.5)

        # Plot maxepochs latest epochs
        for i in self.u_epcs[maxepochs:]:
            # The mask is <= for cumulative and == for non cumulative
            _mask = self.epcs<=i if (self.cfg.cumulative_histogram) else self.epcs==i
            # histogram of masked values
            _ep_hist, _ = np.histogram(self.fval[_mask], bins=self.bin_edges)
            # Plot with label
            plt.plot(self.bin_centers, _ep_hist, "--", color="C%d"%(i+2), label="Epoch %d"%i)

        # Total histogram, with alpha=0.5 to make it lighter
        plt.plot(self.bin_centers, self.hist, color="C0", alpha=0.5)
        # Total histogram only within teh boundaries without alpha
        plt.plot(self.bin_centers[self.hist_mask], self.hist[self.hist_mask], color="C0", label="Total")
        # Show starting value and boundaries
        plt.axvline(self.cfg.startval, linestyle="-.", color="C1", label="Start", alpha=0.5)
        plt.axvline(self.cfg.minval, linestyle="-.", color="C2", label="Boundaries")
        plt.axvline(self.cfg.maxval, linestyle="-.", color="C2")
        # legend, on teh right side outside of plot
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        # Make directory and save fig
        os.makedirs(self.cfg.fig_output_dir / ("epoch%02d"%self.u_epcs[-1]),exist_ok=True)
        plt.savefig(self.cfg.fig_output_dir / ("epoch%02d"%self.u_epcs[-1]) / "hist.png")
        plt.clf()


    def _plot_choices(self,crith,smoothed,choices,nanmask,maxims,minims) -> None:
        plt.axhline(crith, linestyle="--", c="r", alpha=0.5, label="Max height for choices")
        plt.plot(self.bin_centers, self.hist, color="C0", alpha=0.5)
        plt.plot(self.bin_centers[self.hist_mask], self.hist[self.hist_mask], color="C0", label="data")
        plt.plot(self.bin_centers, smoothed, color="C1", label="Smoothed")
        plt.plot(self.bin_centers[nanmask][minims], smoothed[nanmask][minims], "rv", label="Minima")
        plt.plot(self.bin_centers[nanmask][maxims], smoothed[nanmask][maxims], "bv", label="Maxima")
        plt.plot(self.bin_centers[choices], self.hist[choices], "g^", label="Choices")
        plt.legend()
        os.makedirs(self.cfg.fig_output_dir / ("epoch%02d"%self.u_epcs[-1]),exist_ok=True)
        plt.savefig(self.cfg.fig_output_dir / ("epoch%02d"%self.u_epcs[-1]) / "choices.png")
        plt.clf()
