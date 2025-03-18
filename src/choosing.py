# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.signal import find_peaks

from . import utils
from . import inout

from .exceptions import NotEnoughDataError

# Type hints
from typing import Any, Tuple, List
from numpy.typing import NDArray

"""
This module includes a class that holds the loaded data, histograms and some
methods for choosing frames.
"""


class FrameChooser():

    def __init__(self,
                 cfg: Any,
                 fval: NDArray[np.float_],
                 frms: NDArray[np.int_],
                 reps: NDArray[np.int_],
                 epcs: NDArray[np.int_]):
        """
        Takes the config and data as 4 N-length arrays with the rep, fval, epoch
        and frame number (within teh specific epoch and rep) of each datapoint/frame.
        Saves the data and calculates histogram.
        """
        self.cfg = cfg
        self.fval = fval
        self.frms = frms
        self.reps = reps
        self.epcs = epcs
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        self.u_reps = []
        for e in self.u_epcs:
            self.u_reps.append(np.unique(reps[self.epcs == e]))

        self.nextepoch = self.u_epcs[-1]+1
        # make sure this is not an ignored epoch
        while (self.nextepoch in cfg.ignore_epcs):
            self.nextepoch += 1

        self._make_hist()

    @classmethod
    def fromReadData(cls, cfg: Any, load_fval: bool):
        """
        Factory method to easily load data and make the object
        """
        fval, _, frms, reps, epcs = inout.load_data(cfg, load_fval)
        return cls(cfg, fval, frms, reps, epcs)

    def update(self,
               fval: NDArray[np.float_],
               frms: NDArray[np.int_],
               reps: NDArray[np.int_],
               epcs: NDArray[np.int_]) -> None:
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
            self.u_reps.append(np.unique(reps[self.epcs == e]))

        self.make_hist()

    def _make_hist_no_cfg(self,
                          maxbins: int,
                          data_per_bin: int,
                          minval: float,
                          maxval: float,
                          hard_boundaries: bool = False) -> None:
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
        lowest_val = np.min(self.fval)
        self.binmax = min(largest_val, maxval)
        self.binmin = max(lowest_val,  minval)
        # Calculate binsize
        data_within_bounds = np.sum(
            (self.fval > minval) *
            (self.fval < maxval)
        )
        print(f"{data_within_bounds} data points inside boundaries")
        maxbins_calc = data_within_bounds//data_per_bin
        maxbins = min(maxbins_calc, maxbins)
        if (maxbins < 2):
            raise NotEnoughDataError(
                f"Would be using {maxbins} bins.\n"
                "Either there is not enough data or "
                f"data_per_bin is too large (is {data_per_bin}).\n"
                "Also check that boundaries (minval and maxval) are "
                f"set correctly, now using [{minval}, {maxval}].")

        self.binsize = (self.binmax-self.binmin)/maxbins
        print(f"Making {maxbins} bins within boundaries.")
        # Make bin edges
        if (hard_boundaries):
            self.bin_edges = np.linspace(self.binmin,
                                         self.binmax,
                                         maxbins+1)
        else:
            self.bin_edges = np.arange(lowest_val,
                                       largest_val+self.binsize/2,
                                       self.binsize)
        print("Bin size is {:.5g}.".format(
            self.bin_edges[1]-self.bin_edges[0]
        ))
        # Make the histogram
        self.hist, _ = np.histogram(self.fval, bins=self.bin_edges)
        # bin centers are halfway between edges
        self.bin_centers = (self.bin_edges[:-1]+self.bin_edges[1:])/2
        # Make a mask of which bins have higher edge larger than minval
        # AND lower edge lower than maxval
        self.hist_mask = (self.bin_edges[1:] > minval) * \
            (self.bin_edges[:-1] < maxval)

        self.bins_in_bounds = np.sum(self.hist_mask)

        if (not hard_boundaries):
            print(f"Final histogram has {len(self.hist_mask)} bins. "
                  f"{self.bins_in_bounds} of them have their centre "
                  "within the boundaries.")

    def _make_hist(self) -> None:
        self._make_hist_no_cfg(
            self.cfg.maxbins,
            self.cfg.data_per_bin,
            self.cfg.minval,
            self.cfg.maxval
        )

    def make_choices(self, prechoices: int = 0, plot: bool = True) -> Tuple[
            NDArray[np.float_], NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]:
        """
        Uses the histogram to choose bins and returns the bin indices of the choices.
        prechoices is the number of choices already done.
        """
        # Smooth the distribution
        if (self.bins_in_bounds < self.cfg.smooth_window*2):
            raise NotEnoughDataError(
                f"There are only {self.bins_in_bounds} bins within "
                f"boundaries with a smoothing window of {self.cfg.smooth_window}. "
                "For good results there should be at least twice as many bins.\n"
                "Either make more data or set smooth_window in "
                f"the config to at most {self.bins_in_bounds//2}.")
        smoothed = utils.rolling_mean(self.hist)
        # Get a mask to ingore the nan-values in the smoothed distribution
        nanmask = np.isfinite(smoothed)*self.hist_mask
        # Find maxima and minima
        maxims, max_crit = find_peaks(smoothed[nanmask],
                                      **self.cfg.peak_options)
        minims, min_crit = find_peaks(-smoothed[nanmask],
                                      **self.cfg.peak_options)

        maxh = np.max(self.hist[self.hist_mask])
        minh = np.min(self.hist[self.hist_mask])
        # If not ends have sampled been, make sure min height set to zero is
        if (np.max(self.fval) < self.cfg.maxval or np.min(self.fval) > self.cfg.minval):
            minh = 0

        # The criteria for choices is half way between max and min heights
        crith = (maxh-minh)*self.cfg.choice_crit+minh

        choices = []

        # convert from masked indices to unmasked
        indexes = np.arange(len(self.hist))[nanmask]
        for p in minims:
            if (self.hist[nanmask][p] < crith):
                choices.append(indexes[p])

        zero_mask = self.hist_mask*(self.hist != 0)
        indexes = np.arange(len(self.hist))[zero_mask]

        # Helper booleans
        # is either extrema array empty
        empty_extr = maxims.size == 0 or minims.size == 0
        # Is the  first/last extrema a minima or maxima
        max_first = maxims[0] < minims[0]
        max_last = maxims[-1] > minims[-1]
        # Is the  first/last minima over crith
        first_min_ignored = minims.size > 0 and self.hist[nanmask][minims[0]] >= crith
        last_min_ignored = minims.size > 0 and self.hist[nanmask][minims[-1]] >= crith

        # If the right edge is low enough, check if we add it to choices
        if (self.hist[self.hist_mask][-1] < crith):
            # if   no maxima or minima     or last extr is max        or last min ignored
            if (empty_extr or max_last or last_min_ignored):
                choices.append(indexes[-1])
        # Now same for left edge
        if (self.hist[self.hist_mask][0] < crith):
            # if   no maxima or minima     or first extr is max        or first min ignored
            if (empty_extr or max_first or first_min_ignored):
                choices.append(indexes[0])

        # weights go linearily from 1 at no data 0 at at crith
        weights = [(crith-self.hist[c]) for c in choices]
        weights /= np.sum(weights)
        srt_ind = np.argsort(-weights)
        # sort choices
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

        if (plot):
            self._plot_choices(crith, smoothed, choices,
                               nanmask, maxims, minims)
        return self.choose_frames(choices)

    def choose_frames(self, chosen_bins: List[int]) -> Tuple[
            NDArray[np.float_], NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]:
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

        # lists for value, rpoch, rep and frame
        v, e, r, f = [], [], [], []
        for bi, ci in zip(*np.unique(chosen_bins, return_counts=True)):
            if (bi < len(self.bin_edges)-2):
                vals_in_bin = (
                    self.fval >= self.bin_edges[bi])*(self.fval < self.bin_edges[bi+1])
            else:
                vals_in_bin = (
                    self.fval >= self.bin_edges[bi])*(self.fval <= self.bin_edges[bi+1])

            # How many are in bin
            n_in_bin = np.sum(vals_in_bin)
            # Check if enough are to satisfy minchoice criteria
            if (n_in_bin < self.cfg.minchoice):
                # If not, we choose minchoice closest to bin center
                n_in_bin = self.cfg.minchoice
                cnt = (self.bin_edges[bi]+self.bin_edges[bi+1])/2
                # Sort by distance to cnt and choose the n_in_bin first ones
                vals_in_bin = np.argsort(np.abs(self.fval-cnt))[:n_in_bin]

            # Choose ci from bin
            ndx = self.cfg.rng.choice(n_in_bin, size=ci,
                                      replace=self.cfg.allow_choice_duplicates)

            v.extend(self.fval[vals_in_bin][ndx])
            e.extend(self.epcs[vals_in_bin][ndx])
            r.extend(self.reps[vals_in_bin][ndx])
            f.extend(self.frms[vals_in_bin][ndx])

        return v, e, r, f

    def print_choices(self,
                      val: NDArray[np.float_],
                      epc: NDArray[np.int_],
                      rep: NDArray[np.int_],
                      frm: NDArray[np.int_]) -> None:
        print("Final choices:")
        for v, e, r, f in zip(val, epc, rep, frm):
            print(f"frm {f}, rep {r}, epc {e}, fval={v}")

    def plot_hist(self) -> None:
        """
        Make a histogram of the current total data, as well as maxepochs latest
        OR cumulative histograms after each epoch, with labels in legend only for
        maxepoch latest.
        """
        fig, ax = plt.subplots(1)

        maxepochs = -self.cfg.histogram_max_epochs
        if (maxepochs == 1):
            maxepochs = 0
        # if doing cumulative, also plot the epochs previous to maxepochs, but without labels
        if (self.cfg.cumulative_histogram):
            for i in self.u_epcs[:maxepochs]:
                _ep_hist, _ = np.histogram(
                    self.fval[self.epcs <= i], bins=self.bin_edges)
                ax.plot(self.bin_centers, _ep_hist,
                        "--", color="C0", alpha=0.5)

        # Plot maxepochs latest epochs
        for i in self.u_epcs[maxepochs:]:
            # The mask is <= for cumulative and == for non cumulative
            _mask = self.epcs <= i if (
                self.cfg.cumulative_histogram) else self.epcs == i
            # histogram of masked values
            _ep_hist, _ = np.histogram(self.fval[_mask], bins=self.bin_edges)
            # Plot with label
            ax.plot(self.bin_centers, _ep_hist, "--", color="C%d" %
                    (i+2), label="Epoch %d" % i)

        # Total histogram, with alpha=0.5 to make it lighter
        ax.plot(self.bin_centers, self.hist, color="C0", alpha=0.5)
        # Total histogram only within teh boundaries without alpha
        ax.plot(self.bin_centers[self.hist_mask],
                self.hist[self.hist_mask], color="C0", label="Total")
        # Show starting value and boundaries
        ax.axvline(self.cfg.startval, linestyle="-.",
                   color="C1", label="Start", alpha=0.5)
        ax.axvline(self.cfg.minval, linestyle="-.",
                   color="C2", label="Boundaries")
        ax.axvline(self.cfg.maxval, linestyle="-.", color="C2")
        # legend, on teh right side outside of plot
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        fig.set_size_inches(9, 7)
        fig.tight_layout()
        # Make directory and save fig
        outfile = (
            self.cfg.fig_output_dir /
            ("epoch%02d" % self.u_epcs[-1]) /
            "hist.png"
        )
        outfile.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outfile)
        plt.close(fig)

    def _plot_choices(self, crith, smoothed, choices, nanmask, maxims, minims) -> None:
        fig, ax = plt.subplots(1)

        ax.axhline(crith, linestyle="--", c="r", alpha=0.5,
                   label="Max height for choices")
        ax.plot(self.bin_centers, self.hist, color="C0", alpha=0.5)
        ax.plot(self.bin_centers[self.hist_mask],
                self.hist[self.hist_mask], color="C0", label="data")
        ax.plot(self.bin_centers, smoothed, color="C1", label="Smoothed")
        ax.plot(self.bin_centers[nanmask][minims],
                smoothed[nanmask][minims], "rv", label="Minima")
        ax.plot(self.bin_centers[nanmask][maxims],
                smoothed[nanmask][maxims], "bv", label="Maxima")
        ax.plot(self.bin_centers[choices],
                self.hist[choices], "g^", label="Choices")
        ax.legend()
        fig.set_size_inches(8, 7)
        fig.tight_layout()
        outfile = (
            self.cfg.fig_output_dir /
            ("epoch%02d" % self.u_epcs[-1]) /
            "choices.png"
        )
        outfile.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outfile)
        plt.close(fig)
