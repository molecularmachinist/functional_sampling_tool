import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pathlib
import math
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from . import choosing
from . import inout

# Type hints
from typing import Any, Tuple, List, Dict
from numpy.typing import NDArray

colors = ['r', 'b', 'g', 'orange', 'cyan', 'indigo',
          'purple', 'magenta', 'brown', 'k', 'y']
colors += ["sandybrown", "deeppink", "darkslategrey",
           "indianred", "lawngreen", "lightseagreen"]
colors = np.array(colors)


class ClusterChooser(choosing.FrameChooser):

    def __init__(self,
                 cfg: Any,
                 fval: NDArray[np.float_],
                 coords: NDArray[np.float_],
                 frms: NDArray[np.int_],
                 reps: NDArray[np.int_],
                 epcs: NDArray[np.int_]):
        """
        Takes the config and data as 4 N-length arrays with the rep, fval, epoch
        and frame number (within teh specific epoch and rep) of each datapoint/frame.
        Saves the data and calculates histogram.
        """
        self.cfg = cfg
        self.coords = coords
        self.reps = reps
        self.fval = fval
        self.epcs = epcs
        self.frms = frms
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        self.u_reps = []
        for e in self.u_epcs:
            self.u_reps.append(np.unique(reps[self.epcs == e]))

        self.nextepoch = self.u_epcs[-1]+1

        if (len(self.u_epcs) >= self.cfg.epochs_pre_clust):
            self._make_hist()

        self.plain_chooser = choosing.FrameChooser(cfg, fval, frms, reps, epcs)

    @classmethod
    def fromReadData(cls, cfg: Any, load_fval: bool):
        """
        Factory method to easily load data and make the object
        """
        fval, crds, frms, reps, epcs = inout.load_data(cfg, load_fval)
        return cls(cfg, fval, crds, frms, reps, epcs)

    def _make_hist(self) -> None:
        self._make_hist_no_cfg(
            self.cfg.clust_maxbins,
            self.cfg.clust_data_per_bin,
            self.cfg.minval,
            self.cfg.maxval,
            hard_boundaries=True
        )
        self.clust_hist_indexes = np.digitize(self.fval, bins=self.bin_edges)
        self.clust_unique_bins, self.clust_unique_bin_counts = np.unique(
            self.clust_hist_indexes, return_counts=True)

    def make_choices(self, prechoices: int = 0, plot: bool = True) -> Tuple[
            NDArray[np.float_], NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]:
        if (len(self.u_epcs) < self.cfg.epochs_pre_clust):
            return self.plain_chooser.make_choices(prechoices, plot)

        # Allocate array for labels, default value -1 for non labeled frames.
        clusters = np.full(self.fval.shape, -1, dtype=int)

        natoms = self.coords.shape[1]
        dims = self.coords.shape[2]
        clust_results = {}
        # Iterate over bins
        for ci, cc in zip(self.clust_unique_bins, self.clust_unique_bin_counts):
            # check the bin is not "over the edges" and there is enough data to go on
            if (ci != 0 and ci != len(self.bin_edges) and cc >= self.cfg.clust_data_per_bin):
                # A mask of which values are in this bin
                indxs = self.clust_hist_indexes == ci
                nframes = np.sum(indxs, dtype=int)
                # Get the coordinates in this bin
                crds = self.coords[indxs].reshape((nframes, natoms*dims))
                # Do the actual clustering
                clust_results[ci] = make_clusters(
                    coords=crds,
                    maxclust=self.cfg.maxclust,
                    tol=self.cfg.clust_tol,
                    rng=self.cfg.rng
                )
                # Make cluster labels unique over different bins by adding 100 times the bin index to the labels,
                # unless the maxclust is larger than that. In that case we add maxclust times the bin index.
                clust_results[ci]["labels"] += ci*max(100, self.cfg.maxclust)
                clusters[indxs] = clust_results[ci]["labels"]

        # count data in each cluster
        clusts, counts = np.unique(clusters[clusters >= 0], return_counts=True)
        srt_clusts = clusts[np.argsort(counts)]

        nchoices = int(self.cfg.clust_choice_frac*self.cfg.N)

        # Choose either the number appointed by clust_choice_frac, or at maximum half the total clusters
        maxchoice = min(int(self.cfg.clust_choice_frac * self.cfg.N),
                        len(srt_clusts)//2)
        choices = list(srt_clusts[:maxchoice])

        # If we can choose still more, do so by duplicate choosing
        if (len(choices) < nchoices):
            weights = [(np.max(counts)-counts[c == clusts][0])
                       for c in choices]
            weights /= np.sum(weights)
            # Each point has been added once
            len_choice = len(choices)
            # Add more depending on weight. Multiplication is floored, so between 0 and len_choice-1 too few are added
            for i in range(len_choice):
                for j in range(math.floor(weights[i]*(nchoices-len_choice))):
                    choices.append(choices[i])

            for i in range(nchoices-len(choices)):
                choices.append(choices[i])

        print(f"Chose {len(choices)} clusters out of "
              f"{len(clusts)}, with {nchoices} max allowed")
        print(f"{len(np.unique(choices))} unique clusters")

        if (plot):
            # Make sure the plot dir exists
            (self.cfg.fig_output_dir /
             ("epoch%02d" % self.u_epcs[-1]) /
             "clusters"
             ).mkdir(parents=True, exist_ok=True)
            # Make a dictionary of choices and their counts (multiplicities)
            choice_counts = {ci: cc for ci, cc in
                             zip(*np.unique(choices, return_counts=True))}
            # Plot all bins
            for ci in clust_results:
                plot_clust(clust_results[ci],
                           choice_counts,
                           (self.bin_edges[ci-1], self.bin_edges[ci]),
                           plotname=self.cfg.fig_output_dir /
                           ("epoch%02d" % (self.u_epcs[-1])) /
                           "clusters" / ("choices_%d.png" % ci))

        # Choose the actual frames
        fval, epcs, reps, frms = self.choose_frames(choices, clusters)

        v, e, r, f = self.plain_chooser.make_choices(
            prechoices=prechoices+len(choices),
            plot=plot
        )
        fval += v
        epcs += e
        reps += r
        frms += f

        return fval, epcs, reps, frms

    def choose_frames(self, chosen_clusts: List[int], clusters: NDArray[np.int_]) -> Tuple[
            NDArray[np.float_], NDArray[np.int_], NDArray[np.int_], NDArray[np.int_]]:
        """ Input parameters:
                - chosen_clusts : Length N list of the labels of the clusters that have been chosen.
                                  Duplicates (starting from the same bin) are simply many times in the list.
                - clusters      : Array of the cluster labels for all simulated frames
            Output:
                - v : Length N list of the fvals for each new rep
                - e : Length N list of the epochs each new rep comes from
                - r : Length N list of the rep within the epoch each new rep comes from
                - f : Length N list of the frame within the rep each new rep comes from

            Other "inputs":
                Saved already in the object:
                - fval      : An array containing all the fvals from all the simulations so far.
                - epcs      : An array of same shape as fval, with each element being the corresponding
                              epoch each value comes from.
                - reps      : Same as epcs, but with the corresponding rep.
                - frms      : Same as epcs and reps, but with the information of the frame.
        """
        # lists for value, rpoch, rep and frame
        v, e, r, f = [], [], [], []
        for ci in chosen_clusts:
            vals_in_clust = ci == clusters
            ndx = self.cfg.rng.choice(np.sum(vals_in_clust))

            v.append(self.fval[vals_in_clust][ndx])
            e.append(self.epcs[vals_in_clust][ndx])
            r.append(self.reps[vals_in_clust][ndx])
            f.append(self.frms[vals_in_clust][ndx])

        return v, e, r, f

    def plot_hist(self) -> None:
        self.plain_chooser.plot_hist()
        if (len(self.u_epcs) < self.cfg.epochs_pre_clust):
            return
        clust_hist_indexes = np.digitize(self.fval, bins=self.bin_edges)
        clust_unique_bins, clust_unique_bin_counts = np.unique(
            clust_hist_indexes, return_counts=True)

        clust_histx, clust_histy = [], []
        for ci, cc in zip(clust_unique_bins, clust_unique_bin_counts):
            clust_histx.append(np.mean(self.fval[ci == clust_hist_indexes]))
            clust_histy.append(cc)

        fig, ax = plt.subplots(1)

        ax.plot(clust_histx, clust_histy, color="C0")
        ax.axvline(self.cfg.startval, linestyle="-.",
                   color="C1", label="Start", alpha=0.5)
        ax.axvline(self.cfg.minval, linestyle="-.",
                   color="C2", label="Boundaries")
        ax.axvline(self.cfg.maxval, linestyle="-.", color="C2")

        for c in self.bin_edges:
            ax.axvline(c, linestyle="-.", alpha=.1, color="k")

        ax.legend()
        fig.set_size_inches(9, 7)
        fig.tight_layout()
        outfile = (self.cfg.fig_output_dir /
                   ("epoch%02d" % self.u_epcs[-1]) /
                   "hist_clust.png")
        outfile.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outfile)
        plt.close(fig)


def make_clusters(coords: NDArray[np.float_],
                  maxclust: int,
                  tol: float,
                  rng: np.random.Generator) -> Dict[str, NDArray]:
    """
    Takes a shape(n,m) array of coordinates and return shape(n) labels of clusters.
    The labels should be integers between 0 and maxclust-1.

    Input parameters:
        - coords   : shape(n,m) array of the coordinates to cluster.
        - maxclust : maximum allowed number of clusters
        - tol      : the tolerance for the choosing criteria. A float between 0 and 1.
    Output:
        - results : A dictionary with the labels under the \"labels\" key and all data needed for plotting
                    included under different keys.
    """
    shuffled_data = rng.permutation(coords)

    # call pca
    ncomponents = 2
    pca = PCA(n_components=ncomponents, whiten=True)
    pca.fit(shuffled_data)

    print(f"First {ncomponents} PCA components explain {np.sum(pca.explained_variance_ratio_)*100} % of variance")

    # Get the two first PCs for each datapoint
    xpca = pca.transform(coords)

    # how to choose clusters?
    bic = []
    aic = []

    for n_cluster in range(1, maxclust+1):
        gm = GaussianMixture(n_components=n_cluster,
                             random_state=0).fit(xpca)
        bic.append(gm.bic(xpca))
        aic.append(gm.aic(xpca))

    bic = np.array(bic)
    aic = np.array(aic)

    rbic = bic-bic.min()
    rbic /= rbic[0]
    raic = aic-aic.min()
    raic /= raic[0]

    nclust = -1
    # find the first occurence, where raic and rbic are both below tol
    for i, (ra, rb) in enumerate(zip(raic, rbic)):
        if (ra < tol and rb < tol):
            nclust = i+1
            break
    # If nclust is not yet found, take the number with minimum product
    if (nclust < 0):
        nclust = np.argmin(raic*rbic)+1

    print(f"making {nclust} clusters")

    # now train the final GMM
    gm = GaussianMixture(n_components=nclust, random_state=0).fit(xpca)

    results = {}

    results["maxclust"] = maxclust
    results["nclust"] = nclust

    results["bic"] = bic
    results["aic"] = aic
    results["rbic"] = rbic
    results["raic"] = raic
    results["labels"] = gm.predict(xpca)
    results["reduced_dims"] = xpca

    return results


def plot_clust(results: Dict[str, NDArray], choices: Dict[int, int], fval_range: Tuple[float, float], plotname: pathlib.Path):
    """
    Takes a shape(n,m) array of coordinates and return shape(n) labels of clusters.
    The labels should be integers between 0 and maxclust.

    Input parameters:
        - results  : A dictionary with the results from make_clusters.
        - choices  : A dictionary with the chosen labels as keys and multiplicity of choices as values.
        - choices  : A tuple with the boundaries of the bin in function values.
        - plotname : A pathlib.Path where the figure will be saved.
    """
    global colors
    nclust = results["nclust"]
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    axes[0, 0].plot(range(1, results["maxclust"]+1), results["bic"])
    axes[0, 0].plot(range(1, results["maxclust"]+1), results["aic"])
    axes[0, 0].plot(nclust, results["bic"][nclust-1], "rx")
    axes[0, 0].plot(nclust, results["aic"][nclust-1], "rx")

    axes[1, 0].plot(range(1, results["maxclust"]+1), results["rbic"])
    axes[1, 0].plot(range(1, results["maxclust"]+1), results["raic"])
    axes[1, 0].plot(nclust, results["rbic"][nclust-1], "rx")
    axes[1, 0].plot(nclust, results["raic"][nclust-1], "rx")
    fig.tight_layout()

    unique_lbls, inv_lbl_ndxes = np.unique(results["labels"],
                                           return_inverse=True)
    cols = colors[inv_lbl_ndxes % len(colors)]

    axes[0, 1].scatter(*results["reduced_dims"].T,
                       s=1, c=cols)
    handles = []
    for i, lbl in enumerate(unique_lbls):
        label = str(lbl)
        if lbl in choices:
            if (choices[lbl] > 1):
                label += f" ({choices[lbl]})"
            else:
                label += " (*)"
        handles.append(
            mpatches.Patch(color=colors[i % len(colors)], label=label))

    axes[1, 1].legend(handles=handles, loc="upper left")
    axes[1, 1].axis("off")

    fig.suptitle("fval range ({:.5g} to {:.5g})".format(*fval_range))
    fig.tight_layout()

    fig.savefig(plotname)

    plt.clf()
