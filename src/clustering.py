import pathlib
import numpy as np
import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import os
import math
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from . import choosing
from . import inout

# Type hints
from typing import Any, Tuple, Union, List
from numpy.typing import NDArray

colors  = ['r','b','g','orange','cyan','indigo','purple','magenta','brown','k','y']
colors += ["sandybrown", "deeppink","darkslategrey","indianred","lawngreen","lightseagreen"]

rng = np.random.default_rng()


class ClusterChooser(choosing.FrameChooser):


    def __init__(self, cfg: Any, fval: NDArray[np.float_], coords: NDArray[np.float_], frms: NDArray[np.int_], reps: NDArray[np.int_], epcs: NDArray[np.int_]):
        """
        Takes the config and data as 4 N-length arrays with the rep, fval, epoch
        and frame number (within teh specific epoch and rep) of each datapoint/frame.
        Saves the data and calculates histogram.
        """
        self.cfg  = cfg
        self.coords = coords
        self.reps = reps
        self.fval = fval
        self.epcs = epcs
        self.frms = frms
        # Also get the unique epochs and reps per epoch
        self.u_epcs = np.unique(epcs)
        self.u_reps = []
        for e in self.u_epcs:
            self.u_reps.append(np.unique(reps[self.epcs==e]))

        self.nextepoch = self.u_epcs[-1]+1

        if(len(self.u_epcs)>self.cfg.epochs_pre_clust):
            self._make_hist()

        self.plain_chooser = choosing.FrameChooser(cfg, fval, frms, reps, epcs)

    @classmethod
    def fromReadData(cls, cfg: Any, load_fval: bool):
        """
        Factory method to easily load data and make the object
        """
        fval, crds, frms, reps, epcs = inout.load_data(cfg,load_fval)
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
        self.clust_unique_bins,self.clust_unique_bin_counts = np.unique(self.clust_hist_indexes,return_counts=True)


    def make_choices(self, prechoices: int = 0, plot: bool = True) -> Tuple[NDArray[np.float_],NDArray[np.int_],NDArray[np.int_],NDArray[np.int_]]:
        if(len(self.u_epcs)<=self.cfg.epochs_pre_clust):
            return self.plain_chooser.make_choices(prechoices,plot)
        if(plot):
            os.makedirs("figs/epoch%02d"%(self.u_epcs[-1]),exist_ok=True)
            os.makedirs("figs/epoch%02d/clusters"%(self.u_epcs[-1]),exist_ok=True)
        clusters = np.full(self.fval.shape, -1,dtype=int)

        natoms = self.coords.shape[1]
        dims = self.coords.shape[2]
        for ci,cc in zip(self.clust_unique_bins,self.clust_unique_bin_counts):
            if(ci!=0 and ci!=len(self.bin_edges) and cc>=self.cfg.clust_data_per_bin):
                indxs = self.clust_hist_indexes==ci
                nframes = np.sum(indxs, dtype=int)
                crds = self.coords[indxs].reshape((nframes,natoms*dims))
                clusters[indxs] = make_clusters(
                    coords=crds,
                    plot=plot,
                    plotname="figs/epoch%02d/clusters/choices_%d.png"%(self.u_epcs[-1],ci),
                    maxclust=self.cfg.maxclust,
                    tol=self.cfg.clust_tol
                    ) + ci*(self.cfg.maxclust)*2

        clusts,counts = np.unique(clusters[clusters>=0], return_counts=True)
        srt_clusts = clusts[np.argsort(counts)]

        nchoices=int(self.cfg.clust_choice_frac*self.cfg.N)
        maxchoice = min(int(self.cfg.clust_choice_frac*self.cfg.N),len(srt_clusts)//2)
        choices = list(srt_clusts[:maxchoice])
        if(len(choices)<nchoices):
            weights = [(np.max(counts)-counts[c==clusts][0]) for c in choices]
            weights /= np.sum(weights)
            # Each point has been added once
            len_choice = len(choices)
            # Add more depending on weight. Multiplication is floored, so between 0 and len_choice-1 too few are added
            for i in range(len_choice):
                for j in range(math.floor(weights[i]*(nchoices-len_choice))):
                    choices.append(choices[i])

            for i in range(nchoices-len(choices)):
                choices.append(choices[i])

        print(f"Chose {len(choices)} clusters out of {len(clusts)}, with {nchoices} max allowed")
        print(f"{len(np.unique(choices))} unique clusters")
        fval,epcs,reps,frms = self.choose_frames(choices,clusters)

        v,e,r,f = self.plain_chooser.make_choices(prechoices+len(choices),plot)
        fval += v
        epcs += e
        reps += r
        frms += f

        return fval, epcs, reps, frms

    def choose_frames(self, chosen_clusts: List[int], clusters: NDArray[np.int_]) -> Tuple[NDArray[np.float_],NDArray[np.int_],NDArray[np.int_],NDArray[np.int_]]:
        """ Input parameters:
                - chosen_clusts : Length N list of the indices of the clusters that have been chosen.
                                  Duplicates (starting from the same bin) are simply many times in the list.
                - clusters      : Array of the cluster indices for all simulated frames
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


        global rng
        # lists for value, rpoch, rep and frame
        v,e,r,f=[],[],[],[]
        for ci in chosen_clusts:
            vals_in_clust = ci==clusters
            ndx = rng.choice(np.sum(vals_in_clust))

            v.append(self.fval[vals_in_clust][ndx])
            e.append(self.epcs[vals_in_clust][ndx])
            r.append(self.reps[vals_in_clust][ndx])
            f.append(self.frms[vals_in_clust][ndx])

        return v,e,r,f


    def plot_hist(self) -> None:
        self.plain_chooser.plot_hist()
        if(len(self.u_epcs)<=self.cfg.epochs_pre_clust):
            return
        clust_hist_indexes = np.digitize(self.fval, bins=self.bin_edges)
        clust_unique_bins,clust_unique_bin_counts = np.unique(clust_hist_indexes, return_counts=True)

        clust_histx, clust_histy = [],[]
        for ci,cc in zip(clust_unique_bins,clust_unique_bin_counts):
            clust_histx.append(np.mean(self.fval[ci==clust_hist_indexes]))
            clust_histy.append(cc)

        plt.plot(clust_histx, clust_histy, color="C0")
        plt.axvline(self.cfg.startval, linestyle="-.", color="C1", label="Start", alpha=0.5)
        plt.axvline(self.cfg.minval, linestyle="-.", color="C2", label="Boundaries")
        plt.axvline(self.cfg.maxval, linestyle="-.", color="C2")

        for c in self.bin_edges:
            plt.axvline(c, linestyle="-.", alpha=.1, color="k")

        plt.legend()
        os.makedirs("%s/epoch%02d"%(self.cfg.fig_output_dir, self.u_epcs[-1]),exist_ok=True)
        plt.savefig("%s/epoch%02d/hist_clust.png"%(self.cfg.fig_output_dir, self.u_epcs[-1]))
        plt.clf()



def make_clusters(coords: NDArray[np.float_], plot: bool = False, plotname: Union[str,pathlib.Path] = "plot.png", maxclust: int = 15, tol: float = 0.1) ->  NDArray[np.int_]:
    global colors
    """
    Takes a shape(n,m) array of coordinates and return shape(n) labels of clusters.
    """
    global rng

    shuffled_data = rng.permutation(coords)

    # call pca
    ncomponents = 2
    pca = PCA(n_components=ncomponents,whiten=True)
    pca.fit(shuffled_data)

    # plot the explained variance

    print(f"First {ncomponents} PCA components explain {np.sum(pca.explained_variance_ratio_)*100} % of variance")

    # 2D projections on the eigenspace_orig
    xpca = pca.transform(coords)

    first_2 = np.vstack((xpca[:,0],xpca[:,1])).T

    # how to choose clusters?
    bic = []
    aic = []

    for n_cluster in range(1,maxclust+1):
        gm = GaussianMixture(n_components=n_cluster, random_state=0).fit(first_2)
        bic.append(gm.bic(first_2))
        aic.append(gm.aic(first_2))

    bic =np.array(bic)
    aic =np.array(aic)

    rbic  = bic-bic.min()
    rbic /= rbic[0]
    raic  = aic-aic.min()
    raic /= raic[0]

    nclust = -1
    for i,(ra,rb) in enumerate(zip(raic, rbic)):
        if(ra<tol and rb<tol):
            nclust=i+1
            break

    if(nclust<0):
        nclust = np.argmin(raic*rbic)+1


    print(f"making {nclust} clusters")

    # now try with GMM
    gm = GaussianMixture(n_components=nclust, random_state=0).fit(first_2)

    if(plot):
        fig, axes = plt.subplots(2,2,figsize = (8,7))
        axes[0,0].plot(range(1,maxclust+1),bic)
        axes[0,0].plot(range(1,maxclust+1),aic)
        axes[0,0].plot(nclust,bic[nclust-1],"rx")
        axes[0,0].plot(nclust,aic[nclust-1],"rx")


        axes[0,1].plot(range(2,maxclust+1),bic[1:]-bic[:-1])
        axes[0,1].plot(range(2,maxclust+1),aic[1:]-aic[:-1])

        axes[1,0].plot(range(1,maxclust+1),rbic)
        axes[1,0].plot(range(1,maxclust+1),raic)
        axes[1,0].plot(nclust,rbic[nclust-1],"rx")
        axes[1,0].plot(nclust,raic[nclust-1],"rx")
        fig.tight_layout()

        if(nclust<len(colors)):
            cmap = mpl_colors.ListedColormap(colors[:nclust],N=nclust)
        else:
            cmap = mpl_colors.ListedColormap(colors,N=nclust)

        labels = gm.predict(first_2)
        axes[1,1].scatter(first_2[:,0],first_2[:,1],s=1,c=labels, cmap=cmap)

        fig.savefig(plotname)

    plt.clf()
    return labels
