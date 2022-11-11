import pathlib
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt

from .. import inout

# Type hints
from typing import Tuple, Dict
import argparse


def leaf_depths(G: nx.DiGraph,
                node: Tuple[int, int] = (0, 0),
                leafs: Dict[Tuple[int, int], int] = {},
                depth=0):
    """
    Recurse on the directed graph G starting from node, returning a dictionary of all the found leafs,
    as well as their depths, with depth being the depth of the current node.
    """
    if (node == (0, 0) and not G.has_node((0, 0))):
        childs = [n for n in G.nodes if n[0] == 0]
        if (len(childs) == 0):
            childs = [n for n in G.nodes if n[0] == 1]
    else:
        childs = list(G.successors(node))

    if (len(childs) == 0):
        leafs[node] = depth
    else:
        for child in childs:
            leaf_depths(G,
                        node=child,
                        leafs=leafs,
                        depth=depth+1)
    return leafs


def make_graph(data: dict, include_start=True):
    """
    Make a directed graph of all the simulation, with edges from the origin to the child simulations
    """
    G = nx.DiGraph()
    for e in data["fval"]:
        for r in data["fval"][e]:
            # Add the node
            G.add_node((e, r))
            # If data was not read properly, we continue without adding edge
            if (not data["origin"][e][r]):
                continue
            # If include_start==False, we don't do anything for the first epoch
            if ((not include_start) and e == 1):
                continue

            # Otherwise we are good to go
            og = data["origin"][e][r]
            parent = (og["epc"], og["rep"])
            # Prevent errors by making sure the parent is in the graph
            # Should only happen with the first epoch, missing origin files
            # or ignored simulations
            if (not G.has_node(parent)):
                G.add_node(parent)
            # Finally add edge
            G.add_edge((parent), (e, r))
    return G


def plot_graph(fout: pathlib.Path, G: nx.DiGraph, data: dict):
    fig, ax = plt.subplots(1)
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")

    max_frm = 0
    for e, r in pos:
        if (e == 0):
            continue
        max_frm = max(data["start_frm"][e][r], max_frm)

    default_size = np.array((50, 12))

    for parent, child in G.edges:
        posA = np.array(pos[parent]) - (0, default_size[1]/2)
        posB = np.array(pos[child]) + (0, default_size[1]/2)
        posA2 = posA-(0, 1.5*default_size[1])
        posB2 = posB+(0, 1.5*default_size[1])
        path = mpl.path.Path([posA, posA2, posB2, posB],
                             [1, 4, 4, 4])

        arr = mpl.patches.PathPatch(path, fc="none")
        ax.add_patch(arr)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        name="MyBlues", colors=[mpl.colormaps["Blues"](0.05), mpl.colormaps["Blues"](0.75)])

    maxs, mins = (-np.inf, -np.inf), (np.inf, np.inf)
    norm = mpl.colors.Normalize(0, max_frm)
    for node in pos:
        size = default_size.copy()
        if (node[0] == 0):
            size[0] *= 2
        pos[node]-size/2
        maxs = np.maximum(maxs, pos[node]+size)
        mins = np.minimum(mins, pos[node]-size)
        if (node[0] > 0):
            color = cmap(norm(data["start_frm"][node[0]][node[1]]))
        else:
            color = cmap(0)
        rec = mpl.patches.Rectangle(pos[node]-size/2,
                                    *(size),
                                    facecolor=color,
                                    edgecolor="k")
        ax.add_patch(rec)
        ax.text(*pos[node], f"{node[0] if node[0]!=0 else 'start'}",
                ha="center", va="center")

    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    fig.set_size_inches(26, 10)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(
        norm, cmap), ax=ax, location="right", fraction=0.02, aspect=30, pad=0.005)
    cbar.set_label("Number of frames at start of simulation")

    ax.axis("off")

    fig.tight_layout()
    fout.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fout)
    plt.clf()


def ancestry(args: argparse.Namespace) -> None:
    data = inout.load_extract_data(args.cfg, args.doignore)
    G = make_graph(data, include_start=args.include_start)
    plot_graph(args.output, G, data)
