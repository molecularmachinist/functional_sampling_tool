import unittest

import MDAnalysis as mda
import numpy as np

from functional_sampling_tool import _ctransformations
from functional_sampling_tool.transformations import make_whole


from data import make_whole_data


def _wrap_in_subtests(func):
    def subtest_wrapper(self):
        for name in self.systems:
            with self.subTest(name):
                func(self, self.systems[name])
    return subtest_wrapper


class MakeWholeData:
    def __init__(self, name, dat) -> None:
        self.crds = dat["coordinates"]
        self.bonds = dat["bonds"]
        self.dims = dat["dimensions"]
        self.crds_whole = dat["unwrapped_coordinates"]
        self.starters = dat["starters"]
        self.eps = dat["eps"]
        self.name = name
        u = mda.Universe.empty(self.crds.shape[0], trajectory=True)
        u.atoms.positions = self.crds
        u.dimensions = self.dims
        u.add_TopologyAttr('bonds', self.bonds)
        bonds = _ctransformations.traverse_mol(u.atoms.indices,
                                               u.bonds.to_indices(),
                                               self.starters)

        self.made_whole = make_whole(u.trajectory[0],
                                     bonds,
                                     u.atoms.indices).positions


class MakeWhole(unittest.TestCase):
    """
    TestCase for the make_whole function
    """

    def setUp(self):
        systems = {}
        data = make_whole_data()
        for sys in data:
            dat = data[sys]
            systems[sys] = MakeWholeData(sys, dat)

        self.systems = systems

    @_wrap_in_subtests
    def test_0_shape(self, dat: MakeWholeData):
        """
        Test that the shape is unchanged 
        """
        self.assertTupleEqual(
            dat.crds.shape,
            dat.made_whole.shape
        )

    @_wrap_in_subtests
    def test_1_values(self, dat: MakeWholeData):
        """
        Test that the output matches
        """
        diffs = (np.abs(dat.made_whole-dat.crds_whole))
        max_diff = np.max(diffs)
        min_diff = np.min(diffs)
        self.assertLess(
            max_diff, dat.eps,
            f"{max_diff=} ({min_diff=})"
        )


if __name__ == "__main__":
    unittest.main(
        argv=["transformation_tests", "-v"]
    )
