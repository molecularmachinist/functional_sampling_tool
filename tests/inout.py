import copy
import unittest

from functional_sampling_tool import inout
from functional_sampling_tool import exceptions

import numpy as np

from data import epoch_example_dir, load_cfg


class TestCase_check_num(unittest.TestCase):
    """
    TestCase for the inout.check_num function
    """

    def test_epoch(self):
        """
        Test that the single epoch is are correctly found
        """
        self.assertListEqual(
            inout.check_num(epoch_example_dir / "epoch"),
            [1]
        )

    def test_rep(self):
        """
        Test that numbers of reps are correctly found when one is missing from between
        """
        self.assertListEqual(
            inout.check_num(epoch_example_dir / "epoch01" / "rep"),
            [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16]
        )

    def test_no_rep(self):
        """
        Test that empty list is returned when nothing is found
        """
        self.assertListEqual(
            inout.check_num(epoch_example_dir / "rep"),
            []
        )


class TestCase_get_data_from_archive(unittest.TestCase):
    """
    TestCase for the get_data_from_archive function
    """

    def setUp(self):
        self.cfg = load_cfg()

    def test_all_good(self):
        """
        Test that data is loaded correctly when nothing is changed
        """
        load_dir = epoch_example_dir / "epoch01" / "rep10"
        with np.load(load_dir / "fval_data.npz") as npz:
            real_fval = npz["fval"]
            real_crds = npz["crd"]
        fval, crds = inout.get_data_from_archive(load_dir, self.cfg)
        self.assertTrue(
            np.array_equal(
                fval,
                real_fval
            )
        )
        self.assertTrue(
            np.array_equal(
                crds,
                real_crds
            )
        )

    def test_changed_modtime(self):
        """
        Test that LoadError is raised when xtc modtime does not match
        """
        self.assertRaisesRegex(inout.LoadError,
                               "Modification time of .* does not match",
                               inout.get_data_from_archive,
                               epoch_example_dir / "epoch01" / "rep02", self.cfg)

    def test_missing_archive(self):
        """
        Test that FileNotFoundError is raised when archive is missing
        """
        self.assertRaisesRegex(FileNotFoundError,
                               "No such file or directory: '.*fval_data.npz'",
                               inout.get_data_from_archive,
                               epoch_example_dir / "epoch01" / "rep01", self.cfg)

    def test_changed_sel(self):
        """
        Test that LoadError is raised when selections change
        """
        sel = self.cfg.sel
        self.cfg.sel = self.cfg.sel.universe.select_atoms("resid 3")
        self.assertRaisesRegex(inout.LoadError,
                               "Selections in .* do not match",
                               inout.get_data_from_archive,
                               epoch_example_dir / "epoch01" / "rep10", self.cfg)
        self.cfg.sel = sel

        sel = self.cfg.sel_clust
        self.cfg.sel_clust = self.cfg.sel.universe.select_atoms("resid 3")
        self.assertRaisesRegex(inout.LoadError,
                               "Selections in .* do not match",
                               inout.get_data_from_archive,
                               epoch_example_dir / "epoch01" / "rep10", self.cfg)
        self.cfg.sel_clust = sel

    def test_changed_transforms(self):
        """
        Test that LoadError is raised if transformations are turned off
        """
        self.cfg.unwrap_mols = False
        self.assertRaisesRegex(inout.LoadError,
                               "Trajectory transformations changed",
                               inout.get_data_from_archive,
                               epoch_example_dir / "epoch01" / "rep10", self.cfg)
        self.cfg.unwrap_mols = True

        self.cfg.mols_in_box = False
        self.assertRaisesRegex(inout.LoadError,
                               "Trajectory transformations changed",
                               inout.get_data_from_archive,
                               epoch_example_dir / "epoch01" / "rep10", self.cfg)
        self.cfg.mols_in_box = True


if __name__ == "__main__":
    unittest.main()
