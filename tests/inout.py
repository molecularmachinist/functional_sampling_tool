import unittest
import io
import time
from contextlib import redirect_stdout

from functional_sampling_tool import inout
from functional_sampling_tool import exceptions

import numpy as np

from data import epoch_example_dir, set_xtc_mtimes, load_cfg, temp_example_dir, run_in_dir, silence_function


def capture_output_wrapper(func, *args, **kwargs):
    f = io.StringIO()
    with redirect_stdout(f):
        func(*args, **kwargs)
    return f.getvalue()


class CheckNum(unittest.TestCase):
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


class GetDataFromArchive(unittest.TestCase):
    """
    TestCase for the get_data_from_archive function
    """

    def setUp(self):
        set_xtc_mtimes()
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


class CleanLatestEpoch(unittest.TestCase):
    """
    Test case for testing clean_latest_epoch
    """

    def setUp(self):
        self.tmp_dir = temp_example_dir()

    def test_0_success(self):
        """
        Test that the cleaning succeeds when no xtc are present
        """
        self.assertTrue((self.tmp_dir / "none_present1" / "epoch01").is_dir(),
                        msg="If this fails, the test setup failed, not the test")
        wrapper = run_in_dir(self.tmp_dir / "none_present1")
        clean_latest_epoch = wrapper(inout.clean_latest_epoch)
        clean_latest_epoch()
        self.assertFalse((self.tmp_dir / "none_present1" / "epoch01").is_dir())

    def test_1_no_epochs(self):
        """
        Test that the cleaning fails when not in correct folder
        """
        files_before = list((self.tmp_dir / "none_present2").rglob("*"))
        files_before.sort()

        wrapper = run_in_dir(self.tmp_dir / "none_present2" / "epoch01")
        clean_latest_epoch = wrapper(inout.clean_latest_epoch)
        self.assertRaises(exceptions.NoEpochsFoundError,
                          clean_latest_epoch)
        wrapper = run_in_dir(
            self.tmp_dir / "none_present2" / "epoch01" / "rep01"
        )
        clean_latest_epoch = wrapper(inout.clean_latest_epoch)
        self.assertRaises(exceptions.NoEpochsFoundError,
                          clean_latest_epoch)

        files_after = list((self.tmp_dir / "none_present2").rglob("*"))
        files_after.sort()

        self.assertListEqual(files_before, files_after)

    def _check_no_complete(self, tmp_dir):
        files_before = list((tmp_dir).rglob("*"))
        files_before.sort()

        wrapper = run_in_dir(tmp_dir)
        clean_latest_epoch = wrapper(silence_function(
            inout.clean_latest_epoch
        ))
        clean_latest_epoch()

        files_after = list((tmp_dir).rglob("*"))
        files_after.sort()

        self.assertListEqual(files_before, files_after)

    def test_2_has_data(self):
        """
        Test that the cleaning does not complete when xtc files are present
        """
        self._check_no_complete(self.tmp_dir / "all_present")

    def test_3_single_xtc(self):
        """
        Test that the cleaning does not complete when a single xtc file is present
        """
        self._check_no_complete(self.tmp_dir / "one_present1")

    def test_4_force_xtc(self):
        """
        Test that the cleaning does complete when a single xtc file is present, but force is True
        """
        expected_output = "Found .*, but will still force clean\n" \
            "Data in epoch[0-9][0-9]+ will be lost forever\n" \
            "You have 5 seconds to cancel with ctrl-C"

        tmp_dir = self.tmp_dir / "one_present2"
        self.assertTrue((tmp_dir / "epoch01").is_dir(),
                        msg="If this fails, the test setup failed, not the test")
        wrapper = run_in_dir(tmp_dir)
        clean_latest_epoch = wrapper(inout.clean_latest_epoch)
        t1 = time.monotonic()
        output = capture_output_wrapper(clean_latest_epoch, force=True)
        elapsed = time.monotonic()-t1
        self.assertGreaterEqual(elapsed, 5.0,
                                msg="The function should give 5 seconds to cancel")
        self.assertRegex(output, expected_output)
        self.assertFalse((tmp_dir / "epoch01").is_dir())


if __name__ == "__main__":
    unittest.main(
        argv=["inout_tests", "-v"]
    )
