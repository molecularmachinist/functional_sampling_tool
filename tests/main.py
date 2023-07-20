import unittest
from transformations import MakeWhole
from inout import CheckNum, GetDataFromArchive, CleanLatestEpoch


if __name__ == "__main__":
    unittest.main(
        argv=["fst_tests", "-v"]
    )
