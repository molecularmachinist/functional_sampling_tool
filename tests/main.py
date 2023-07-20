import unittest
import transformations
import inout


if __name__ == "__main__":
    unittest.main(
        defaultTest=["transformations", "inout"],
        argv=["fst_tests", "-v"]
    )
