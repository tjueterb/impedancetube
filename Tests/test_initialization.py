import numpy as np
import sys

try:
    # if directory is the root directory:
    # adding '.' to the path seems to be necessary for debugging this file in VS Code
    sys.path.append('.')
    import impedancetube as imp

except:
    # if directory is test directory:
    sys.path.append('..')
    import impedancetube as imp


def test_two_load_method():
    msm = imp.Measurement()

    assert(np.allclose(msm.rho, 1.2019904485758144))
    assert(np.allclose(msm.c, 343.23719141193016))


if __name__ == "__main__":
    test_two_load_method()
    print("\nEverything passed!\n")
