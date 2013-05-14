import nose
import sys
import os

import env

#os.environ["LD_LIBRARY_PATH"] = os.path.join(env.OPEN_MS_BUILD_DIR, "lib")
#sys.path.insert(0, ".")

import nose
import sys
import os

import env

#os.environ["LD_LIBRARY_PATH"] = os.path.join(env.OPEN_MS_BUILD_DIR, "lib")
#sys.path.insert(0, ".")

nose.run_exit(argv=["-v", "-w", "tests/memoryleaktests", "-s"])
