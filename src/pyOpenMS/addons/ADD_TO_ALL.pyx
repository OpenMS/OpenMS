
from libc.stddef cimport ptrdiff_t
from libc.stdint cimport *
from libcpp.map cimport map as libcpp_map
cimport numpy as np
import numpy as np
ctypedef libcpp_vector[ double ] _DoubleList
ctypedef libcpp_vector[ int ] _IntList


# ArrayWrapper classes (ArrayWrapperFloat, ArrayWrapperDouble, etc.) are now
# provided by autowrap. They support the Python buffer protocol for efficient
# numpy array conversion. See autowrap/data_files/autowrap/ArrayWrappers.pyx