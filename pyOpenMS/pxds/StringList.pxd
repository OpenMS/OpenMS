from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
# see OpenMS/DATASTRUCTURES/ListUtils.h
ctypedef libcpp_vector[ String ] StringList

