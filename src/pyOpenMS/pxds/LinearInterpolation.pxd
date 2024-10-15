from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *

cdef extern from "<OpenMS/ML/INTERPOLATION/LinearInterpolation.h>" namespace "OpenMS::Math":
    
    cdef cppclass LinearInterpolation[KeyType,ValueType]:
        # wrap-doc:
            #  Provides access to linearly interpolated values (and
            #  derivatives) from discrete data points.  Values beyond the given range
            #  of data points are implicitly taken as zero.
            #  
            #  The input is just a vector of values ("Data").  These are interpreted
            #  as the y-coordinates at the x-coordinate positions 0,...,data_.size-1.
            #  
            #  The interpolated data can also be scaled and shifted in
            #  the x-dimension by an affine mapping.  That is, we have "inside" and
            #  "outside" x-coordinates.  The affine mapping can be specified in two
            #  ways:
            #  - using setScale() and setOffset(),
            #  - using setMapping()
            #  
            #  By default the identity mapping (scale=1, offset=0) is used.
            #  
            #  Using the value() and derivative() methods you can sample linearly
            #  interpolated values for a given x-coordinate position of the data and
            #  the derivative of the data

        # wrap-instances:
        #  LinearInterpolation := LinearInterpolation[double, double]
        LinearInterpolation() except + nogil 
        LinearInterpolation(LinearInterpolation &) except + nogil 
        ValueType value(KeyType arg_pos) except + nogil  # wrap-doc:Returns the interpolated value
        void addValue(KeyType arg_pos, ValueType arg_value) except + nogil  # wrap-doc:Performs linear resampling. The `arg_value` is split up and added to the data points around `arg_pos`
        ValueType derivative(KeyType arg_pos) except + nogil  # wrap-doc:Returns the interpolated derivative

        libcpp_vector[ValueType]  getData() except + nogil  # wrap-doc:Returns the internal random access container from which interpolated values are being sampled
        void setData(libcpp_vector[ValueType] & data) except + nogil  # wrap-doc:Assigns data to the internal random access container from which interpolated values are being sampled
        bool empty() except + nogil  # wrap-doc:Returns `true` if getData() is empty
        KeyType key2index(KeyType pos) except + nogil  # wrap-doc:The transformation from "outside" to "inside" coordinates
        KeyType index2key(KeyType pos) except + nogil  # wrap-doc:The transformation from "inside" to "outside" coordinates
        KeyType  getScale() except + nogil  # wrap-doc:"Scale" is the difference (in "outside" units) between consecutive entries in "Data"
        void setScale(KeyType & scale) except + nogil  # wrap-doc:"Scale" is the difference (in "outside" units) between consecutive entries in "Data"
        KeyType  getOffset() except + nogil  # wrap-doc:"Offset" is the point (in "outside" units) which corresponds to "Data[0]"
        void setOffset(KeyType & offset) except + nogil  # wrap-doc:"Offset" is the point (in "outside" units) which corresponds to "Data[0]"
        void setMapping(KeyType & scale, KeyType & inside, KeyType & outside) except + nogil # TODO
        void setMapping(KeyType & inside_low, KeyType & outside_low, KeyType & inside_high, KeyType & outside_high) except + nogil # TODO
        KeyType  getInsideReferencePoint() except + nogil  
        KeyType  getOutsideReferencePoint() except + nogil 
        KeyType supportMin() except + nogil 
        KeyType supportMax() except + nogil 
        LinearInterpolation(KeyType scale, KeyType offset) except + nogil 
