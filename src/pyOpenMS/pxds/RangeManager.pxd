from Types cimport *
from DPosition cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/KERNEL/RangeManager.h>" namespace "OpenMS":
    
    cdef cppclass RangeManagerRtMzInt "OpenMS::RangeManager<RangeRT, RangeMZ, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerRtMzInt() except + nogil 
        RangeManagerRtMzInt(RangeManagerRtMzInt  &) except + nogil 
        double getMinRT() except + nogil  # wrap-doc:Returns the minimum RT
        double getMaxRT() except + nogil  # wrap-doc:Returns the maximum RT
        double getMinMZ() except + nogil  # wrap-doc:Returns the minimum m/z
        double getMaxMZ() except + nogil  # wrap-doc:Returns the maximum m/z
        double getMinIntensity() except + nogil  # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() except + nogil  # wrap-doc:Returns the maximum intensity
        void clearRanges() except + nogil  # wrap-doc:Resets all range dimensions as empty

    cdef cppclass RangeManagerMzInt "OpenMS::RangeManager<RangeMZ, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerMzInt() except + nogil 
        RangeManagerMzInt(RangeManagerMzInt  &) except + nogil 
        double getMinMZ() except + nogil  # wrap-doc:Returns the minimum m/z
        double getMaxMZ() except + nogil  # wrap-doc:Returns the maximum m/z
        double getMinIntensity() except + nogil  # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() except + nogil  # wrap-doc:Returns the maximum intensity
        void clearRanges() except + nogil  # wrap-doc:Resets all range dimensions as empty

    cdef cppclass RangeManagerRtInt "OpenMS::RangeManager<RangeRT, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerRtInt() except + nogil 
        RangeManagerRtInt(RangeManagerRtInt  &) except + nogil 
        double getMinRT() except + nogil  # wrap-doc:Returns the minimum RT
        double getMaxRT() except + nogil  # wrap-doc:Returns the maximum RT
        double getMinIntensity() except + nogil  # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() except + nogil  # wrap-doc:Returns the maximum intensity
        void clearRanges() except + nogil  # wrap-doc:Resets all range dimensions as empty
