from Types cimport *
from DPosition cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/KERNEL/RangeManager.h>" namespace "OpenMS":
    
    cdef cppclass RangeManagerRtMzInt "OpenMS::RangeManager<RangeRT, RangeMZ, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerRtMzInt() nogil except +
        RangeManagerRtMzInt(RangeManagerRtMzInt  &) nogil except +
        double getMinRT() nogil except + # wrap-doc:Returns the minimum RT
        double getMaxRT() nogil except + # wrap-doc:Returns the maximum RT
        double getMinMZ() nogil except + # wrap-doc:Returns the minimum m/z
        double getMaxMZ() nogil except + # wrap-doc:Returns the maximum m/z
        double getMinIntensity() nogil except + # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() nogil except + # wrap-doc:Returns the maximum intensity
        void clearRanges() nogil except + # wrap-doc:Resets all range dimensions as empty

    cdef cppclass RangeManagerMzInt "OpenMS::RangeManager<RangeMZ, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerMzInt() nogil except +
        RangeManagerMzInt(RangeManagerMzInt  &) nogil except +
        double getMinMZ() nogil except + # wrap-doc:Returns the minimum m/z
        double getMaxMZ() nogil except + # wrap-doc:Returns the maximum m/z
        double getMinIntensity() nogil except + # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() nogil except + # wrap-doc:Returns the maximum intensity
        void clearRanges() nogil except + # wrap-doc:Resets all range dimensions as empty

    cdef cppclass RangeManagerRtInt "OpenMS::RangeManager<RangeRT, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerRtInt() nogil except +
        RangeManagerRtInt(RangeManagerRtInt  &) nogil except +
        double getMinRT() nogil except + # wrap-doc:Returns the minimum RT
        double getMaxRT() nogil except + # wrap-doc:Returns the maximum RT
        double getMinIntensity() nogil except + # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() nogil except + # wrap-doc:Returns the maximum intensity
        void clearRanges() nogil except + # wrap-doc:Resets all range dimensions as empty
