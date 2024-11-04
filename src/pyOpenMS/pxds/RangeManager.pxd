from Types cimport *
from DPosition cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/KERNEL/RangeManager.h>" namespace "OpenMS":

    cdef cppclass RangeBase "OpenMS::RangeBase":
        # wrap-ignore
        # no-pxd-import

        RangeBase() except + nogil
        RangeBase(RangeBase &) except + nogil

        void setMin(double min)
        void setMax(double max)
        double getMin() 
        double getMax()
        void extend(double value) except + nogil  # wrap-doc:Extend the range such that it includes the given @p value
        bool contains(double value) except + nogil  # wrap-doc:Is value within [min, max]? 


    cdef cppclass RangeRT (RangeBase) :
        # wrap-inherits:
        #   RangeBase

        RangeRT() except + nogil
        RangeRT(RangeRT &) except + nogil #compiler

        void setMinRT(double min)
        void setMaxRT(double max)
        double getMinRT()
        double getMaxRT()
        void extendRT(double value) except + nogil  # wrap-doc:Extend the range such that it includes the given @p value
        bool containsRT(double value) except + nogil  # wrap-doc:Is value within [min, max]? 


    cdef cppclass RangeMZ (RangeBase) :
        # wrap-inherits:
        #   RangeBase

        RangeMZ() except + nogil
        RangeMZ(RangeMZ &) except + nogil #compiler

        void setMinMZ(double min)
        void setMaxMZ(double max)
        double getMinMZ()
        double getMaxMZ()
        void extendMZ(double value) except + nogil  # wrap-doc:Extend the range such that it includes the given @p value
        bool containsMZ(double value) except + nogil  # wrap-doc:Is value within [min, max]? 


    cdef cppclass RangeIntensity (RangeBase) :
        # wrap-inherits:
        #   RangeBase

        RangeIntensity() except + nogil
        RangeIntensity(RangeIntensity &) except + nogil

        void setMinIntensity(double min)
        void setMaxIntensity(double max)        
        double getMinIntensity() except + nogil  # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() except + nogil  # wrap-doc:Returns the maximum intensity
        void extendIntensity(double value) except + nogil  # wrap-doc:Extend the range such that it includes the given @p value
        bool containsIntensity(double value) except + nogil  # wrap-doc:Is value within [min, max]? 
 

    cdef cppclass RangeMobility (RangeBase) :
        # wrap-inherits:
        #   RangeBase

        RangeMobility() except + nogil
        RangeMobility(RangeMobility &) except + nogil

        void setMinMobility(double min)
        void setMaxMobility(double max)        
        double getMinMobility() except + nogil  # wrap-doc:Returns the minimum Mobility
        double getMaxMobility() except + nogil  # wrap-doc:Returns the maximum Mobility
        void extendMobility(double value) except + nogil  # wrap-doc:Extend the range such that it includes the given @p value
        bool containsMobility(double value) except + nogil  # wrap-doc:Is value within [min, max]? 


    cdef cppclass RangeManagerRtMzInt "OpenMS::RangeManager<RangeRT, RangeMZ, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerRtMzInt() except + nogil 
        RangeManagerRtMzInt(RangeManagerRtMzInt &) except + nogil 

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
        RangeManagerMzInt(RangeManagerMzInt &) except + nogil 

        double getMinMZ() except + nogil  # wrap-doc:Returns the minimum m/z
        double getMaxMZ() except + nogil  # wrap-doc:Returns the maximum m/z
        double getMinIntensity() except + nogil  # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() except + nogil  # wrap-doc:Returns the maximum intensity
        void clearRanges() except + nogil  # wrap-doc:Resets all range dimensions as empty

    cdef cppclass RangeManagerRtInt "OpenMS::RangeManager<RangeRT, RangeIntensity>":
        # wrap-ignore
        # no-pxd-import
        RangeManagerRtInt() except + nogil
        RangeManagerRtInt(RangeManagerRtInt &) except + nogil

        double getMinRT() except + nogil # wrap-doc:Returns the minimum RT
        double getMaxRT() except + nogil # wrap-doc:Returns the maximum RT
        double getMinIntensity() except + nogil # wrap-doc:Returns the minimum intensity
        double getMaxIntensity() except + nogil # wrap-doc:Returns the maximum intensity
        void clearRanges() except + nogil # wrap-doc:Resets all range dimensions as empty


