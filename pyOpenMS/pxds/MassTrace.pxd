from Types cimport *
from String cimport *
from ConvexHull2D cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/KERNEL/MassTrace.h>" namespace "OpenMS":

    cdef cppclass MassTrace:

        MassTrace()  nogil except +
        MassTrace(MassTrace &) nogil except + # wrap-ignore

        Size getSize()
        String getLabel()
        void setLabel(String label)

        DoubleReal getCentroidMZ()
        DoubleReal getCentroidRT()
        DoubleReal getCentroidSD()
        DoubleReal getFWHM()
        DoubleReal getTraceLength()
        libcpp_pair[Size,Size] getFWHMborders()
        libcpp_vector[DoubleReal] getSmoothedIntensities()
        DoubleReal getScanTime()

        DoubleReal computeSmoothedPeakArea()
        DoubleReal computePeakArea()
        Size       findMaxByIntPeak(bool)
        Size       estimateFWHM(bool)
        DoubleReal computeFwhmArea()
        DoubleReal computeFwhmAreaSmooth()
        DoubleReal computeFwhmAreaRobust()
        DoubleReal computeFwhmAreaSmoothRobust()
        DoubleReal getIntensity(bool)
        DoubleReal getMaxIntensity(bool)

        ConvexHull2D getConvexhull()

        void setCentroidSD(DoubleReal &tmp_sd)
        void setSmoothedIntensities(libcpp_vector[ double ] &db_vec)
        void updateSmoothedMaxRT()
        void updateWeightedMeanRT()
        void updateSmoothedWeightedMeanRT()
        void updateMedianRT()
        void updateMedianMZ()
        void updateMeanMZ()
        void updateWeightedMeanMZ()
        void updateWeightedMZsd()

