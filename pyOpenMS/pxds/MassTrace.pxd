from Types cimport *
from String cimport *
from ConvexHull2D cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/KERNEL/MassTrace.h>" namespace "OpenMS":

    cdef cppclass Kernel_MassTrace "OpenMS::MassTrace":

        Kernel_MassTrace()  nogil except +
        Kernel_MassTrace(Kernel_MassTrace &) nogil except + # wrap-ignore

        Size getSize() nogil except +
        String getLabel() nogil except +
        void setLabel(String label) nogil except +

        DoubleReal getCentroidMZ() nogil except +
        DoubleReal getCentroidRT() nogil except +
        DoubleReal getCentroidSD() nogil except +
        DoubleReal getFWHM() nogil except +
        DoubleReal getTraceLength() nogil except +
        libcpp_pair[Size,Size] getFWHMborders() nogil except +
        libcpp_vector[DoubleReal] getSmoothedIntensities() nogil except +
        DoubleReal getScanTime() nogil except +

        DoubleReal computeSmoothedPeakArea() nogil except +
        DoubleReal computePeakArea() nogil except +
        Size       findMaxByIntPeak(bool) nogil except +
        Size       estimateFWHM(bool) nogil except +
        DoubleReal computeFwhmArea() nogil except +
        DoubleReal computeFwhmAreaSmooth() nogil except +
        DoubleReal computeFwhmAreaRobust() nogil except +
        DoubleReal computeFwhmAreaSmoothRobust() nogil except +
        DoubleReal getIntensity(bool) nogil except +
        DoubleReal getMaxIntensity(bool) nogil except +

        ConvexHull2D getConvexhull() nogil except +

        void setCentroidSD(DoubleReal &tmp_sd) nogil except +
        void setSmoothedIntensities(libcpp_vector[ double ] &db_vec) nogil except +
        void updateSmoothedMaxRT() nogil except +
        void updateWeightedMeanRT() nogil except +
        void updateSmoothedWeightedMeanRT() nogil except +
        void updateMedianRT() nogil except +
        void updateMedianMZ() nogil except +
        void updateMeanMZ() nogil except +
        void updateWeightedMeanMZ() nogil except +
        void updateWeightedMZsd() nogil except +

