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

        double getCentroidMZ() nogil except +
        double getCentroidRT() nogil except +
        double getCentroidSD() nogil except +
        double getFWHM() nogil except +
        double getTraceLength() nogil except +
        libcpp_pair[Size,Size] getFWHMborders() nogil except +
        libcpp_vector[double] getSmoothedIntensities() nogil except +
        double getScanTime() nogil except +

        double computeSmoothedPeakArea() nogil except +
        double computePeakArea() nogil except +
        Size       findMaxByIntPeak(bool) nogil except +
        Size       estimateFWHM(bool) nogil except +
        double computeFwhmArea() nogil except +
        double computeFwhmAreaSmooth() nogil except +
        double computeFwhmAreaRobust() nogil except +
        double computeFwhmAreaSmoothRobust() nogil except +
        double getIntensity(bool) nogil except +
        double getMaxIntensity(bool) nogil except +

        ConvexHull2D getConvexhull() nogil except +

        void setCentroidSD(double &tmp_sd) nogil except +
        void setSmoothedIntensities(libcpp_vector[ double ] &db_vec) nogil except +
        void updateSmoothedMaxRT() nogil except +
        void updateWeightedMeanRT() nogil except +
        void updateSmoothedWeightedMeanRT() nogil except +
        void updateMedianRT() nogil except +
        void updateMedianMZ() nogil except +
        void updateMeanMZ() nogil except +
        void updateWeightedMeanMZ() nogil except +
        void updateWeightedMZsd() nogil except +

