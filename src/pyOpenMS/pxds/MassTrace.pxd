from Types cimport *
from String cimport *
from ConvexHull2D cimport *
from Peak2D cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector


cdef extern from "<OpenMS/KERNEL/MassTrace.h>" namespace "OpenMS":

    cdef cppclass Kernel_MassTrace "OpenMS::MassTrace":

        Kernel_MassTrace() nogil except +
        Kernel_MassTrace(Kernel_MassTrace &) nogil except +
        Kernel_MassTrace(const libcpp_vector[ Peak2D ] &trace_peaks) nogil except +

        # public members
        double fwhm_mz_avg

        Size getSize() nogil except + # wrap-doc:Returns the number of peaks contained in the mass trace
        String getLabel() nogil except + # wrap-doc:Returns label of mass trace
        void setLabel(String label) nogil except + # wrap-doc:Sets label of mass trace

        double getCentroidMZ() nogil except + # wrap-doc:Returns the centroid m/z
        double getCentroidRT() nogil except + # wrap-doc:Returns the centroid RT
        double getCentroidSD() nogil except + # wrap-doc:Returns the centroid SD
        double getFWHM() nogil except + # wrap-doc:Returns FWHM
        double getTraceLength() nogil except + # wrap-doc:Returns the length of the trace (as difference in RT)
        libcpp_pair[Size,Size] getFWHMborders() nogil except + # wrap-doc:Returns FWHM boarders
        libcpp_vector[double] getSmoothedIntensities() nogil except + # wrap-doc:Returns smoothed intensities (empty if no smoothing was explicitly done beforehand!)
        double getAverageMS1CycleTime() nogil except + # wrap-doc:Returns average scan time of mass trace

        double computeSmoothedPeakArea() nogil except + # wrap-doc:Sums all non-negative (smoothed!) intensities in the mass trace
        double computePeakArea() nogil except + # wrap-doc:Sums intensities of all peaks in the mass trace
        Size findMaxByIntPeak(bool) nogil except + # wrap-doc:Returns the index of the mass trace's highest peak within the MassTrace container (based either on raw or smoothed intensities)
        Size estimateFWHM(bool) nogil except + # wrap-doc:Estimates FWHM of chromatographic peak in seconds (based on either raw or smoothed intensities)
        double computeFwhmArea() nogil except +# TODO
        double computeFwhmAreaSmooth() nogil except + # wrap-doc:Computes chromatographic peak area within the FWHM range.
        # double computeFwhmAreaRobust() nogil except +
        # double computeFwhmAreaSmoothRobust() nogil except +
        double getIntensity(bool) nogil except + # wrap-doc:Returns the intensity
        double getMaxIntensity(bool) nogil except + # wrap-doc:Returns the max intensity

        ConvexHull2D getConvexhull() nogil except + # wrap-doc:Returns the mass trace's convex hull

        void setCentroidSD(double &tmp_sd) nogil except +
        void setSmoothedIntensities(libcpp_vector[ double ] &db_vec) nogil except + # wrap-doc:Sets smoothed intensities (smoothing is done externally, e.g. by LowessSmoothing)
        void updateSmoothedMaxRT() nogil except +
        void updateWeightedMeanRT() nogil except + # wrap-doc:Compute & update centroid RT as a intensity-weighted mean of RTs
        void updateSmoothedWeightedMeanRT() nogil except +
        void updateMedianRT() nogil except + # wrap-doc:Compute & update centroid RT as median position of intensities
        void updateMedianMZ() nogil except + # wrap-doc:Compute & update centroid m/z as median of m/z values
        void updateMeanMZ() nogil except + # wrap-doc:Compute & update centroid m/z as mean of m/z values
        void updateWeightedMeanMZ() nogil except + # wrap-doc:Compute & update centroid m/z as weighted mean of m/z values
        void updateWeightedMZsd() nogil except +
            # wrap-doc:
                #   Compute & update m/z standard deviation of mass trace as weighted mean of m/z values
                #   -----
                #   Make sure to call update(Weighted)(Mean|Median)MZ() first! <br>
                #   use getCentroidSD() to get result


        void setQuantMethod(MT_QUANTMETHOD method) nogil except + # wrap-doc:Determine if area or median is used for quantification
        MT_QUANTMETHOD getQuantMethod() nogil except + # wrap-doc:Check if area or median is used for quantification
        
cdef extern from "<OpenMS/KERNEL/MassTrace.h>" namespace "OpenMS::MassTrace":

    cdef enum MT_QUANTMETHOD:
        MT_QUANT_AREA,
        MT_QUANT_MEDIAN,
        SIZE_OF_MT_QUANTMETHOD
