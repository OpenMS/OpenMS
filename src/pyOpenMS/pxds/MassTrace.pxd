from Types cimport *
from String cimport *
from ConvexHull2D cimport *
from Peak2D cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector


cdef extern from "<OpenMS/KERNEL/MassTrace.h>" namespace "OpenMS":

    cdef cppclass Kernel_MassTrace "OpenMS::MassTrace":

        Kernel_MassTrace() except + nogil 
        Kernel_MassTrace(Kernel_MassTrace &) except + nogil 
        Kernel_MassTrace(const libcpp_vector[ Peak2D ] &trace_peaks) except + nogil 

        # public members
        double fwhm_mz_avg

        Size getSize() except + nogil  # wrap-doc:Returns the number of peaks contained in the mass trace
        String getLabel() except + nogil  # wrap-doc:Returns label of mass trace
        void setLabel(String label) except + nogil  # wrap-doc:Sets label of mass trace

        double getCentroidMZ() except + nogil  # wrap-doc:Returns the centroid m/z
        double getCentroidRT() except + nogil  # wrap-doc:Returns the centroid RT
        double getCentroidSD() except + nogil  # wrap-doc:Returns the centroid SD
        double getFWHM() except + nogil  # wrap-doc:Returns FWHM
        double getTraceLength() except + nogil  # wrap-doc:Returns the length of the trace (as difference in RT)
        libcpp_pair[Size,Size] getFWHMborders() except + nogil  # wrap-doc:Returns FWHM boarders
        libcpp_vector[double] getSmoothedIntensities() except + nogil  # wrap-doc:Returns smoothed intensities (empty if no smoothing was explicitly done beforehand!)
        double getAverageMS1CycleTime() except + nogil  # wrap-doc:Returns average scan time of mass trace

        double computeSmoothedPeakArea() except + nogil  # wrap-doc:Sums all non-negative (smoothed!) intensities in the mass trace
        double computePeakArea() except + nogil  # wrap-doc:Sums intensities of all peaks in the mass trace
        double computeIntensitySum() except + nogil  # wrap-doc:Sum all peak intensities in the mass trace
        Size findMaxByIntPeak(bool) except + nogil  # wrap-doc:Returns the index of the mass trace's highest peak within the MassTrace container (based either on raw or smoothed intensities)
        Size estimateFWHM(bool) except + nogil  # wrap-doc:Estimates FWHM of chromatographic peak in seconds (based on either raw or smoothed intensities)
        double computeFwhmArea() except + nogil # TODO
        double computeFwhmAreaSmooth() except + nogil  # wrap-doc:Computes chromatographic peak area within the FWHM range.
        # double computeFwhmAreaRobust() except + nogil 
        # double computeFwhmAreaSmoothRobust() except + nogil 
        double getIntensity(bool) except + nogil  # wrap-doc:Returns the intensity
        double getMaxIntensity(bool) except + nogil  # wrap-doc:Returns the max intensity

        ConvexHull2D getConvexhull() except + nogil  # wrap-doc:Returns the mass trace's convex hull

        void setCentroidSD(double &tmp_sd) except + nogil 
        void setSmoothedIntensities(libcpp_vector[ double ] &db_vec) except + nogil  # wrap-doc:Sets smoothed intensities (smoothing is done externally, e.g. by LowessSmoothing)
        void updateSmoothedMaxRT() except + nogil 
        void updateWeightedMeanRT() except + nogil  # wrap-doc:Compute & update centroid RT as a intensity-weighted mean of RTs
        void updateSmoothedWeightedMeanRT() except + nogil 
        void updateMedianRT() except + nogil  # wrap-doc:Compute & update centroid RT as median position of intensities
        void updateMedianMZ() except + nogil  # wrap-doc:Compute & update centroid m/z as median of m/z values
        void updateMeanMZ() except + nogil  # wrap-doc:Compute & update centroid m/z as mean of m/z values
        void updateWeightedMeanMZ() except + nogil  # wrap-doc:Compute & update centroid m/z as weighted mean of m/z values
        void updateWeightedMZsd() except + nogil 
        # wrap-doc:
            #  Compute & update m/z standard deviation of mass trace as weighted mean of m/z values
            #  
            #  Make sure to call update(Weighted)(Mean|Median)MZ() first! <br>
            #  use getCentroidSD() to get result


        void setQuantMethod(MT_QUANTMETHOD method) except + nogil  # wrap-doc:Determine if area or median is used for quantification
        MT_QUANTMETHOD getQuantMethod() except + nogil  # wrap-doc:Check if area or median is used for quantification
        
cdef extern from "<OpenMS/KERNEL/MassTrace.h>" namespace "OpenMS::MassTrace":

    cdef enum MT_QUANTMETHOD:
        MT_QUANT_AREA,
        MT_QUANT_MEDIAN,
        SIZE_OF_MT_QUANTMETHOD
