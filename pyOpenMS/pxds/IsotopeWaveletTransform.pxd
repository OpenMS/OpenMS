from libcpp cimport bool
from String cimport *
from IsotopeWavelet cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
# from LinearRegression cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeWaveletTransform[PeakT]:
        # wrap-instances:
        #   IsotopeWaveletTransform := IsotopeWaveletTransform[Peak1D]

        # IsotopeWaveletTransform() nogil except +
        IsotopeWaveletTransform(IsotopeWaveletTransform) nogil except + #wrap-ignore
        IsotopeWaveletTransform(DoubleReal min_mz, DoubleReal max_mz, UInt max_charge, Size max_scan_size, bool use_cuda, bool hr_data, String intenstype) nogil except +
        void getTransform(MSSpectrum[ PeakT ] &c_trans, MSSpectrum[ PeakT ] &c_ref, UInt c) nogil except +
        void getTransformHighRes(MSSpectrum[ PeakT ] &c_trans, MSSpectrum[ PeakT ] &c_ref, UInt c) nogil except +
        void identifyCharge(MSSpectrum[ PeakT ] &candidates, MSSpectrum[ PeakT ] &ref, UInt scan_index, UInt c, DoubleReal ampl_cutoff, bool check_PPMs) nogil except +
        void initializeScan(MSSpectrum[ PeakT ] &c_ref, UInt c) nogil except +
        void updateBoxStates(MSExperiment[ PeakT, ChromatogramPeak ] &map_, Size scan_index, UInt RT_interleave, UInt RT_votes_cutoff, Int front_bound, Int end_bound) nogil except +
        # void mergeFeatures(IsotopeWaveletTransform[ PeakT ] *later_iwt, UInt RT_interleave, UInt RT_votes_cutoff) nogil except +
        FeatureMap[ Feature ] mapSeeds2Features(MSExperiment[ PeakT, ChromatogramPeak] &map_, UInt RT_votes_cutoff) nogil except +
        ## std::multimap[ DoubleReal, Box ] getClosedBoxes() nogil except +
        ## DoubleReal getLinearInterpolation(typename MSSpectrum[ PeakT ]::const_iterator &left_iter, DoubleReal mz_pos, typename MSSpectrum[ PeakT ]::const_iterator &right_iter) nogil except +
        DoubleReal getLinearInterpolation(DoubleReal mz_a, DoubleReal intens_a, DoubleReal mz_pos, DoubleReal mz_b, DoubleReal intens_b) nogil except +
        DoubleReal getSigma() nogil except +
        void setSigma(DoubleReal sigma) nogil except +
        void computeMinSpacing(MSSpectrum[ PeakT ] &c_ref) nogil except +
        DoubleReal getMinSpacing() nogil except +
        Size getMaxScanSize() nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>" namespace "OpenMS::IsotopeWaveletTransform":
    
    cdef cppclass TransSpectrum "OpenMS::IsotopeWaveletTransform::TransSpectrum":
        TransSpectrum() nogil except +
        TransSpectrum(TransSpectrum) nogil except + #wrap-ignore
        # POINTER #  TransSpectrum(MSSpectrum[ PeakType ] * reference) nogil except +
        void destroy() nogil except +
        DoubleReal getRT() nogil except +
        DoubleReal getMZ(UInt i) nogil except +
        DoubleReal getRefIntensity(UInt i) nogil except +
        DoubleReal getTransIntensity(UInt i) nogil except +
        void setTransIntensity(UInt i, DoubleReal intens) nogil except +
        Size size() nogil except +
        # POINTER # MSSpectrum[ PeakType ] * getRefSpectrum() nogil except +
        # POINTER # MSSpectrum[ PeakType ] * getRefSpectrum() nogil except +
        # NAMESPACE # MSSpectrum[ PeakType ]::const_iterator MZBegin(DoubleReal mz) nogil except +
        # NAMESPACE # MSSpectrum[ PeakType ]::const_iterator MZEnd(DoubleReal mz) nogil except +
        # NAMESPACE # MSSpectrum[ PeakType ]::const_iterator end() nogil except +
        # NAMESPACE # MSSpectrum[ PeakType ]::const_iterator begin() nogil except +

