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
        IsotopeWaveletTransform(double min_mz, double max_mz, UInt max_charge, Size max_scan_size, bool hr_data, String intenstype) nogil except +

        void getTransform(MSSpectrum &c_trans, MSSpectrum &c_ref, UInt c) nogil except +
        void getTransformHighRes(MSSpectrum &c_trans, MSSpectrum &c_ref, UInt c) nogil except +
        void identifyCharge(MSSpectrum &candidates, MSSpectrum &ref, UInt scan_index, UInt c, double ampl_cutoff, bool check_PPMs) nogil except +
        void initializeScan(MSSpectrum &c_ref, UInt c) nogil except +
        void updateBoxStates(MSExperiment &map_, Size scan_index, UInt RT_interleave, UInt RT_votes_cutoff, Int front_bound, Int end_bound) nogil except +
        # void mergeFeatures(IsotopeWaveletTransform[ PeakT ] *later_iwt, UInt RT_interleave, UInt RT_votes_cutoff) nogil except +
        FeatureMap mapSeeds2Features(MSExperiment &map_, UInt RT_votes_cutoff) nogil except +
        ## std::multimap[ double, Box ] getClosedBoxes() nogil except +
        ## double getLinearInterpolation(typename MSSpectrum::const_iterator &left_iter, double mz_pos, typename MSSpectrum::const_iterator &right_iter) nogil except +
        double getLinearInterpolation(double mz_a, double intens_a, double mz_pos, double mz_b, double intens_b) nogil except +
        double getSigma() nogil except +
        void setSigma(double sigma) nogil except +
        void computeMinSpacing(MSSpectrum &c_ref) nogil except +
        double getMinSpacing() nogil except +
        Size getMaxScanSize() nogil except +

# TODO C++ compiler errors
# pyopenms/pyopenms.cpp: error: template argument 1 is invalid
# on boost::shared_ptr<OpenMS::IsotopeWaveletTransform::TransSpectrum> inst;
# cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>" namespace "OpenMS::IsotopeWaveletTransform":

#     cdef cppclass TransSpectrum "OpenMS::IsotopeWaveletTransform::TransSpectrum":
#         TransSpectrum() nogil except +
#         TransSpectrum(TransSpectrum) nogil except + #wrap-ignore

#         # POINTER #  TransSpectrum(MSSpectrum * reference) nogil except +
#         void destroy() nogil except +
#         double getRT() nogil except +
#         double getMZ(UInt i) nogil except +
#         double getRefIntensity(UInt i) nogil except +
#         double getTransIntensity(UInt i) nogil except +
#         void setTransIntensity(UInt i, double intens) nogil except +
#         Size size() nogil except +
#         # POINTER # MSSpectrum * getRefSpectrum() nogil except +
#         # POINTER # MSSpectrum * getRefSpectrum() nogil except +
#         # NAMESPACE # MSSpectrum::const_iterator MZBegin(double mz) nogil except +
#         # NAMESPACE # MSSpectrum::const_iterator MZEnd(double mz) nogil except +
#         # NAMESPACE # MSSpectrum::const_iterator end() nogil except +
#         # NAMESPACE # MSSpectrum::const_iterator begin() nogil except +
#

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>" namespace "OpenMS::IsotopeWaveletTransform":

    cdef cppclass BoxElement "OpenMS::IsotopeWaveletTransform::BoxElement":
        BoxElement(BoxElement) nogil except + #wrap-ignore
        double mz
        UInt c
        double score
        double intens
        double ref_intens
        double RT
        UInt RT_index
        UInt MZ_begin
        UInt MZ_end

