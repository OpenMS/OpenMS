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

        IsotopeWaveletTransform(double min_mz, double max_mz, UInt max_charge, Size max_scan_size, bool hr_data, String intenstype) nogil except +
        IsotopeWaveletTransform(IsotopeWaveletTransform &) nogil except + # compiler

        void getTransform(MSSpectrum &c_trans, MSSpectrum &c_ref, UInt c) nogil except +
            # wrap-doc:
                #   Computes the isotope wavelet transform of charge state `c`
                #   -----
                #   :param c_trans: The transform
                #   :param c_ref: The reference spectrum
                #   :param c: The charge state minus 1 (e.g. c=2 means charge state 3) at which you want to compute the transform

        void getTransformHighRes(MSSpectrum &c_trans, MSSpectrum &c_ref, UInt c) nogil except +
            # wrap-doc:
                #   Computes the isotope wavelet transform of charge state `c`
                #   -----
                #   :param c_trans: The transform
                #   :param c_ref: The reference spectrum
                #   :param c: The charge state minus 1 (e.g. c=2 means charge state 3) at which you want to compute the transform

        void identifyCharge(MSSpectrum &candidates, MSSpectrum &ref, UInt scan_index, UInt c, double ampl_cutoff, bool check_PPMs) nogil except +
            # wrap-doc:
                #   Given an isotope wavelet transformed spectrum 'candidates', this function assigns to every significant
                #   pattern its corresponding charge state and a score indicating the reliability of the prediction. The result of this
                #   process is stored internally. Important: Before calling this function, apply updateRanges() to the original map
                #   -----
                #   :param candidates: A isotope wavelet transformed spectrum. Entry "number i" in this vector must correspond to the
                #   charge-"(i-1)"-transform of its mass signal. (This is exactly the output of the function `getTransforms`.)
                #   :param ref: The reference scan (the untransformed raw data) corresponding to `candidates`
                #   :param c: The corresponding charge state minus 1 (e.g. c=2 means charge state 3)
                #   :param scan_index: The index of the scan (w.r.t. to some map) currently under consideration
                #   :param ampl_cutoff: The thresholding parameter. This parameter is the only (and hence a really important)
                #   parameter of the isotope wavelet transform. On the basis of `ampl_cutoff` the program tries to distinguish between
                #   noise and signal. Please note that it is not a "simple" hard thresholding parameter in the sense of drawing a virtual
                #   line in the spectrum, which is then used as a guillotine cut. Maybe you should play around a bit with this parameter to
                #   get a feeling about its range. For peptide mass fingerprints on small data sets (like single MALDI-scans e.g.), it
                #   makes sense to start `ampl_cutoff=0` or even `ampl_cutoff=-1`,
                #   indicating no thresholding at all. Note that also ampl_cutoff=0 triggers (a moderate) thresholding based on the
                #   average intensity in the wavelet transform
                #   :param check_PPMs: If enabled, the algorithm will check each monoisotopic mass candidate for its plausibility
                #   by computing the ppm difference between this mass and the averagine model

        void initializeScan(MSSpectrum &c_ref, UInt c) nogil except +# TODO
        void updateBoxStates(MSExperiment &map_, Size scan_index, UInt RT_interleave, UInt RT_votes_cutoff, Int front_bound, Int end_bound) nogil except +
            # wrap-doc:
                #   A function keeping track of currently open and closed sweep line boxes
                #   This function is used by the isotope wavelet feature finder and must be called for each processed scan
                #   -----
                #   :param map: The original map containing the data set to be analyzed
                #   :param scan_index: The index of the scan currently under consideration w.r.t. its MS map
                #   This information is necessary to sweep across the map after each scan has been evaluated
                #   :param RT_votes_cutoff: See the IsotopeWaveletFF class

        # void mergeFeatures(IsotopeWaveletTransform[ PeakT ] *later_iwt, UInt RT_interleave, UInt RT_votes_cutoff) nogil except +
        FeatureMap mapSeeds2Features(MSExperiment &map_, UInt RT_votes_cutoff) nogil except +
            # wrap-doc:
                #   Filters the candidates further more and maps the internally used data structures to the OpenMS framework
                #   -----
                #   :param map: The original map containing the data set to be analyzed
                #   :param max_charge: The maximal charge state under consideration
                #   :param RT_votes_cutoff: See the IsotopeWaveletFF class

        ## std::multimap[ double, Box ] getClosedBoxes() nogil except +
        ## double getLinearInterpolation(typename MSSpectrum::const_iterator &left_iter, double mz_pos, typename MSSpectrum::const_iterator &right_iter) nogil except +
        double getLinearInterpolation(double mz_a, double intens_a, double mz_pos, double mz_b, double intens_b) nogil except +
            # wrap-doc:
                #   Computes a linear (intensity) interpolation
                #   -----
                #   :param mz_a: The m/z value of the point left to the query
                #   :param intens_a: The intensity value of the point left to the query
                #   :param mz_pos: The query point
                #   :param mz_b: The m/z value of the point right to the query
                #   :param intens_b: The intensity value of the point left to the query

        double getSigma() nogil except +# TODO
        void setSigma(double sigma) nogil except +# TODO
        void computeMinSpacing(MSSpectrum &c_ref) nogil except +# TODO
        double getMinSpacing() nogil except +# TODO
        Size getMaxScanSize() nogil except +# TODO

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
