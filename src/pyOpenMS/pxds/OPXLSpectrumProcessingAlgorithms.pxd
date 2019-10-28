from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool
from MSSpectrum cimport *
from MSExperiment cimport *
from DataArrays cimport *
from SimpleTSGXLMS cimport SimplePeak


cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>" namespace "OpenMS":

    cdef cppclass OPXLSpectrumProcessingAlgorithms:

        OPXLSpectrumProcessingAlgorithms(OPXLSpectrumProcessingAlgorithms) nogil except +
        OPXLSpectrumProcessingAlgorithms() nogil except +

        MSSpectrum mergeAnnotatedSpectra(MSSpectrum& first_spectrum,
                                         MSSpectrum& second_spectrum) nogil except +

        MSExperiment preprocessSpectra(MSExperiment& exp,
                                        double fragment_mass_tolerance,
                                        bool fragment_mass_tolerance_unit_ppm,
                                        Size peptide_min_size,
                                        Int min_precursor_charge,
                                        Int max_precursor_charge,
                                        libcpp_vector[ size_t ]& discarded_spectra,
                                        bool deisotope,
                                        bool labeled) nogil except +

        void getSpectrumAlignmentFastCharge(libcpp_vector[ libcpp_pair[ size_t, size_t ] ]& alignment,
                                            double fragment_mass_tolerance,
                                            bool fragment_mass_tolerance_unit_ppm,
                                            const MSSpectrum& theo_spectrum,
                                            const MSSpectrum& exp_spectrum,
                                            const IntegerDataArray& theo_charges,
                                            const IntegerDataArray& exp_charges,
                                            FloatDataArray& ppm_error_array,
                                            double intensity_cutoff) nogil except +

        void getSpectrumAlignmentSimple(libcpp_vector[ libcpp_pair[ size_t, size_t ] ]& alignment,
                                            double fragment_mass_tolerance,
                                            bool fragment_mass_tolerance_unit_ppm,
                                            const libcpp_vector[ SimplePeak ]& theo_spectrum,
                                            const MSSpectrum& exp_spectrum,
                                            const IntegerDataArray& exp_charges) nogil except +
