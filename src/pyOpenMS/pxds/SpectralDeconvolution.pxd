from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from MSSpectrum cimport *
from DefaultParamHandler cimport *
from DeconvolvedSpectrum cimport *
from FLASHHelperClasses cimport *
from PeakGroup cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>" namespace "OpenMS":

    cdef cppclass SpectralDeconvolution(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        # wrap-doc:
        #  Spectral deconvolution algorithm for top-down MS.
        #  From MSSpectrum, this class outputs DeconvolvedSpectrum.
        #  Deconvolution takes three steps:
        #    i) decharging and select candidate masses - speed up via binning
        #    ii) collecting isotopes from the candidate masses and deisotoping
        #    iii) scoring and filter out low scoring masses

        # Constructors
        SpectralDeconvolution() except + nogil
        SpectralDeconvolution(SpectralDeconvolution &) except + nogil

        # Main deconvolution method
        void performSpectrumDeconvolution(MSSpectrum & spec, int scan_number, PeakGroup & precursor_peak_group) except + nogil
        # wrap-doc:
        #  Main deconvolution function that generates the deconvolved spectrum.
        #  :param spec: The original spectrum
        #  :param scan_number: Scan number from input spectrum
        #  :param precursor_peak_group: Precursor peak group (for MS2+)

        # Result access
        DeconvolvedSpectrum getDeconvolvedSpectrum() except + nogil
        # wrap-doc:Return the deconvolved spectrum after performSpectrumDeconvolution is called

        # Averagine access
        PrecalAveragine getAveragine() except + nogil
        # wrap-doc:Get calculated averagine. Call after calculateAveragine is called.

        void setAveragine(PrecalAveragine & avg) except + nogil
        # wrap-doc:Set the precalculated averagine

        # Configuration
        void setTargetMasses(libcpp_vector[double] & masses, bool exclude) except + nogil
        # wrap-doc:Set targeted or excluded masses for targeted deconvolution. Masses are targeted or excluded in all ms levels.

        void calculateAveragine(bool use_RNA_averagine) except + nogil
        # wrap-doc:Precalculate averagine (for predefined mass bins) to speed up averagine generation

        void setToleranceEstimation() except + nogil
        # wrap-doc:When estimating tolerance, set max_mass_dalton_tolerance to a large value

        void setTargetDecoyType(TargetDecoyType target_decoy_type, DeconvolvedSpectrum & target_dspec) except + nogil
        # wrap-doc:Set target decoy type for the SpectralDeconvolution run


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>" namespace "OpenMS::SpectralDeconvolution":

    int getNominalMass(double mass) except + nogil  # wrap-attach:SpectralDeconvolution
    # wrap-doc:Convert double mass to nominal mass (integer)

    float getCosine(libcpp_vector[float] & a, int a_start, int a_end, IsotopeDistribution & b, int offset, int min_iso_len) except + nogil  # wrap-attach:SpectralDeconvolution
    # wrap-doc:Calculate cosine between two vectors with optimization parameters

    float getIsotopeCosineAndIsoOffset(double mono_mass, libcpp_vector[float] & per_isotope_intensities, int & offset, PrecalAveragine & avg, int iso_int_shift, int window_width, libcpp_vector[double] & excluded_masses) except + nogil  # wrap-attach:SpectralDeconvolution
    # wrap-doc:Examine intensity distribution over isotope indices and determine most plausible isotope index

    int min_iso_size  # wrap-attach:SpectralDeconvolution
    # wrap-doc:Minimum isotopologue count in a peak group (=2)
