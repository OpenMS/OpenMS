from Types cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from TheoreticalSpectrumGenerator cimport *
from SpectrumAlignment cimport *
from Peak1D cimport *
from PeptideHit cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/CHEMISTRY/SpectrumAnnotator.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAnnotator(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        SpectrumAnnotator() nogil except +
        # wrap-doc:
                #   Annotates spectra from identifications and theoretical spectra or
                #   identifications from spectra and theoretical spectra matching
                #   with various options

        SpectrumAnnotator(SpectrumAnnotator &) nogil except +

        void annotateMatches(MSSpectrum & spec, PeptideHit & ph, TheoreticalSpectrumGenerator & tg, SpectrumAlignment & sa) nogil except +
        # wrap-doc:
                #   Adds ion match annotation to the `spec` input spectrum
                #   -----
                #   :param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
                #   :param ph: A spectrum identifications to be used for the annotation, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
                #   :param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
                #   :param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance

        void addIonMatchStatistics(PeptideIdentification & pi, MSSpectrum & spec, TheoreticalSpectrumGenerator & tg, SpectrumAlignment & sa) nogil except +
        # wrap-doc:
                #   Adds ion match statistics to `pi` PeptideIdentifcation
                #   -----
                #   :param pi: A spectrum identifications to be annotated, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
                #   :param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
                #   :param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
                #   :param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
