from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_utf8_string # triggers input conversion provider
from StringList cimport *
from MSSpectrum cimport *


cdef extern from "<OpenMS/CHEMISTRY/Tagger.h>" namespace "OpenMS":

    cdef cppclass Tagger:
        # wrap-doc:
                #   Constructor for Tagger
                #   -----
                #   The parameter `max_charge_` should be >= `min_charge_`
                #   Also `max_tag_length` should be >= `min_tag_length`
                #   -----
                #   :param min_tag_length: The minimal sequence tag length
                #   :param ppm: The tolerance for matching residue masses to peak delta masses
                #   :param max_tag_length: The maximal sequence tag length
                #   :param min_charge: Minimal fragment charge considered for each sequence tag
                #   :param max_charge: Maximal fragment charge considered for each sequence tag
                #   :param fixed_mods: A list of modification names. The modified residues replace the unmodified versions
                #   :param var_mods: A list of modification names. The modified residues are added as additional entries to the list of residues

        Tagger(Tagger &) nogil except + # compiler

        Tagger(size_t min_tag_length,
               double ppm,
               size_t max_tag_length,
               size_t min_charge,
               size_t max_charge,
               const StringList& fixed_mods,
               const StringList& var_mods) nogil except +

        void getTag(const libcpp_vector[ double ]& mzs,
                    libcpp_vector[ libcpp_utf8_string ]& tags) nogil except +
        # wrap-doc:
                #   Generate tags from mass vector `mzs`
                #   -----
                #   The parameter `tags` is filled with one string per sequence tag
                #   It uses the standard residues from ResidueDB including
                #   the fixed and variable modifications given to the constructor
                #   -----
                #   :param mzs: A vector of mz values, containing the mz values from a centroided fragment spectrum
                #   :param tags: The vector of tags, that is filled with this function

        void getTag(const MSSpectrum& spec,
                    libcpp_vector[ libcpp_utf8_string ]& tags) nogil except +
         # wrap-doc:
                #   Generate tags from an MSSpectrum
                #   -----
                #   The parameter `tags` is filled with one string per sequence tag
                #   It uses the standard residues from ResidueDB including
                #   the fixed and variable modifications given to the constructor
                #   -----
                #   :param spec: A centroided fragment spectrum
                #   :param tags: The vector of tags, that is filled with this function

        void setMaxCharge(size_t max_charge) nogil except +
        # wrap-doc:
                #   Change the maximal charge considered by the tagger
                #   -----
                #   Allows to change the maximal considered charge e.g. based on a spectra
                #   precursor charge without calling the constructor multiple times
                #   -----
                #   :param max_charge: The new maximal charge
