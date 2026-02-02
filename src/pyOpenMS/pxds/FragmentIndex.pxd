from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool

from Types cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/FragmentIndex.h>" namespace "OpenMS::FragmentIndex":

    cdef cppclass FragmentIndex_Peptide "OpenMS::FragmentIndex::Peptide":
        # wrap-attach:
        #   FragmentIndex
        # wrap-doc:
        #  Compact descriptor of a peptide instance held by the FragmentIndex.
        #
        #  Attributes:
        #    protein_idx: Index into the FASTA entries vector identifying the source protein
        #    modification_idx_: Index into modification variants for the given unmodified subsequence
        #    sequence_: Pair of (start, length) within the source protein sequence
        #    precursor_mz_: Mono-isotopic m/z at charge 1 (M+H)+

        FragmentIndex_Peptide(UInt32 protein_idx, UInt32 modification_idx, libcpp_pair[uint16_t, uint16_t] sequence, float precursor_mz) except + nogil
        FragmentIndex_Peptide(FragmentIndex_Peptide &) except + nogil  # compiler

        UInt32 protein_idx
        UInt32 modification_idx_
        libcpp_pair[uint16_t, uint16_t] sequence_
        float precursor_mz_

    cdef cppclass FragmentIndex_SpectrumMatch "OpenMS::FragmentIndex::SpectrumMatch":
        # wrap-attach:
        #   FragmentIndex
        # wrap-doc:
        #  Match between a query peak and an entry in the database.
        #
        #  Attributes:
        #    num_matched_: Number of peak-fragment hits
        #    precursor_charge_: The precursor charge used for the search
        #    isotope_error_: The isotope error used for the search
        #    peptide_idx_: Index of the matched peptide in the fragment index

        FragmentIndex_SpectrumMatch() except + nogil  # compiler
        FragmentIndex_SpectrumMatch(FragmentIndex_SpectrumMatch &) except + nogil  # compiler

        uint32_t num_matched_
        uint16_t precursor_charge_
        int16_t isotope_error_
        size_t peptide_idx_

    cdef cppclass FragmentIndex_SpectrumMatchesTopN "OpenMS::FragmentIndex::SpectrumMatchesTopN":
        # wrap-attach:
        #   FragmentIndex
        # wrap-doc:
        #  Container for SpectrumMatch results from a query.

        FragmentIndex_SpectrumMatchesTopN() except + nogil  # compiler
        FragmentIndex_SpectrumMatchesTopN(FragmentIndex_SpectrumMatchesTopN &) except + nogil  # compiler

        libcpp_vector[FragmentIndex_SpectrumMatch] hits_

        void clear() except + nogil  # wrap-doc:Clears all hits

    cdef cppclass FragmentIndex_Hit "OpenMS::FragmentIndex::Hit":
        # wrap-attach:
        #   FragmentIndex
        # wrap-doc:
        #  A match between a single query peak and a database fragment.

        FragmentIndex_Hit(UInt32 peptide_idx, float fragment_mz) except + nogil
        FragmentIndex_Hit(FragmentIndex_Hit &) except + nogil  # compiler

        UInt32 peptide_idx
        float fragment_mz

cdef extern from "<OpenMS/ANALYSIS/ID/FragmentIndex.h>" namespace "OpenMS":

    cdef cppclass FragmentIndex(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        # wrap-doc:
        #  Fragment-index-based database for fast peptide-spectrum matching.
        #
        #  Generates from a set of FASTA entries a 2D data structure storing all
        #  theoretical b and y ion masses, organized by fragment mass and precursor mass.
        #  Supports both closed search and open search (large precursor tolerance) modes.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    from pyopenms import FragmentIndex, FASTAFile, FASTAEntry, MSSpectrum
        #
        #    # Load FASTA database
        #    entries = []
        #    FASTAFile().load("database.fasta", entries)
        #
        #    # Build index
        #    fi = FragmentIndex()
        #    params = fi.getParameters()
        #    params.setValue("precursor:mass_tolerance", 10.0)
        #    params.setValue("precursor:mass_tolerance_unit", "ppm")
        #    fi.setParameters(params)
        #    fi.build(entries)
        #
        #    # Query a spectrum
        #    sms = FragmentIndex.SpectrumMatchesTopN()
        #    fi.querySpectrum(spectrum, sms)
        #    for hit in sms.hits_:
        #      peptide = fi.getPeptides()[hit.peptide_idx_]

        FragmentIndex() except + nogil  # wrap-doc:Default constructor
        FragmentIndex(FragmentIndex &) except + nogil  # compiler

        bool isBuild() except + nogil
            # wrap-doc:Returns true if the index has been built and is ready for queries

        libcpp_vector[FragmentIndex_Peptide] getPeptides() except + nogil
            # wrap-doc:Returns a copy of all peptides in the index

        void build(libcpp_vector[FASTAEntry] & fasta_entries) except + nogil
            # wrap-doc:Builds the fragment index from FASTA entries. Must be called before querying

        void clear() except + nogil
            # wrap-doc:Clears the index and resets to unbuilt state

        void querySpectrum(MSSpectrum & spectrum, FragmentIndex_SpectrumMatchesTopN & sms) except + nogil
            # wrap-doc:
            #  Queries a complete experimental spectrum against the database.
            #  Loops over all precursor charges from min to max.
            #  :param spectrum: The experimental MS/MS spectrum to query
            #  :param sms: Output container for the top spectrum matches
