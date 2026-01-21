from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp.string cimport string as libcpp_string
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProFormaData.h>" namespace "OpenMS":

    cdef enum ProFormaWriteMode "OpenMS::ProFormaWriteMode":
        # wrap-doc:
        #  Write mode for ProForma string serialization.
        #  LOSSLESS preserves original formatting, CANONICAL produces normalized output.
        LOSSLESS
        CANONICAL

    cdef enum AASequenceConversionPolicy "OpenMS::AASequenceConversionPolicy":
        # wrap-doc:
        #  Conversion policy for transforming Peptidoform to AASequence.
        #  FAIL_ON_LOSS fails on any unrepresentable construct.
        #  DROP_UNLOCALISED drops unlocalised/labile/global modifications.
        #  BEST_EFFORT converts as much as possible, skipping unsupported.
        FAIL_ON_LOSS
        DROP_UNLOCALISED
        BEST_EFFORT

    cdef enum ConversionIssueType "OpenMS::ConversionIssueType":
        UNRESOLVED_MOD
        UNLOCALISED_MOD
        LABILE_MOD
        GLOBAL_MOD
        AMBIGUOUS_MOD
        AMBIGUOUS_REGION
        MODIFIED_RANGE
        CROSS_LINK
        MULTIPLE_CHAINS
        ALTERNATIVE_MODS
        UNSUPPORTED_FEATURE

    cdef cppclass ConversionIssue "OpenMS::ConversionIssue":
        # wrap-doc:
        #  Description of a conversion issue from Peptidoform to AASequence
        ConversionIssue() except + nogil
        ConversionIssue(ConversionIssue &) except + nogil
        ConversionIssueType type
        String description
        size_t position

    cdef enum CvDatabase "OpenMS::CvDatabase":
        UNIMOD
        MOD
        RESID
        XLMOD
        GNO

    cdef cppclass CvAccession "OpenMS::CvAccession":
        # wrap-doc:
        #  Controlled vocabulary accession for a modification (e.g., UNIMOD:35)
        CvAccession() except + nogil
        CvAccession(CvAccession &) except + nogil
        CvDatabase database
        String accession

    cdef cppclass NamedMod "OpenMS::NamedMod":
        # wrap-doc:
        #  Named modification with optional CV prefix hint
        NamedMod() except + nogil
        NamedMod(NamedMod &) except + nogil
        String name

    cdef cppclass MassDelta "OpenMS::MassDelta":
        # wrap-doc:
        #  Mass delta modification with optional source hint
        MassDelta() except + nogil
        MassDelta(MassDelta &) except + nogil
        double mass
        String original_text

    cdef cppclass FormulaTag "OpenMS::FormulaTag":
        # wrap-doc:
        #  Chemical formula modification tag
        FormulaTag() except + nogil
        FormulaTag(FormulaTag &) except + nogil
        String formula_string

    cdef cppclass InfoTag "OpenMS::InfoTag":
        # wrap-doc:
        #  Info tag for arbitrary text annotations
        InfoTag() except + nogil
        InfoTag(InfoTag &) except + nogil
        String text

    cdef cppclass Label "OpenMS::Label":
        # wrap-doc:
        #  Label for cross-links, branches, or ambiguous grouping
        Label() except + nogil
        Label(Label &) except + nogil
        String identifier

    cdef cppclass SequenceElement "OpenMS::SequenceElement":
        # wrap-doc:
        #  A single amino acid with its modifications
        SequenceElement() except + nogil
        SequenceElement(SequenceElement &) except + nogil
        char amino_acid
        libcpp_vector[Modification] modifications

    cdef cppclass Modification "OpenMS::Modification":
        # wrap-doc:
        #  A modification with one or more alternative tags
        Modification() except + nogil
        Modification(Modification &) except + nogil

    cdef cppclass UnlocalisedMod "OpenMS::UnlocalisedMod":
        # wrap-doc:
        #  Unlocalised modification with optional occurrence count
        UnlocalisedMod() except + nogil
        UnlocalisedMod(UnlocalisedMod &) except + nogil
        libcpp_vector[Modification] modifications

    cdef cppclass LabileModification "OpenMS::LabileModification":
        # wrap-doc:
        #  Labile modification that may be lost during fragmentation
        LabileModification() except + nogil
        LabileModification(LabileModification &) except + nogil
        Modification modification

    cdef cppclass IsotopeReplacement "OpenMS::IsotopeReplacement":
        # wrap-doc:
        #  Isotope replacement for stable isotope labeling
        IsotopeReplacement() except + nogil
        IsotopeReplacement(IsotopeReplacement &) except + nogil
        String isotope

    cdef cppclass Peptidoform "OpenMS::Peptidoform":
        # wrap-doc:
        #  A single peptidoform (one peptide chain) with modifications
        #  Represents a complete peptide chain including global modifications,
        #  unlocalised modifications, labile modifications, terminal modifications,
        #  and the amino acid sequence with modifications.
        #
        #  Note: The following fields are intentionally not exposed and should be
        #  accessed via the ProForma parser/writer API functions:
        #  - name: optional peptidoform name (v2.1 extension)
        #  - sequence: the amino acid sequence with modifications (complex type)
        #  - global_mods: global modifications like isotope labels (complex type)
        Peptidoform() except + nogil
        Peptidoform(Peptidoform &) except + nogil
        libcpp_vector[UnlocalisedMod] unlocalised_mods
        libcpp_vector[LabileModification] labile_mods
        libcpp_vector[Modification] n_term_mods
        libcpp_vector[Modification] c_term_mods

    cdef cppclass PeptidoformIon "OpenMS::PeptidoformIon":
        # wrap-doc:
        #  A peptidoform ion (one or more chains with optional charge)
        #
        #  Note: The following fields are intentionally not exposed and should be
        #  accessed via the ProForma parser/writer API functions:
        #  - name: optional ion name (v2.1 extension)
        #  - charge: charge state specification (complex variant type)
        PeptidoformIon() except + nogil
        PeptidoformIon(PeptidoformIon &) except + nogil
        libcpp_vector[Peptidoform] chains
