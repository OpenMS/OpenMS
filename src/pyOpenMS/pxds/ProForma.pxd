from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp.string cimport string as libcpp_string
from Types cimport *
from String cimport *
from AASequence cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProForma.h>" namespace "OpenMS":

    # =========================================================================
    # Enums
    # =========================================================================

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

    cdef enum CvDatabase "OpenMS::CvDatabase":
        UNIMOD
        MOD
        RESID
        XLMOD
        GNO

    # =========================================================================
    # Data Structures
    # =========================================================================

    cdef cppclass ConversionIssue "OpenMS::ConversionIssue":
        # wrap-doc:
        #  Description of a conversion issue from Peptidoform to AASequence
        ConversionIssue() except + nogil
        ConversionIssue(ConversionIssue &) except + nogil
        ConversionIssueType type
        String description
        size_t position

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

    # =========================================================================
    # ProFormaParser Class
    # =========================================================================

    cdef cppclass ProFormaParser "OpenMS::ProFormaParser":
        # wrap-doc:
        #  Recursive descent parser for ProForma v2 peptidoform notation
        #
        #  This class parses ProForma strings into an Abstract Syntax Tree (AST) representation.
        #  ProForma is a standard notation for representing peptidoforms (peptides with modifications).
        #
        #  All methods are static. Use ProFormaParser.parse() to parse a ProForma string.
        #
        #  Usage example:
        #
        #  .. code-block:: python
        #
        #    pf = ProFormaParser.parse("EM[UNIMOD:35]K")
        #    # pf now contains the parsed Peptidoform AST
        #    s = ProFormaParser.toString(pf, ProFormaWriteMode.LOSSLESS)
        #    # s is "EM[UNIMOD:35]K"
        #
        #  Dummy class to attach ProFormaParser namespace functions.
        #  This class should not be instantiated directly.
        ProFormaParser() except + nogil  # wrap-ignore
        ProFormaParser(ProFormaParser &) except + nogil  # wrap-ignore

# =============================================================================
# ProFormaParser Static Methods
# =============================================================================

cdef extern from "<OpenMS/CHEMISTRY/ProForma.h>" namespace "OpenMS::ProFormaParser":

    # Static parse methods
    Peptidoform parse(const String& input) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Parse a ProForma string into a Peptidoform AST

    PeptidoformIon parseIon(const String& input) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Parse a ProForma string into a PeptidoformIon AST (with charge state)

    # Static toString methods
    String toString(const Peptidoform& pf, ProFormaWriteMode mode) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Convert a Peptidoform AST back to ProForma string notation with specified mode

    String toStringIon "OpenMS::ProFormaParser::toString" (const PeptidoformIon& pfi, ProFormaWriteMode mode) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Convert a PeptidoformIon AST back to ProForma string notation with specified mode

    # Modification resolution
    void resolveModifications(Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Resolve all modifications in a Peptidoform using ModificationsDB

    # AASequence conversion
    AASequence toAASequence(const Peptidoform& pf, AASequenceConversionPolicy policy) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Convert a Peptidoform to an OpenMS AASequence

    Peptidoform fromAASequence(const AASequence& seq) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Create a Peptidoform from an OpenMS AASequence

    bool isRepresentableAsAASequence(const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Check if a Peptidoform can be fully represented as an AASequence

    libcpp_vector[ConversionIssue] getAASequenceConversionIssues(const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Get a list of all issues that would arise during AASequence conversion

    # Mass calculation methods
    bool canCalculateMass(const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Check if mass can be calculated for a Peptidoform

    bool canCalculateMassIon "OpenMS::ProFormaParser::canCalculateMass" (const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Check if mass can be calculated for a PeptidoformIon

    libcpp_vector[ConversionIssue] getMassCalculationIssues(const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Get issues preventing mass calculation for a Peptidoform

    libcpp_vector[ConversionIssue] getMassCalculationIssuesIon "OpenMS::ProFormaParser::getMassCalculationIssues" (const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Get issues preventing mass calculation for a PeptidoformIon

    double getMonoWeight(const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Calculate monoisotopic mass of a Peptidoform in Daltons

    double getMonoWeightIon "OpenMS::ProFormaParser::getMonoWeight" (const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Calculate monoisotopic mass of a PeptidoformIon in Daltons

    double getMZ(const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Calculate m/z for a PeptidoformIon using its charge state

    double getMZCharge "OpenMS::ProFormaParser::getMZ" (const Peptidoform& pf, int charge) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Calculate m/z for a Peptidoform at a given charge state

    # Theoretical spectrum generation
    bool canGenerateSpectrum(const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Check if a theoretical spectrum can be generated for a Peptidoform

    bool canGenerateSpectrumIon "OpenMS::ProFormaParser::canGenerateSpectrum" (const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Check if a theoretical spectrum can be generated for a PeptidoformIon

    libcpp_vector[ConversionIssue] getSpectrumGenerationIssues(const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Get issues preventing spectrum generation for a Peptidoform

    libcpp_vector[ConversionIssue] getSpectrumGenerationIssuesIon "OpenMS::ProFormaParser::getSpectrumGenerationIssues" (const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Get issues preventing spectrum generation for a PeptidoformIon

    MSSpectrum generateSpectrum(const Peptidoform& pf, int min_charge, int max_charge, const libcpp_string& ion_types, bool add_losses, bool add_metainfo) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Generate theoretical MS/MS spectrum for a Peptidoform. ion_types uses chars a,b,c,x,y,z for ion series, M for precursor, I for immonium (e.g. "by" or "abyM")

    MSSpectrum generateSpectrumIon "OpenMS::ProFormaParser::generateSpectrum" (const PeptidoformIon& pfi, int min_charge, int max_charge, const libcpp_string& ion_types, bool add_losses, bool add_metainfo) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Generate theoretical MS/MS spectrum for a PeptidoformIon (supports cross-linked peptides). ion_types uses chars a,b,c,x,y,z for ion series, M for precursor, I for immonium

# =============================================================================
# JSON Serialization Functions
# =============================================================================

cdef extern from "<OpenMS/CHEMISTRY/ProForma.h>" namespace "OpenMS":

    String peptidoformToJSON "OpenMS::toJSON" (const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Convert Peptidoform to JSON string representation

    Peptidoform peptidoformFromJSON(const String& json_str) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Construct Peptidoform from JSON string

    String peptidoformIonToJSON "OpenMS::toJSON" (const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Convert PeptidoformIon to JSON string representation

    PeptidoformIon peptidoformIonFromJSON(const String& json_str) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Construct PeptidoformIon from JSON string
