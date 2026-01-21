from libcpp cimport bool
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from ProFormaData cimport *
from AASequence cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProFormaParser.h>" namespace "OpenMS":

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

cdef extern from "<OpenMS/CHEMISTRY/ProFormaParser.h>" namespace "OpenMS::ProFormaParser":

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


# JSON serialization functions (free functions in OpenMS namespace)
# These are declared in ProFormaData.h and don't require nlohmann headers
cdef extern from "<OpenMS/CHEMISTRY/ProFormaData.h>" namespace "OpenMS":

    String peptidoformToJSON "OpenMS::toJSON" (const Peptidoform& pf) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Convert Peptidoform to JSON string representation

    Peptidoform peptidoformFromJSON(const String& json_str) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Construct Peptidoform from JSON string

    String peptidoformIonToJSON "OpenMS::toJSON" (const PeptidoformIon& pfi) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Convert PeptidoformIon to JSON string representation

    PeptidoformIon peptidoformIonFromJSON(const String& json_str) except + nogil  # wrap-attach:ProFormaParser wrap-doc:Construct PeptidoformIon from JSON string
