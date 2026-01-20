from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from ProFormaData cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProFormaParser.h>" namespace "OpenMS":

    cdef cppclass ProFormaParser "OpenMS::ProFormaParser":
        # wrap-doc:
        #   Recursive descent parser for ProForma v2 peptidoform notation
        #
        #   This class parses ProForma strings into an Abstract Syntax Tree (AST) representation.
        #   ProForma is a standard notation for representing peptidoforms (peptides with modifications).
        #
        #   All methods are static. Use ProFormaParser.parse() to parse a ProForma string.
        #
        #   Usage example:
        #       pf = ProFormaParser.parse("EM[UNIMOD:35]K")
        #       # pf now contains the parsed Peptidoform AST
        #       s = ProFormaParser.toString(pf)
        #       # s is "EM[UNIMOD:35]K"
        # wrap-ignore
        pass

cdef extern from "<OpenMS/CHEMISTRY/ProFormaParser.h>" namespace "OpenMS::ProFormaParser":

    # Static parse methods
    Peptidoform parse(const String& input) except + nogil
    # wrap-doc:
    #   Parse a ProForma string into a Peptidoform AST

    PeptidoformIon parseIon(const String& input) except + nogil
    # wrap-doc:
    #   Parse a ProForma string into a PeptidoformIon AST (with charge state)

    # Static toString method for Peptidoform
    String toString(const Peptidoform& pf, ProFormaWriteMode mode) except + nogil
    # wrap-doc:
    #   Convert a Peptidoform AST back to ProForma string notation with specified mode

    # Modification resolution
    void resolveModifications(Peptidoform& pf) except + nogil
    # wrap-doc:
    #   Resolve all modifications in a Peptidoform using ModificationsDB

    # AASequence conversion
    AASequence toAASequence(const Peptidoform& pf, AASequenceConversionPolicy policy) except + nogil
    # wrap-doc:
    #   Convert a Peptidoform to an OpenMS AASequence

    Peptidoform fromAASequence(const AASequence& seq) except + nogil
    # wrap-doc:
    #   Create a Peptidoform from an OpenMS AASequence

    bool isRepresentableAsAASequence(const Peptidoform& pf) except + nogil
    # wrap-doc:
    #   Check if a Peptidoform can be fully represented as an AASequence

    libcpp_vector[ConversionIssue] getAASequenceConversionIssues(const Peptidoform& pf) except + nogil
    # wrap-doc:
    #   Get a list of all issues that would arise during AASequence conversion


# JSON serialization functions (free functions in OpenMS namespace)
# These are declared in ProFormaData.h and don't require nlohmann headers
cdef extern from "<OpenMS/CHEMISTRY/ProFormaData.h>" namespace "OpenMS":

    String peptidoformToJSON "OpenMS::toJSON" (const Peptidoform& pf) except + nogil
    # wrap-doc:
    #   Convert Peptidoform to JSON string representation

    Peptidoform peptidoformFromJSON(const String& json_str) except + nogil
    # wrap-doc:
    #   Construct Peptidoform from JSON string

    String peptidoformIonToJSON "OpenMS::toJSON" (const PeptidoformIon& pfi) except + nogil
    # wrap-doc:
    #   Convert PeptidoformIon to JSON string representation

    PeptidoformIon peptidoformIonFromJSON(const String& json_str) except + nogil
    # wrap-doc:
    #   Construct PeptidoformIon from JSON string
