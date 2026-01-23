from libcpp cimport bool
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from Types cimport *
from String cimport *
from EmpiricalFormula cimport *
from AASequence cimport *
from PeptideHit cimport *

cdef extern from "<OpenMS/CHEMISTRY/MzPAF.h>" namespace "OpenMS":

    cdef enum MzPAFIonSeries "OpenMS::MzPAFIonSeries":
        # wrap-doc:
        #  Ion series types for mzPAF peak annotations.
        #  A, B, C are N-terminal ions; X, Y, Z are C-terminal ions.
        #  PRECURSOR, IMMONIUM, INTERNAL, REPORTER, FORMULA, NAMED are special ion types.
        A
        B
        C
        X
        Y
        Z
        PRECURSOR
        IMMONIUM
        INTERNAL
        REPORTER
        FORMULA
        NAMED
        UNKNOWN

    cdef enum MzPAFDeltaUnit "OpenMS::MzPAFDeltaUnit":
        # wrap-doc:
        #  Unit for mass delta values: DALTON for Da, PPM for parts per million.
        DALTON
        PPM

    cdef enum MzPAFWriteMode "OpenMS::MzPAFWriteMode":
        # wrap-doc:
        #  Write mode for mzPAF serialization.
        #  LOSSLESS preserves original formatting, CANONICAL normalizes output.
        LOSSLESS
        CANONICAL

    cdef cppclass MzPAFNeutralLoss "OpenMS::MzPAFNeutralLoss":
        # wrap-doc:
        #  Neutral loss in an mzPAF annotation (e.g., -H2O, -NH3)
        MzPAFNeutralLoss() except + nogil
        MzPAFNeutralLoss(MzPAFNeutralLoss &) except + nogil
        EmpiricalFormula formula
        String original_text
        bool operator==(MzPAFNeutralLoss) except + nogil

    cdef cppclass MzPAFMassDelta "OpenMS::MzPAFMassDelta":
        # wrap-doc:
        #  Mass delta in an mzPAF annotation (e.g., /0.001, /-1.4ppm)
        MzPAFMassDelta() except + nogil
        MzPAFMassDelta(MzPAFMassDelta &) except + nogil
        double value
        MzPAFDeltaUnit unit
        String original_text
        bool operator==(MzPAFMassDelta) except + nogil

    cdef cppclass MzPAFAnnotation "OpenMS::MzPAFAnnotation":
        # wrap-doc:
        #  A single mzPAF peak annotation.
        #
        #  Represents one annotation for a peak in mzPAF (Peak Annotation Format),
        #  the HUPO-PSI standard for fragment ion annotations.
        #
        #  Examples:
        #  - y4 - Simple y-ion at position 4
        #  - b2-H2O - b-ion with neutral loss
        #  - y4^2 - Doubly charged y-ion
        #  - y4/0.001*0.75 - With mass delta and confidence
        #  - IY - Immonium ion (tyrosine)
        #  - m3:6 - Internal fragment
        #  - r[TMT127N] - Reporter ion
        MzPAFAnnotation() except + nogil
        MzPAFAnnotation(MzPAFAnnotation &) except + nogil
        MzPAFIonSeries ion_series
        libcpp_vector[MzPAFNeutralLoss] neutral_losses
        bool isValid() except + nogil
        bool operator==(MzPAFAnnotation) except + nogil

    cdef cppclass MzPAFPeakAnnotations "OpenMS::MzPAFPeakAnnotations":
        # wrap-doc:
        #  Multiple mzPAF annotations for a single peak (comma-separated alternatives)
        MzPAFPeakAnnotations() except + nogil
        MzPAFPeakAnnotations(MzPAFPeakAnnotations &) except + nogil
        libcpp_vector[MzPAFAnnotation] annotations
        bool empty() except + nogil
        size_t size() except + nogil
        bool operator==(MzPAFPeakAnnotations) except + nogil

    cdef cppclass MzPAF "OpenMS::MzPAF":
        # wrap-doc:
        #  Parser and writer for mzPAF (Peak Annotation Format) notation.
        #
        #  mzPAF is the HUPO-PSI standard string notation for fragment ion
        #  peak annotations, used in mzSpecLib and other spectral library formats.
        #
        #  Example usage:
        #
        #  .. code-block:: python
        #
        #    # Parse a single annotation
        #    ann = MzPAF.parse("y4^2-H2O/0.001*0.75")
        #
        #    # Parse multiple comma-separated annotations
        #    anns = MzPAF.parseMultiple("b2,y4^2")
        #
        #    # Convert back to string
        #    s = MzPAF.toString(ann, MzPAFWriteMode.CANONICAL)
        #
        #    # Check if a string is mzPAF format
        #    if MzPAF.isMzPAFFormat("y4^2"):
        #        print("Valid mzPAF")
        #
        #  Dummy class to attach MzPAF namespace functions.
        MzPAF() except + nogil  # wrap-ignore
        MzPAF(MzPAF &) except + nogil  # wrap-ignore

cdef extern from "<OpenMS/CHEMISTRY/MzPAF.h>" namespace "OpenMS::MzPAF":

    # Static parse methods
    MzPAFAnnotation parse(const String& input) except + nogil  # wrap-attach:MzPAF wrap-doc:Parse an mzPAF string into a single annotation

    MzPAFPeakAnnotations parseMultiple(const String& input) except + nogil  # wrap-attach:MzPAF wrap-doc:Parse an mzPAF string with multiple comma-separated annotations

    # Static toString methods
    String toString(const MzPAFAnnotation& ann, MzPAFWriteMode mode) except + nogil  # wrap-attach:MzPAF wrap-doc:Convert an annotation to mzPAF string with specified mode

    String toStringMultiple "OpenMS::MzPAF::toString" (const MzPAFPeakAnnotations& anns, MzPAFWriteMode mode) except + nogil  # wrap-attach:MzPAF wrap-doc:Convert multiple annotations to comma-separated mzPAF string

    # PeakAnnotation integration
    PeptideHit_PeakAnnotation toPeakAnnotation(const MzPAFAnnotation& mzpaf, double mz, double intensity) except + nogil  # wrap-attach:MzPAF wrap-doc:Create a PeptideHit.PeakAnnotation from mzPAF data

    MzPAFPeakAnnotations fromPeakAnnotation(const PeptideHit_PeakAnnotation& peak_annotation) except + nogil  # wrap-attach:MzPAF wrap-doc:Parse mzPAF annotations from a PeptideHit.PeakAnnotation

    # Utilities
    bool isMzPAFFormat(const String& annotation) except + nogil  # wrap-attach:MzPAF wrap-doc:Check if a string appears to be in mzPAF format

    char ionSeriesToChar(MzPAFIonSeries series) except + nogil  # wrap-attach:MzPAF wrap-doc:Get the character for an ion series (a, b, c, x, y, z, etc.)

    # Note: charToIonSeries not exposed - output reference parameter incompatible with Cython bindings
