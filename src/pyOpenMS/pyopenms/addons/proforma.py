"""
Addon methods for ProForma Peptidoform and PeptidoformIon classes.

These provide convenient instance methods that wrap the static ProForma methods,
making the API more Pythonic and compatible with legacy tests.
"""

from __future__ import annotations

from . import addon


# =============================================================================
# Peptidoform convenience methods
# =============================================================================


@addon("Peptidoform", "fromString")
@classmethod
def peptidoform_from_string(cls, proforma_string: str):
    """
    Parse a ProForma string and return a Peptidoform.

    Parameters
    ----------
    proforma_string : str
        A ProForma v2 notation string (e.g., "EM[UNIMOD:35]K").

    Returns
    -------
    Peptidoform
        The parsed peptidoform AST.

    Raises
    ------
    RuntimeError
        If the string cannot be parsed.

    Examples
    --------
    >>> pf = Peptidoform.fromString("EM[UNIMOD:35]K")
    >>> pf.toString()
    'EM[UNIMOD:35]K'
    """
    # Access ProForma.parse through the module
    import pyopenms

    return pyopenms.ProForma.parse(proforma_string)


@addon("Peptidoform", "fromAASequence")
@classmethod
def peptidoform_from_aasequence(cls, aas):
    """
    Create a Peptidoform from an AASequence.

    Parameters
    ----------
    aas : AASequence
        An OpenMS AASequence object.

    Returns
    -------
    Peptidoform
        A Peptidoform representing the same peptide.

    Examples
    --------
    >>> aas = AASequence.fromString("PEPTIDE")
    >>> pf = Peptidoform.fromAASequence(aas)
    """
    import pyopenms

    return pyopenms.ProForma.fromAASequence(aas)


@addon("Peptidoform", "fromJSON")
@classmethod
def peptidoform_from_json(cls, json_str: str):
    """
    Create a Peptidoform from a JSON string.

    Parameters
    ----------
    json_str : str
        A JSON string representing a Peptidoform.

    Returns
    -------
    Peptidoform
        The deserialized Peptidoform.
    """
    import pyopenms

    return pyopenms.ProForma.peptidoformFromJSON(json_str)


@addon("Peptidoform")
def toString(self, mode=None) -> str:
    """
    Convert to ProForma string notation.

    Parameters
    ----------
    mode : ProForma.WriteMode, optional
        LOSSLESS (default) preserves original formatting,
        CANONICAL produces normalized output.

    Returns
    -------
    str
        The ProForma string representation.
    """
    import pyopenms

    if mode is None:
        mode = pyopenms.ProForma.WriteMode.LOSSLESS
    return pyopenms.ProForma.toString(self, mode)


@addon("Peptidoform")
def toAASequence(self, policy=None):
    """
    Convert to an OpenMS AASequence.

    Parameters
    ----------
    policy : ProForma.ConversionPolicy, optional
        Conversion policy. Defaults to BEST_EFFORT.

    Returns
    -------
    AASequence
        The converted sequence.

    Raises
    ------
    RuntimeError
        If policy is FAIL_ON_LOSS and conversion cannot be lossless.
    """
    import pyopenms

    if policy is None:
        policy = pyopenms.ProForma.ConversionPolicy.BEST_EFFORT
    return pyopenms.ProForma.toAASequence(self, policy)


@addon("Peptidoform")
def isRepresentableAsAASequence(self) -> bool:
    """Check if this peptidoform can be fully represented as an AASequence."""
    import pyopenms

    return pyopenms.ProForma.isRepresentableAsAASequence(self)


@addon("Peptidoform")
def getAASequenceConversionIssues(self):
    """
    Get conversion issues for AASequence conversion.

    Returns
    -------
    list of ConversionIssue
        Issues that would arise during conversion.
    """
    import pyopenms

    return pyopenms.ProForma.getAASequenceConversionIssues(self)


@addon("Peptidoform")
def canCalculateMass(self) -> bool:
    """Check if mass can be calculated for this peptidoform."""
    import pyopenms

    return pyopenms.ProForma.canCalculateMass(self)


@addon("Peptidoform")
def getMassCalculationIssues(self):
    """
    Get issues preventing mass calculation.

    Returns
    -------
    list of ConversionIssue
        Issues preventing mass calculation.
    """
    import pyopenms

    return pyopenms.ProForma.getMassCalculationIssues(self)


@addon("Peptidoform")
def getMonoWeight(self) -> float:
    """
    Calculate monoisotopic mass in Daltons.

    Returns
    -------
    float
        Monoisotopic mass.
    """
    import pyopenms

    return pyopenms.ProForma.getMonoWeight(self)


@addon("Peptidoform")
def getMZ(self, charge: int) -> float:
    """
    Calculate m/z at a given charge state.

    Parameters
    ----------
    charge : int
        Charge state.

    Returns
    -------
    float
        The m/z value.
    """
    import pyopenms

    return pyopenms.ProForma.getMZCharge(self, charge)


@addon("Peptidoform")
def canGenerateSpectrum(self) -> bool:
    """Check if a theoretical spectrum can be generated."""
    import pyopenms

    return pyopenms.ProForma.canGenerateSpectrum(self)


@addon("Peptidoform")
def getSpectrumGenerationIssues(self):
    """
    Get issues preventing spectrum generation.

    Returns
    -------
    list of ConversionIssue
        Issues preventing spectrum generation.
    """
    import pyopenms

    return pyopenms.ProForma.getSpectrumGenerationIssues(self)


@addon("Peptidoform")
def generateSpectrum(
    self,
    min_charge: int = 1,
    max_charge: int = 2,
    ion_types: str = "by",
    add_losses: bool = False,
    add_metainfo: bool = False,
):
    """
    Generate a theoretical MS/MS spectrum.

    Parameters
    ----------
    min_charge : int
        Minimum fragment charge state (default: 1).
    max_charge : int
        Maximum fragment charge state (default: 2).
    ion_types : str
        Ion types to generate: a,b,c,x,y,z for ion series,
        M for precursor, I for immonium (default: "by").
    add_losses : bool
        Include neutral losses (default: False).
    add_metainfo : bool
        Add peak annotations (default: False).

    Returns
    -------
    MSSpectrum
        The theoretical spectrum.
    """
    import pyopenms

    return pyopenms.ProForma.generateSpectrum(
        self, min_charge, max_charge, ion_types, add_losses, add_metainfo
    )


@addon("Peptidoform")
def __repr__(self) -> str:
    """Return a detailed representation."""
    seq = self.toString()
    try:
        mass = self.getMonoWeight() if self.canCalculateMass() else None
        if mass is not None:
            return f"Peptidoform(sequence={seq!r}, mono_mass={mass:.4f})"
    except Exception:
        pass
    return f"Peptidoform(sequence={seq!r})"


@addon("Peptidoform")
def __str__(self) -> str:
    """Return the ProForma string representation."""
    return self.toString()


# =============================================================================
# PeptidoformIon convenience methods
# =============================================================================


@addon("PeptidoformIon", "fromString")
@classmethod
def peptidoform_ion_from_string(cls, proforma_string: str):
    """
    Parse a ProForma string and return a PeptidoformIon.

    Parameters
    ----------
    proforma_string : str
        A ProForma v2 notation string with optional charge (e.g., "PEPTIDE/2").

    Returns
    -------
    PeptidoformIon
        The parsed peptidoform ion AST.
    """
    import pyopenms

    return pyopenms.ProForma.parseIon(proforma_string)


@addon("PeptidoformIon", "fromJSON")
@classmethod
def peptidoform_ion_from_json(cls, json_str: str):
    """
    Create a PeptidoformIon from a JSON string.

    Parameters
    ----------
    json_str : str
        A JSON string representing a PeptidoformIon.

    Returns
    -------
    PeptidoformIon
        The deserialized PeptidoformIon.
    """
    import pyopenms

    return pyopenms.ProForma.peptidoformIonFromJSON(json_str)


@addon("PeptidoformIon")
def toString(self, mode=None) -> str:
    """
    Convert to ProForma string notation.

    Parameters
    ----------
    mode : ProForma.WriteMode, optional
        LOSSLESS (default) preserves original formatting,
        CANONICAL produces normalized output.

    Returns
    -------
    str
        The ProForma string representation.
    """
    import pyopenms

    if mode is None:
        mode = pyopenms.ProForma.WriteMode.LOSSLESS
    return pyopenms.ProForma.toStringIon(self, mode)


@addon("PeptidoformIon")
def canCalculateMass(self) -> bool:
    """Check if mass can be calculated for this peptidoform ion."""
    import pyopenms

    return pyopenms.ProForma.canCalculateMassIon(self)


@addon("PeptidoformIon")
def getMassCalculationIssues(self):
    """
    Get issues preventing mass calculation.

    Returns
    -------
    list of ConversionIssue
        Issues preventing mass calculation.
    """
    import pyopenms

    return pyopenms.ProForma.getMassCalculationIssuesIon(self)


@addon("PeptidoformIon")
def getMonoWeight(self) -> float:
    """
    Calculate monoisotopic mass in Daltons.

    Returns
    -------
    float
        Monoisotopic mass.
    """
    import pyopenms

    return pyopenms.ProForma.getMonoWeightIon(self)


@addon("PeptidoformIon")
def getMZ(self) -> float:
    """
    Calculate m/z using the embedded charge state.

    Returns
    -------
    float
        The m/z value.
    """
    import pyopenms

    return pyopenms.ProForma.getMZ(self)


@addon("PeptidoformIon")
def canGenerateSpectrum(self) -> bool:
    """Check if a theoretical spectrum can be generated."""
    import pyopenms

    return pyopenms.ProForma.canGenerateSpectrumIon(self)


@addon("PeptidoformIon")
def getSpectrumGenerationIssues(self):
    """
    Get issues preventing spectrum generation.

    Returns
    -------
    list of ConversionIssue
        Issues preventing spectrum generation.
    """
    import pyopenms

    return pyopenms.ProForma.getSpectrumGenerationIssuesIon(self)


@addon("PeptidoformIon")
def generateSpectrum(
    self,
    min_charge: int = 1,
    max_charge: int = 2,
    ion_types: str = "by",
    add_losses: bool = False,
    add_metainfo: bool = False,
):
    """
    Generate a theoretical MS/MS spectrum.

    Parameters
    ----------
    min_charge : int
        Minimum fragment charge state (default: 1).
    max_charge : int
        Maximum fragment charge state (default: 2).
    ion_types : str
        Ion types to generate (default: "by").
    add_losses : bool
        Include neutral losses (default: False).
    add_metainfo : bool
        Add peak annotations (default: False).

    Returns
    -------
    MSSpectrum
        The theoretical spectrum.
    """
    import pyopenms

    return pyopenms.ProForma.generateSpectrumIon(
        self, min_charge, max_charge, ion_types, add_losses, add_metainfo
    )


@addon("PeptidoformIon")
def __len__(self) -> int:
    """Return the number of chains in this ion."""
    return len(self.chains)


@addon("PeptidoformIon")
def __repr__(self) -> str:
    """Return a detailed representation."""
    return f"PeptidoformIon(chains={len(self.chains)})"


@addon("PeptidoformIon")
def __str__(self) -> str:
    """Return the ProForma string representation."""
    return self.toString()
