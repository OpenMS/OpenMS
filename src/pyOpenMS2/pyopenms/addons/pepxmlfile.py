"""Pure Python addon for PepXMLFile backward compatibility."""
import warnings
from . import addon
from .idxmlfile import _fill_container


@addon("PepXMLFile")
def load(self, filename, protein_ids=None, peptide_ids=None):
    """Load identifications from a pepXML file."""
    if protein_ids is None and peptide_ids is None:
        return self._load_internal(filename)

    if isinstance(peptide_ids, list):
        warnings.warn(
            "Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5. "
            "Use PeptideIdentificationList instead.",
            DeprecationWarning,
            stacklevel=2,
        )

    result = self._load_internal(filename)
    proteins, peptides = result

    if protein_ids is not None:
        _fill_container(protein_ids, proteins)
    if peptide_ids is not None:
        _fill_container(peptide_ids, peptides)


@addon("PepXMLFile")
def store(self, filename, protein_ids, peptide_ids):
    """Store identifications to a pepXML file."""
    if isinstance(peptide_ids, list):
        warnings.warn(
            "Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5. "
            "Use PeptideIdentificationList instead.",
            DeprecationWarning,
            stacklevel=2,
        )
    self._store_internal(filename, list(protein_ids), list(peptide_ids))
