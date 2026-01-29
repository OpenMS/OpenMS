"""Pure Python addon for MzIdentMLFile backward compatibility."""
import warnings
from . import addon
from .idxmlfile import _fill_container


@addon("MzIdentMLFile")
def load(self, filename, protein_ids=None, peptide_ids=None):
    """Load identifications from an mzIdentML file."""
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


@addon("MzIdentMLFile")
def store(self, filename, protein_ids, peptide_ids):
    """Store identifications to an mzIdentML file."""
    if isinstance(peptide_ids, list):
        warnings.warn(
            "Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5. "
            "Use PeptideIdentificationList instead.",
            DeprecationWarning,
            stacklevel=2,
        )
    self._store_internal(filename, list(protein_ids), list(peptide_ids))
