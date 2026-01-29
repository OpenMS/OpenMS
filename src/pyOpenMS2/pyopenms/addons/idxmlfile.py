"""Pure Python addon for IdXMLFile backward compatibility.

Also provides shared helpers for MzIdentMLFile and PepXMLFile addons.
"""
import warnings
from . import addon


def _fill_container(container, items):
    """Fill a list or PeptideIdentificationList with items."""
    if isinstance(container, list):
        container[:] = list(items)
    else:
        # PeptideIdentificationList or similar C++ container
        container.clear()
        for item in items:
            container.push_back(item)


def _load_with_compat(self, filename, protein_ids=None, peptide_ids=None, stacklevel=3):
    """Shared load logic with backward-compatible list filling."""
    result = self._load_internal(filename)

    if protein_ids is None and peptide_ids is None:
        return result

    if isinstance(peptide_ids, list):
        warnings.warn(
            "Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5. "
            "Use PeptideIdentificationList instead.",
            DeprecationWarning,
            stacklevel=stacklevel,
        )

    proteins, peptides = result

    if protein_ids is not None:
        _fill_container(protein_ids, proteins)
    if peptide_ids is not None:
        _fill_container(peptide_ids, peptides)

    return result


def _store_with_compat(self, filename, protein_ids, peptide_ids, stacklevel=3):
    """Shared store logic with deprecation warning for list arguments."""
    if isinstance(peptide_ids, list):
        warnings.warn(
            "Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5. "
            "Use PeptideIdentificationList instead.",
            DeprecationWarning,
            stacklevel=stacklevel,
        )
    self._store_internal(filename, list(protein_ids), list(peptide_ids))


@addon("IdXMLFile")
def load(self, filename, protein_ids=None, peptide_ids=None):
    """Load identifications from an idXML file.

    Supports both the new API (returns tuple) and the old API
    (fills provided lists in place).
    """
    return _load_with_compat(self, filename, protein_ids, peptide_ids)


@addon("IdXMLFile")
def store(self, filename, protein_ids, peptide_ids, document_id=""):
    """Store identifications to an idXML file."""
    _store_with_compat(self, filename, protein_ids, peptide_ids)
