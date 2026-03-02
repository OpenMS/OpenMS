"""Pure Python addon for PepXMLFile backward compatibility."""
from . import addon
from .idxmlfile import _load_with_compat, _store_with_compat


@addon("PepXMLFile")
def load(self, filename, protein_ids=None, peptide_ids=None):
    """Load identifications from a pepXML file."""
    return _load_with_compat(self, filename, protein_ids, peptide_ids)


@addon("PepXMLFile")
def store(self, filename, protein_ids, peptide_ids):
    """Store identifications to a pepXML file."""
    _store_with_compat(self, filename, protein_ids, peptide_ids)
