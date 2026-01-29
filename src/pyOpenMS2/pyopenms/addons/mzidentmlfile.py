"""Pure Python addon for MzIdentMLFile backward compatibility."""
from . import addon
from .idxmlfile import _load_with_compat, _store_with_compat


@addon("MzIdentMLFile")
def load(self, filename, protein_ids=None, peptide_ids=None):
    """Load identifications from an mzIdentML file."""
    return _load_with_compat(self, filename, protein_ids, peptide_ids)


@addon("MzIdentMLFile")
def store(self, filename, protein_ids, peptide_ids):
    """Store identifications to an mzIdentML file."""
    _store_with_compat(self, filename, protein_ids, peptide_ids)
