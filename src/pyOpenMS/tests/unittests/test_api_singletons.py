"""
Tests for database singleton API completeness (Task 1).

Verifies that all database singletons from the Cython 3.5 API
are present and functional in the nanobind build.
"""
import pytest
import pyopenms


class TestCrossLinksDB:
    """CrossLinksDB is entirely missing from nanobind."""

    def test_exists(self):
        assert hasattr(pyopenms, 'CrossLinksDB')

    def test_get_instance(self):
        db = pyopenms.CrossLinksDB.getInstance()
        assert db is not None

    def test_has(self):
        db = pyopenms.CrossLinksDB.getInstance()
        assert isinstance(db.has("DSS"), bool)

    def test_get_number_of_modifications(self):
        db = pyopenms.CrossLinksDB.getInstance()
        n = db.getNumberOfModifications()
        assert isinstance(n, int)
        assert n > 0

    def test_get_modification_by_index(self):
        db = pyopenms.CrossLinksDB.getInstance()
        mod = db.getModification(0)
        assert mod is not None

    def test_get_all_search_modifications(self):
        db = pyopenms.CrossLinksDB.getInstance()
        mods = db.getAllSearchModifications()
        assert isinstance(mods, list)

    def test_search_modifications(self):
        db = pyopenms.CrossLinksDB.getInstance()
        # Should be callable (exact args may vary)
        assert callable(getattr(db, 'searchModifications', None))

    def test_search_modifications_by_diff_mono_mass(self):
        db = pyopenms.CrossLinksDB.getInstance()
        assert callable(getattr(db, 'searchModificationsByDiffMonoMass', None))

    def test_get_best_modification_by_diff_mono_mass(self):
        db = pyopenms.CrossLinksDB.getInstance()
        assert callable(getattr(db, 'getBestModificationByDiffMonoMass', None))

    def test_is_instantiated(self):
        _ = pyopenms.CrossLinksDB.getInstance()
        assert pyopenms.CrossLinksDB.isInstantiated()

    def test_read_from_obo_file_removed(self):
        """readFromOBOFile was removed in the ModificationsDB/RibonucleotideDB refactoring."""
        db = pyopenms.CrossLinksDB.getInstance()
        assert not hasattr(db, 'readFromOBOFile')


class TestModificationsDBGaps:
    """ModificationsDB exists but is missing some methods."""

    def test_search_modifications_by_diff_mono_mass(self):
        db = pyopenms.ModificationsDB.getInstance()
        assert hasattr(db, 'searchModificationsByDiffMonoMass')
        # Should return results for a common mass shift
        results = db.searchModificationsByDiffMonoMass(42.010565, 0.01)
        assert isinstance(results, (list, set))

    def test_add_modification(self):
        db = pyopenms.ModificationsDB.getInstance()
        assert hasattr(db, 'addModification')


class TestProteaseDBGaps:
    """ProteaseDB exists but is missing getEnzymeByRegEx and hasRegEx."""

    def test_get_enzyme_by_regex(self):
        db = pyopenms.ProteaseDB.getInstance()
        assert hasattr(db, 'getEnzymeByRegEx')

    def test_has_regex(self):
        db = pyopenms.ProteaseDB.getInstance()
        assert hasattr(db, 'hasRegEx')


class TestRibonucleotideDBGaps:
    """RibonucleotideDB exists but is missing getRibonucleotideAlternatives."""

    def test_get_ribonucleotide_alternatives(self):
        db = pyopenms.RibonucleotideDB.getInstance()
        assert hasattr(db, 'getRibonucleotideAlternatives')
