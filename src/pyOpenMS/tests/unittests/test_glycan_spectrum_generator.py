"""Glycan fragmentation bindings and site-specific retention regressions."""

import pytest
import pyopenms as oms


def composition(**counts):
    result = oms.GlycanComposition()
    result.components = list(counts.items())
    return result


def test_diagnostics_and_metadata():
    generator_type = oms.TheoreticalGlycanSpectrumGenerator
    options = generator_type.Options()
    options.add_b_ions = False
    options.add_y_ions = False
    generator = generator_type(options)
    fragments = generator.get_fragments(composition(HexNAc=2, NeuAc=1))
    ions = {f.name: f for f in fragments}
    assert ions["glycan:diagnostic:HexNAc"].get_mz() == pytest.approx(204.08665, abs=1e-5)
    assert ions["glycan:diagnostic:Neu5Ac"].get_mz() == pytest.approx(292.10269, abs=1e-5)
    assert all(f.charge == 1 for f in fragments)
    spectrum = generator_type.to_spectrum(list(reversed(fragments)))
    assert spectrum.isSorted()
    assert spectrum.size() == len(fragments)
    assert len(spectrum.getStringDataArrays()[0]) == len(fragments)
    assert len(spectrum.getIntegerDataArrays()[0]) == len(fragments)
    # get_options returns a copy, so mutations cannot bypass C++ validation.
    modified_options = generator.get_options()
    modified_options.max_charge = 0
    assert generator.get_options().max_charge == 2
    with pytest.raises(Exception):
        generator.set_options(modified_options)


def test_tree_and_formula_residues():
    tree = oms.GlycanStructure()
    root = tree.add_monosaccharide("HexNAc")
    branch = tree.add_monosaccharide("Hex", root, "beta1-4")
    tree.add_monosaccharide("Fuc", branch)
    assert tree.get_nodes()[branch].parent == root
    assert tree.get_nodes()[branch].linkage == "beta1-4"
    fragments = oms.TheoreticalGlycanSpectrumGenerator().get_fragments(tree)
    internal = [f for f in fragments if f.root_cleavage == branch and f.branch_cleavages == [2]]
    assert internal
    assert internal[0].get_mz() > 0
    tag = oms.FormulaTag()
    tag.formula_string = "C6H10O5"
    custom = oms.GlycanComposition()
    custom.components = [(tag, 1)]
    assert oms.TheoreticalGlycanSpectrumGenerator().get_fragments(custom)


def test_glycopeptide_retention():
    generator_type = oms.TheoreticalGlycanSpectrumGenerator
    generator = generator_type()
    peptide = oms.AASequence.fromString("ANST")
    glycan = composition(HexNAc=2, Hex=3)
    hcd = generator.get_glycopeptide_fragments(peptide, glycan, 1, generator_type.FragmentationMethod.HCD)
    etd = generator.get_glycopeptide_fragments(peptide, glycan, 1, generator_type.FragmentationMethod.ETD)
    assert any(f.name == "peptide:b2;glycan=HexNAc1;site=N2" for f in hcd)
    assert all(f.ion_type == generator_type.IonType.PEPTIDE for f in etd)
    assert any(f.name == "peptide:c2;glycan=Hex3HexNAc2;site=N2" for f in etd)
    assert all(f.attachment_position is None for f in etd if f.name == "peptide:c1")
    with pytest.raises(Exception):
        generator.get_glycopeptide_fragments(peptide, glycan, 4, generator_type.FragmentationMethod.HCD)


def test_custom_stub_and_limits():
    generator_type = oms.TheoreticalGlycanSpectrumGenerator
    options = generator_type.Options()
    retention = generator_type.PeptideRetention()
    retention.stripped = False
    retention.stubs = [composition(HexNAc=1, Fuc=1)]
    options.peptide_retention = {"b": retention}
    generator = generator_type(options)
    fragments = generator.get_glycopeptide_fragments(
        oms.AASequence.fromString("ANST"), composition(HexNAc=2, Fuc=1), 1,
        generator_type.FragmentationMethod.HCD,
    )
    assert any(f.name == "peptide:b2;glycan=Fuc1HexNAc1;site=N2" for f in fragments)
    options.max_fragments = 1
    generator.set_options(options)
    with pytest.raises(Exception):
        generator.get_fragments(composition(HexNAc=1))


@pytest.mark.parametrize("parents", [(5, 0, 1, 2, 3), (5, 0, 0, 0, 0), (5, 0, 0, 1, 2)])
@pytest.mark.parametrize("max_cleavages", [1, 2, 3])
def test_structural_fragments_match_exhaustive_graph_cuts(parents, max_cleavages):
    """Independent oracle: remove bonds and recover all connected components."""
    generator_type = oms.TheoreticalGlycanSpectrumGenerator
    tree = oms.GlycanStructure()
    symbols = ["HexNAc", "Hex", "Hex", "Fuc", "NeuAc"]
    # Use OpenMS's isotope masses; the oracle independently checks connectivity.
    formulas = ["C8H13NO5", "C6H10O5", "C6H10O5", "C6H10O4", "C11H17NO8"]
    masses = [oms.EmpiricalFormula(formula).getMonoWeight() for formula in formulas]
    water_mass = oms.EmpiricalFormula("H2O").getMonoWeight()
    for index, (symbol, parent) in enumerate(zip(symbols, parents)):
        tree.add_monosaccharide(symbol, None if index == 0 else parent)
    options = generator_type.Options()
    options.max_charge = 1
    options.add_diagnostic_ions = False
    options.max_cleavages = max_cleavages
    actual = generator_type(options).get_fragments(tree)
    expected = {}
    # Node 5 is the reducing end (H2O); edge i joins residue i to parents[i].
    for mask in range(1, 1 << 5):
        if mask.bit_count() > max_cleavages:
            continue
        groups = [{index} for index in range(6)]
        for index, parent in enumerate(parents):
            if mask & (1 << index):
                continue
            left = next(group for group in groups if index in group)
            right = next(group for group in groups if parent in group)
            if left is not right:
                left.update(right)
                groups.remove(right)
        for group in groups:
            reducing = 5 in group
            residues = group - {5}
            root = None if reducing else min(residues)
            branches = tuple(index for index, parent in enumerate(parents) if parent in group and index not in group)
            key = (generator_type.IonType.Y if reducing else generator_type.IonType.B, root, branches)
            expected[key] = sum(masses[index] for index in residues) + (water_mass if reducing else 0)
    observed = {(ion.ion_type, ion.root_cleavage, tuple(ion.branch_cleavages)): ion.neutral_mass for ion in actual}
    assert observed.keys() == expected.keys()
    assert len(observed) == len(actual)
    for key, mass in expected.items():
        assert observed[key] == pytest.approx(mass, abs=1e-5)


def test_terminal_sites_and_ethcd():
    generator_type = oms.TheoreticalGlycanSpectrumGenerator
    generator = generator_type()
    peptide = oms.AASequence.fromString("NST")
    glycan = composition(HexNAc=1)
    for site, prefix, suffix in [(0, "b1", "y1"), (2, "y1", "b1")]:
        ions = generator.get_glycopeptide_fragments(peptide, glycan, site, generator_type.FragmentationMethod.ETHCD)
        assert any(f.name.startswith(f"peptide:{prefix};glycan=HexNAc1;") for f in ions)
        unmodified = [f for f in ions if f.name == f"peptide:{suffix}"]
        assert unmodified and all(f.attachment_position is None for f in unmodified)
    charged = oms.FormulaTag()
    charged.formula_string = "C6H10O5"
    charged.charge = 1
    invalid = oms.GlycanComposition()
    invalid.components = [(charged, 1)]
    with pytest.raises(Exception):
        generator.get_fragments(invalid)
