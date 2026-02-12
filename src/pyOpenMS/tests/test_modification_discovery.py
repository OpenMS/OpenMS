"""
End-to-end test for synthetic modification discovery using PeptideSearchEngineFIAlgorithm.

Strategy: Generate theoretical MS/MS fragment spectra from unmodified+Carbamidomethyl
peptides, but shift precursor m/z by known modification masses. This mimics real open
search data where fragment ions match the unmodified backbone and the precursor mass
shift reveals the modification.

The search engine should discover all injected modifications via delta mass analysis.
"""
import random

import pytest

from pyopenms import (
    AASequence,
    FASTAEntry,
    ModifiedPeptideGenerator,
    MSExperiment,
    MSSpectrum,
    Param,
    Peak1D,
    PeptideIdentificationList,
    PeptideSearchEngineFIAlgorithm,
    Precursor,
    ProteaseDigestion,
    ProteinIdentification,
    TheoreticalSpectrumGenerator,
)


# Synthetic protein sequences (designed with good amino acid diversity)
PROTEIN_DB = [
    ("P01", "Protein01",
     "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
     "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"
     "GITNSEYRQWDLKAPMFHCVSITGNREYWDKLMPAHFQCSTVINEYRWDLK"
     "APMHSCFTGQNVIREYWDKLMSPAHCFQNTSGIVREYWDKLHMPASCFQGN"),
    ("P02", "Protein02",
     "MKAILNHVGSTFREDWQCPYLKMISGDTFNHRVAWQECPLKYMTGISNHFR"
     "DVEWAQCPLKTMIYGSNHFRDVEWAQCPKLIMTGSYNHFRDVEWAQCKPLIM"
     "TGSYHFNRDVEWAQCPKLITMGSYHFNRDVEWAQCPKLITMGSYHNFRDVEW"
     "AQCPKLITMGSYHFNRDVEWAQCPKLITMGSYHFNRDVEWAQCPKLITMGSY"),
    ("P03", "Protein03",
     "MGHYIKLTPNRESWDVAFQCKMHGYILKTPNRESWDVAFQCKMHGLIYTKP"
     "NRESWDVAFQCKHMGIYLTKPNRESWDVAFQCKMHGIYLKTNPRESWDVAFQ"
     "CKMHGYLIKTPNRESWDVAFQCKMHGIYLTKPNRESWDVAFQCKMHGIYLTK"
     "PNRESWDVAFQCKMHGIYLTKPNRESWDVAFQCKMHGIYLTKPNRESWDVAF"),
    ("P04", "Protein04",
     "MSVDNKTHFRGECAWYPILQMSDKTNHFRGEVAWCYQPILKMSDETKNHFRG"
     "VAWCEQYPILKMSDTKENHFRGVAWCEYQPILKMSDETKHNFRGVACWEYQPI"
     "LKMSDETKNHFRGVAWCEYQPILKMSDETKNHFRGVAWCEYQPILKMSDETNK"
     "HFRGVAWCEYQPILKMSDETKNHFRGVAWCEYQPILKMSDETKNHFRGVAWCE"),
    ("P05", "Protein05",
     "MAKLFGYNRSTECWDIPHQVMKALYFGNRSTWECDIPHQVKMALGFYNRSTWE"
     "CDIPHQVKMALFGYNRSTEWCDIPHQVKMAFGLYNRSTWECDIPHQVKMALFY"
     "GNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQ"
     "VKMALFGYNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQVKMALFGYNRSTE"),
    ("P06", "Protein06",
     "MTGYLSKFHERNDICWPAQVMTGLYSKHFERNDICWPAQVMTGLYSKFHRNDE"
     "ICWPAQVKMTGYLSKFHERNDICWPAQVKMTGLYSKHFERNDICWPAQVKMTG"
     "LYKSFHERNDICWPAQVKMTGLYSKHFERNDICWPAQVKMTGLYSKFHERNDIC"
     "WPAQVKMTGLYSKFHERNDICWPAQVKMTGLYSKFHERNDICWPAQVKMTGLY"),
    ("P07", "Protein07",
     "MDIKHWNRSYPLCTGEFAQVKMDIHKWNRSYPLCTGFAEQVKMDIKHWNRSYP"
     "LCTGFEAQVKMDIHKWNRSYPLCTGEFAQVKMDIKHWNRSYPLCTGFEAQVKM"
     "DIHKWNRSYPLCTGEFAQVKMDIHKWNRSYPLCTGFEAQVKMDIHKWNRSYPL"
     "CTGFEAQVKMDIHKWNRSYPLCTGFEAQVKMDIHKWNRSYPLCTGFEAQVKMD"),
    ("P08", "Protein08",
     "MEYKFADLHGSNTCRWQPIVKMEYFKADLHGSNTCRWQPVIKMEYFKADLHGS"
     "NTCRWQPIVKMEYFKADLHGSNTRWCQPIVKMEYFDKALGHSNTRWCQPIVKM"
     "EYFKADLGHSNTCRWQPIVKMEYFKADLHGSNTCRWQPIVKMEYFKADLHGSN"
     "TCRWQPIVKMEYFKADLHGSNTCRWQPIVKMEYFKADLHGSNTCRWQPIVKME"),
]


# Modifications to inject: (name, mass_shift_Da)
TEST_MODIFICATIONS = [
    ("Oxidation (M)",            15.9949),
    ("Phospho (S)",              79.9663),
    ("Phospho (T)",              79.9663),
    ("Acetyl (K)",               42.0106),
    ("Deamidated (N)",            0.9840),
    ("Methyl (R)",               14.0157),
    ("Dimethyl (K)",             28.0314),
    ("Carbamyl (K)",             43.0058),
    ("Dioxidation (M)",          31.9898),
    ("HexNAc (S)",              203.0794),
    ("Acetyl (Protein N-term)",  42.0106),
    ("Formyl (Protein N-term)",  27.9949),
    ("Unknown (artificial)",    123.4560),  # artificial mass not matching any known modification
]


def _build_fasta_db():
    """Create FASTA entries from protein database."""
    entries = []
    for ident, desc, seq in PROTEIN_DB:
        entry = FASTAEntry()
        entry.identifier = ident
        entry.description = desc
        entry.sequence = seq
        entries.append(entry)
    return entries


def _digest_proteins(fasta_entries):
    """Digest proteins with Trypsin and apply Carbamidomethyl(C)."""
    digester = ProteaseDigestion()
    digester.setEnzyme("Trypsin")
    digester.setMissedCleavages(1)

    fixed_mods = ModifiedPeptideGenerator.getModifications([b"Carbamidomethyl (C)"])

    all_peptides = []
    for entry in fasta_entries:
        protein = AASequence.fromString(entry.sequence)
        peptides = []
        digester.digest(protein, peptides, 7, 40)
        for pep in peptides:
            ModifiedPeptideGenerator.applyFixedModifications(fixed_mods, pep)
            all_peptides.append(pep)
    return all_peptides


def _generate_spectra(peptides, modifications, target_per_mod=300, seed=42):
    """Generate MS/MS spectra with shifted precursor m/z."""
    rng = random.Random(seed)

    tsg = TheoreticalSpectrumGenerator()
    tsg_param = tsg.getParameters()
    tsg_param.setValue("add_first_prefix_ion", "true")
    tsg_param.setValue("add_metainfo", "true")
    tsg.setParameters(tsg_param)

    exp = MSExperiment()
    rt = 100.0
    spec_count = 0

    for mod_name, mod_mass in modifications:
        indices = list(range(len(peptides)))
        rng.shuffle(indices)

        created = 0
        for idx in indices:
            if created >= target_per_mod:
                break
            pep = peptides[idx]
            if pep.size() < 8:
                continue

            charge = rng.randint(2, 4)

            # Generate theoretical spectrum from unmodified+Carbamidomethyl peptide
            spec = MSSpectrum()
            tsg.getSpectrum(spec, pep, 1, min(charge - 1, 2))
            spec.sortByPosition()
            if spec.size() < 10:
                continue

            spec.setMSLevel(2)
            spec.setRT(rt)
            rt += 0.1

            # Shift precursor m/z by modification mass
            unmod_mz = pep.getMZ(charge)
            shifted_mz = unmod_mz + mod_mass / charge

            prec = Precursor()
            prec.setMZ(shifted_mz)
            prec.setCharge(charge)
            spec.setPrecursors([prec])
            spec.setNativeID(f"spectrum={spec_count}")

            exp.addSpectrum(spec)
            created += 1
            spec_count += 1

    return exp


def test_synthetic_modification_discovery():
    """End-to-end test: generate modified peptides, create spectra, run open search, verify discovery."""

    # 1. Build database, digest, generate spectra
    fasta_db = _build_fasta_db()
    peptides = _digest_proteins(fasta_db)
    assert len(peptides) > 50, f"Expected >50 peptides, got {len(peptides)}"

    spectra = _generate_spectra(peptides, TEST_MODIFICATIONS, target_per_mod=300)
    assert spectra.size() > 1000, f"Expected >1000 spectra, got {spectra.size()}"

    # 2. Configure open search
    algo = PeptideSearchEngineFIAlgorithm()
    p = algo.getParameters()
    p.setValue("precursor:mass_tolerance", 500.0)
    p.setValue("precursor:mass_tolerance_unit", "Da")
    p.setValue("fragment:mass_tolerance", 20.0)
    p.setValue("fragment:mass_tolerance_unit", "ppm")
    p.setValue("modifications:fixed", ["Carbamidomethyl (C)"])
    p.setValue("modifications:variable", [])
    p.setValue("decoys", "false")
    p.setValue("peptide:min_size", 7)
    p.setValue("peptide:max_size", 40)
    p.setValue("peptide:missed_cleavages", 1)
    p.setValue("report:top_hits", 1)
    algo.setParameters(p)

    # 3. Run search with modification analysis
    result = algo.searchWithModificationAnalysis(spectra, fasta_db, "")

    # 4. Basic checks
    ExitCodes = PeptideSearchEngineFIAlgorithm.PeptideSearchEngineFIAlgorithm_ExitCodes
    assert result.exit_code == ExitCodes.EXECUTION_OK
    assert result.is_open_search is True
    assert result.peptide_ids.size() > 500

    # 5. Check delta mass bins for expected modifications
    dm_entries = result.modification_analysis.delta_mass_stats.entries
    assert len(dm_entries) > 0, "No delta mass entries found"
    assert result.modification_analysis.delta_mass_stats.total_psms > 0

    # Build lookup of delta masses
    dm_lookup = {e.delta_mass: e.count for e in dm_entries}

    # Expected delta masses that must be found (mass, label)
    must_find = [
        (15.9949, "Oxidation"),
        (79.9663, "Phospho"),
        (42.0106, "Acetyl"),
        (0.9840,  "Deamidated"),
        (14.0157, "Methyl"),
        (28.0314, "Dimethyl"),
        (43.0058, "Carbamyl"),
        (31.9898, "Dioxidation"),
        (203.0794, "HexNAc"),
        (27.9949, "Formyl"),
    ]

    for expected_mass, label in must_find:
        found = False
        count = 0
        for dm, cnt in dm_lookup.items():
            if abs(dm - expected_mass) < 0.05:
                found = True
                count = cnt
                break
        assert found, f"Missing delta mass bin for {label} (~{expected_mass} Da)"
        assert count >= 50, f"Low count for {label}: {count}"

    # 6. Verify that the artificial unknown modification (123.456 Da) is present
    #    in the delta mass histogram but NOT mapped to a known modification.
    found_unknown = False
    for e in dm_entries:
        if abs(e.delta_mass - 123.456) < 0.05:
            found_unknown = True
            assert e.count >= 50, f"Low count for unknown mod: {e.count}"
            assert not e.is_known_modification, "Unknown mod should not be mapped as known"
            break
    assert found_unknown, "Missing delta mass bin for artificial unknown modification (~123.456 Da)"


def test_closed_search_baseline():
    """Closed search: verify standard search works with in-memory overload."""
    fasta_db = _build_fasta_db()[:1]  # single protein

    tsg = TheoreticalSpectrumGenerator()
    tsg_param = tsg.getParameters()
    tsg_param.setValue("add_first_prefix_ion", "true")
    tsg_param.setValue("add_metainfo", "true")
    tsg.setParameters(tsg_param)

    test_seqs = [
        "VLGFHQR",
        "M(Oxidation)PNASTIC(Carbamidomethyl)YWDLK",
    ]

    exp = MSExperiment()
    rt = 100.0
    for seq_str in test_seqs:
        seq = AASequence.fromString(seq_str)
        charge = 2
        spec = MSSpectrum()
        tsg.getSpectrum(spec, seq, 1, 1)
        spec.sortByPosition()
        spec.setMSLevel(2)
        spec.setRT(rt)
        rt += 1.0

        prec = Precursor()
        prec.setMZ(seq.getMZ(charge))
        prec.setCharge(charge)
        spec.setPrecursors([prec])
        spec.setNativeID(f"spectrum={exp.size()}")
        exp.addSpectrum(spec)

    algo = PeptideSearchEngineFIAlgorithm()
    p = algo.getParameters()
    p.setValue("precursor:mass_tolerance", 10.0)
    p.setValue("precursor:mass_tolerance_unit", "ppm")
    p.setValue("fragment:mass_tolerance", 20.0)
    p.setValue("fragment:mass_tolerance_unit", "ppm")
    p.setValue("modifications:fixed", ["Carbamidomethyl (C)"])
    p.setValue("modifications:variable", ["Oxidation (M)"])
    p.setValue("decoys", "false")
    p.setValue("peptide:min_size", 7)
    p.setValue("peptide:max_size", 40)
    p.setValue("peptide:missed_cleavages", 1)
    algo.setParameters(p)

    prot_ids = []
    pep_ids = PeptideIdentificationList()
    ec = algo.search(exp, fasta_db, prot_ids, pep_ids)

    ExitCodes = PeptideSearchEngineFIAlgorithm.PeptideSearchEngineFIAlgorithm_ExitCodes
    assert ec == ExitCodes.EXECUTION_OK
    assert pep_ids.size() > 0
    assert len(prot_ids) == 1
    assert prot_ids[0].getSearchEngine() == "PeptideDataBaseSearchFI"


def test_empty_spectra():
    """Empty spectra should produce zero PSMs."""
    fasta_db = _build_fasta_db()[:1]

    algo = PeptideSearchEngineFIAlgorithm()
    p = algo.getParameters()
    p.setValue("decoys", "false")
    algo.setParameters(p)

    prot_ids = []
    pep_ids = PeptideIdentificationList()
    ec = algo.search(MSExperiment(), fasta_db, prot_ids, pep_ids)

    ExitCodes = PeptideSearchEngineFIAlgorithm.PeptideSearchEngineFIAlgorithm_ExitCodes
    assert ec == ExitCodes.EXECUTION_OK
    assert pep_ids.size() == 0
