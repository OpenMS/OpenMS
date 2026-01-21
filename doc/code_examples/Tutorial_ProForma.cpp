// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_ProForma]

// This tutorial demonstrates how to use the ProForma v2 parser and writer
// for representing peptides with modifications.
//
// ProForma is a standardized notation for proteoforms developed by the
// HUPO Proteomics Standards Initiative (PSI).

#include <OpenMS/CHEMISTRY/ProFormaParser.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // ===========================================
  // Basic Parsing
  // ===========================================
  cout << "=== Basic Parsing ===" << endl;

  // Parse a simple peptide
  Peptidoform pf = ProFormaParser::parse("PEPTIDE");
  cout << "Simple: " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Parse with UNIMOD modification (Oxidation on M)
  pf = ProFormaParser::parse("PEM[UNIMOD:35]TIDE");
  cout << "UNIMOD: " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Parse with named modification
  pf = ProFormaParser::parse("PES[Phospho]TIDE");
  cout << "Named:  " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Parse with mass delta
  pf = ProFormaParser::parse("PEM[+15.9949]TIDE");
  cout << "Delta:  " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;
  cout << endl;

  // ===========================================
  // Terminal Modifications
  // ===========================================
  cout << "=== Terminal Modifications ===" << endl;

  // N-terminal acetylation
  pf = ProFormaParser::parse("[Acetyl]-PEPTIDE");
  cout << "N-term: " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // C-terminal amidation
  pf = ProFormaParser::parse("PEPTIDE-[Amidated]");
  cout << "C-term: " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Both terminals modified
  pf = ProFormaParser::parse("[Acetyl]-PEPTIDE-[Amidated]");
  cout << "Both:   " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;
  cout << endl;

  // ===========================================
  // Write Modes: LOSSLESS vs CANONICAL
  // ===========================================
  cout << "=== Write Modes ===" << endl;

  // Parse with limited precision mass delta
  pf = ProFormaParser::parse("EM[+15.99]K");

  // Lossless preserves original formatting
  cout << "Lossless:  " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Canonical normalizes to 4 decimal places
  cout << "Canonical: " << ProFormaParser::toString(pf, ProFormaWriteMode::CANONICAL) << endl;
  cout << endl;

  // ===========================================
  // Advanced Features
  // ===========================================
  cout << "=== Advanced Features ===" << endl;

  // Global isotope labeling (SILAC)
  pf = ProFormaParser::parse("<13C>PEPTIDE");
  cout << "Isotope label: " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Global modification on specific residues (TMT)
  pf = ProFormaParser::parse("<[TMT6plex]@K,N-term>PEPTIDEK");
  cout << "TMT labeling:  " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Unlocalised modification (position unknown)
  pf = ProFormaParser::parse("[Phospho]?STYPEPTIDE");
  cout << "Unlocalised:   " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;

  // Labile modification (lost during fragmentation)
  pf = ProFormaParser::parse("{Glycan:Hex}PEPTIDE");
  cout << "Labile:        " << ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS) << endl;
  cout << endl;

  // ===========================================
  // Multi-chain Peptides and Charge States
  // ===========================================
  cout << "=== Multi-chain and Charge States ===" << endl;

  // Parse with charge state (returns PeptidoformIon)
  PeptidoformIon pfi = ProFormaParser::parseIon("PEPTIDE/2");
  cout << "Charged ion - chains: " << pfi.chains.size();
  if (pfi.charge.has_value() && std::holds_alternative<int>(pfi.charge.value()))
  {
    cout << ", charge: " << std::get<int>(pfi.charge.value());
  }
  cout << endl;

  // Cross-linked peptides
  pfi = ProFormaParser::parseIon("PEPK[+138.068#XL1]IDE//ANOK[#XL1]THER/3");
  cout << "Cross-linked - chains: " << pfi.chains.size() << endl;
  cout << endl;

  // ===========================================
  // Mass Calculation
  // ===========================================
  cout << "=== Mass Calculation ===" << endl;

  pf = ProFormaParser::parse("PEM[UNIMOD:35]TIDE");
  ProFormaParser::resolveModifications(pf);

  if (ProFormaParser::canCalculateMass(pf))
  {
    double mass = ProFormaParser::getMonoWeight(pf);
    cout << "Monoisotopic mass: " << mass << " Da" << endl;

    double mz = ProFormaParser::getMZ(pf, 2);
    cout << "m/z at charge 2+: " << mz << endl;
  }
  cout << endl;

  // ===========================================
  // AASequence Conversion
  // ===========================================
  cout << "=== AASequence Conversion ===" << endl;

  pf = ProFormaParser::parse("PEM[UNIMOD:35]TIDE");

  if (ProFormaParser::isRepresentableAsAASequence(pf))
  {
    ProFormaParser::resolveModifications(pf);
    AASequence aas = ProFormaParser::toAASequence(pf, AASequenceConversionPolicy::FAIL_ON_LOSS);
    cout << "AASequence: " << aas.toString() << endl;
    cout << "Mono mass:  " << aas.getMonoWeight() << " Da" << endl;
  }

  // Convert AASequence back to ProForma
  AASequence aas = AASequence::fromString("PEM(Oxidation)TIDE");
  pf = ProFormaParser::fromAASequence(aas);
  cout << "ProForma:   " << ProFormaParser::toString(pf, ProFormaWriteMode::CANONICAL) << endl;
  cout << endl;

  // ===========================================
  // Theoretical Spectrum Generation
  // ===========================================
  cout << "=== Spectrum Generation ===" << endl;

  pf = ProFormaParser::parse("PEPTIDE");
  ProFormaParser::resolveModifications(pf);

  if (ProFormaParser::canGenerateSpectrum(pf))
  {
    // Generate b/y ions at charge 1-2
    MSSpectrum spec = ProFormaParser::generateSpectrum(pf, 1, 2, "by", false, true);
    cout << "Generated spectrum with " << spec.size() << " peaks" << endl;
  }

  // Cross-linked peptide spectrum
  pfi = ProFormaParser::parseIon("PEPK[+138.068#XL1]IDE//ANOK[#XL1]THER/3");
  ProFormaParser::resolveModifications(pfi.chains[0]);
  ProFormaParser::resolveModifications(pfi.chains[1]);

  if (ProFormaParser::canGenerateSpectrum(pfi))
  {
    MSSpectrum spec = ProFormaParser::generateSpectrum(pfi, 1, 2, "aby", false, true);
    cout << "Cross-linked spectrum with " << spec.size() << " peaks" << endl;
  }
  cout << endl;

  // ===========================================
  // JSON Serialization
  // ===========================================
  cout << "=== JSON Serialization ===" << endl;

  pf = ProFormaParser::parse("[Acetyl]-PEM[UNIMOD:35]TIDE");
  String json_str = toJSON(pf);
  cout << "JSON length: " << json_str.size() << " bytes" << endl;

  // Roundtrip: deserialize and compare
  Peptidoform restored = peptidoformFromJSON(json_str);
  String original = ProFormaParser::toString(pf, ProFormaWriteMode::LOSSLESS);
  String roundtrip = ProFormaParser::toString(restored, ProFormaWriteMode::LOSSLESS);
  cout << "Roundtrip OK: " << (original == roundtrip ? "yes" : "no") << endl;
  cout << endl;

  cout << "Tutorial completed successfully!" << endl;

  return 0;
}

//! [doxygen_snippet_ProForma]
