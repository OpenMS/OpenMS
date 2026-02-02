// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_ProForma]

// This tutorial demonstrates how to use the ProForma v2 parser and writer
// for representing peptides with modifications.
//
// ProForma is a standardized notation for proteoforms developed by the
// HUPO Proteomics Standards Initiative (PSI).

#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

// Type aliases for ProForma nested types
using Peptidoform = ProForma::Peptidoform;
using PeptidoformIon = ProForma::PeptidoformIon;
using WriteMode = ProForma::WriteMode;
using ConversionPolicy = ProForma::ConversionPolicy;

int main()
{
  // ===========================================
  // Basic Parsing
  // ===========================================
  cout << "=== Basic Parsing ===" << endl;

  // Parse a simple peptide
  Peptidoform pf = ProForma::parse("PEPTIDE");
  cout << "Simple: " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Parse with UNIMOD modification (Oxidation on M)
  pf = ProForma::parse("PEM[UNIMOD:35]TIDE");
  cout << "UNIMOD: " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Parse with named modification
  pf = ProForma::parse("PES[Phospho]TIDE");
  cout << "Named:  " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Parse with mass delta
  pf = ProForma::parse("PEM[+15.9949]TIDE");
  cout << "Delta:  " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;
  cout << endl;

  // ===========================================
  // Terminal Modifications
  // ===========================================
  cout << "=== Terminal Modifications ===" << endl;

  // N-terminal acetylation
  pf = ProForma::parse("[Acetyl]-PEPTIDE");
  cout << "N-term: " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // C-terminal amidation
  pf = ProForma::parse("PEPTIDE-[Amidated]");
  cout << "C-term: " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Both terminals modified
  pf = ProForma::parse("[Acetyl]-PEPTIDE-[Amidated]");
  cout << "Both:   " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;
  cout << endl;

  // ===========================================
  // Write Modes: LOSSLESS vs CANONICAL
  // ===========================================
  cout << "=== Write Modes ===" << endl;

  // Parse with limited precision mass delta
  pf = ProForma::parse("EM[+15.99]K");

  // Lossless preserves original formatting
  cout << "Lossless:  " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Canonical normalizes to 4 decimal places
  cout << "Canonical: " << ProForma::toString(pf, WriteMode::CANONICAL) << endl;
  cout << endl;

  // ===========================================
  // Advanced Features
  // ===========================================
  cout << "=== Advanced Features ===" << endl;

  // Global isotope labeling (SILAC)
  pf = ProForma::parse("<13C>PEPTIDE");
  cout << "Isotope label: " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Global modification on specific residues (TMT)
  pf = ProForma::parse("<[TMT6plex]@K,N-term>PEPTIDEK");
  cout << "TMT labeling:  " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Unlocalised modification (position unknown)
  pf = ProForma::parse("[Phospho]?STYPEPTIDE");
  cout << "Unlocalised:   " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;

  // Labile modification (lost during fragmentation)
  pf = ProForma::parse("{Glycan:Hex}PEPTIDE");
  cout << "Labile:        " << ProForma::toString(pf, WriteMode::LOSSLESS) << endl;
  cout << endl;

  // ===========================================
  // Multi-chain Peptides and Charge States
  // ===========================================
  cout << "=== Multi-chain and Charge States ===" << endl;

  // Parse with charge state (returns PeptidoformIon)
  PeptidoformIon pfi = ProForma::parseIon("PEPTIDE/2");
  cout << "Charged ion - chains: " << pfi.chains.size();
  if (pfi.charge.has_value() && std::holds_alternative<int>(pfi.charge.value()))
  {
    cout << ", charge: " << std::get<int>(pfi.charge.value());
  }
  cout << endl;

  // Cross-linked peptides
  pfi = ProForma::parseIon("PEPK[+138.068#XL1]IDE//ANOK[#XL1]THER/3");
  cout << "Cross-linked - chains: " << pfi.chains.size() << endl;
  cout << endl;

  // ===========================================
  // Mass Calculation
  // ===========================================
  cout << "=== Mass Calculation ===" << endl;

  pf = ProForma::parse("PEM[UNIMOD:35]TIDE");
  ProForma::resolveModifications(pf);

  if (ProForma::canCalculateMass(pf))
  {
    double mass = ProForma::getMonoWeight(pf);
    cout << "Monoisotopic mass: " << mass << " Da" << endl;

    double mz = ProForma::getMZ(pf, 2);
    cout << "m/z at charge 2+: " << mz << endl;
  }
  cout << endl;

  // ===========================================
  // AASequence Conversion
  // ===========================================
  cout << "=== AASequence Conversion ===" << endl;

  pf = ProForma::parse("PEM[UNIMOD:35]TIDE");

  if (ProForma::isRepresentableAsAASequence(pf))
  {
    ProForma::resolveModifications(pf);
    AASequence aas = ProForma::toAASequence(pf, ConversionPolicy::FAIL_ON_LOSS);
    cout << "AASequence: " << aas.toString() << endl;
    cout << "Mono mass:  " << aas.getMonoWeight() << " Da" << endl;
  }

  // Convert AASequence back to ProForma
  AASequence aas = AASequence::fromString("PEM(Oxidation)TIDE");
  pf = ProForma::fromAASequence(aas);
  cout << "ProForma:   " << ProForma::toString(pf, WriteMode::CANONICAL) << endl;
  cout << endl;

  // ===========================================
  // Theoretical Spectrum Generation
  // ===========================================
  cout << "=== Spectrum Generation ===" << endl;

  pf = ProForma::parse("PEPTIDE");
  ProForma::resolveModifications(pf);

  if (ProForma::canGenerateSpectrum(pf))
  {
    // Generate b/y ions at charge 1-2
    MSSpectrum spec = ProForma::generateSpectrum(pf, 1, 2, "by", false, true);
    cout << "Generated spectrum with " << spec.size() << " peaks" << endl;
  }

  // Cross-linked peptide spectrum
  pfi = ProForma::parseIon("PEPK[+138.068#XL1]IDE//ANOK[#XL1]THER/3");
  ProForma::resolveModifications(pfi.chains[0]);
  ProForma::resolveModifications(pfi.chains[1]);

  if (ProForma::canGenerateSpectrum(pfi))
  {
    MSSpectrum spec = ProForma::generateSpectrum(pfi, 1, 2, "aby", false, true);
    cout << "Cross-linked spectrum with " << spec.size() << " peaks" << endl;
  }
  cout << endl;

  // ===========================================
  // JSON Serialization
  // ===========================================
  cout << "=== JSON Serialization ===" << endl;

  pf = ProForma::parse("[Acetyl]-PEM[UNIMOD:35]TIDE");
  String json_str = ProForma::toJSON(pf);
  cout << "JSON length: " << json_str.size() << " bytes" << endl;

  // Roundtrip: deserialize and compare
  Peptidoform restored = ProForma::peptidoformFromJSON(json_str);
  String original = ProForma::toString(pf, WriteMode::LOSSLESS);
  String roundtrip = ProForma::toString(restored, WriteMode::LOSSLESS);
  cout << "Roundtrip OK: " << (original == roundtrip ? "yes" : "no") << endl;
  cout << endl;

  cout << "Tutorial completed successfully!" << endl;

  return 0;
}

//! [doxygen_snippet_ProForma]
