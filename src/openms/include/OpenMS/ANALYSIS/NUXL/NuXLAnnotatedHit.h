// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideHit.h> // for PeakAnnotation
#include <OpenMS/DATASTRUCTURES/StringView.h>

namespace OpenMS
{
/// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
/// floats need to be initialized to zero as default
class NuXLAnnotatedHit
{
  public:
  /*
     Slim indices/views to lookup the actual sequence
   */
  StringView sequence;
  SignedSize peptide_mod_index = 0; // enumeration index of the non-NA peptide modification
  Size NA_mod_index = 0; // index of the NA modification
  Size NA_adduct_amb_index = 0; // store index the entry in the set of ambiguous precursor adducts (e.g, C-NH3 vs. U-H2O)

  int isotope_error = 0; // wheter the hit has been matched with isotopic misassignment

  static constexpr const char UNKNOWN_NUCLEOTIDE = '?';
  char cross_linked_nucleotide = UNKNOWN_NUCLEOTIDE;

  /**
       The main score (score) is a linear combination of (weighted) subscores

       For the fast score (ignoring all shifted peaks) we calculate: 
         score = 1.0 * total_loss_score 
               + 1.0 * total_MIC 
               + 0.333 * mass_error_p;
 
       For the all-ion score we calculate:
         peptides:
	      score = -6.486416409280039 
               + 4.059968526608637   * ah.total_MIC         
               + 0.5842539236790404  * ah.modds
               + 0.21721652155697285 * ah.total_loss_score
               + 1.9988345415208777  * ah.mass_error_p;
         XLs:
	       score = -6.648631037190969
               + 0.4688059636415974  * ah.Morph
               + 4.0386886051238     * ah.MIC         
               + 0.5446999629799386  * ah.modds
               + 0.25318342707227187 * ah.total_loss_score
               + 0.12472562244230834 * ah.partial_loss_score
               + 1.2107674392113372  * ah.mass_error_p
               + 2.3319284783288805  * ah.pl_MIC;
  */
  float score = 0;

  /**
       Normalized precursor mass error score.
       Mass error is assumed normally distributed with:
         - mean = 0
         - sd = sqrt(precursor_mass_tolerance) => variance = precursor tolerance
  */
  float mass_error_p = 0;

  //
  // Scores exclusively calculated from peaks without nucleotide shifts:
  //
  
  /**
      The total loss score is the X!Tandem HyperScore calculated from b-,y-ions 
      without any nucleotide shift.
  */
  float total_loss_score = 0;

  /**
      The matched ion current in immonium (immonium_score) and precursor ions (precursor_score) 
      without any nucleotide shift.

      see DOI: 10.1021/pr3007045 A Systematic Investigation into the Nature of Tryptic HCD Spectra
      imY = EmpiricalFormula("C8H10NO").getMonoWeight(); // 85%
      imW = EmpiricalFormula("C10H11N2").getMonoWeight(); // 84%
      imF = EmpiricalFormula("C8H10N").getMonoWeight(); // 84%
      imL = EmpiricalFormula("C5H12N").getMonoWeight(); // I/L 76%
      imH = EmpiricalFormula("C5H8N3").getMonoWeight(); // 70%
      imC = EmpiricalFormula("C2H6NS").getMonoWeight(); // CaC 61%
      imK1 = EmpiricalFormula("C5H13N2").getMonoWeight(); // 2%
      imP = EmpiricalFormula("C4H8N").getMonoWeight(); //?
      imQ = 101.0715; // 52%
      imE = 102.0555; // 37%
      imM = 104.0534; // 3%
  */
  float immonium_score = 0;
  float precursor_score = 0;

  /**
      The matched ion current (MIC), average fragment error (err), and morpheus score (Morph) are calculated 
      for b-,y-,a-ions without nucleotide shift. Morph is just the number of matched peaks + the fraction of MIC
  */
  float MIC = 0;
  float err = 0;
  float Morph = 0;

  /**
     The match odds (-log10) of observing this number of b-,a-, and y-ions assuming a uniform distribution of noise peaks.      
  */
  float modds = 0;

  //
  // Scores exclusively calculated from nucleotide shifted peaks:
  //
 
  /**
      The partial loss score is the X!Tandem HyperScore calculated from b-,a-, and y-ions 
      with nucleotide shifts. Matches from b- and a-ions are combined, i.e. a matching a_n-ion is counted as b_n-ion.
      For a precursor with charge N, all fragment ion charges up to N-1 are considered.

      Calculation of HyperScore:
      yFact = logfactorial_(y_ion_count);
      bFact = logfactorial_(b_ion_count);
      hyperScore = log1p(dot_product) + yFact + bFact;
  */
  float partial_loss_score = 0;

  /**
      The matched ion current (pl_MIC) of ladder ions, average fragment error (pl_err), and morpheus score (pl_Morph) are calculated 
      from b-,y-,a-ions with nucleotide shift.
      Morph: number of matched peaks + the fraction of MIC
  */
  float pl_MIC = 0;
  float pl_err = 0.0;
  float pl_Morph = 0;

  /*
     The match odds (-log10) of observing this number of b-,a-, and y-ions with nucleotide shifts assuming a uniform distribution of noise peaks.      
  */
  float pl_modds = 0;

  /*
     The MIC of precursor with all nucleotide shifts.
     Three variants: No additional loss, loss of water, and loss ammonia.
     Charge states considered: 1..N (precursor charge)
  */
  float pl_pc_MIC = 0;

  /**
      The matched ion current calculated from immonium ions with nucleotide shifts.
      Only singly charged immonium ions are considered.

      imY = EmpiricalFormula("C8H10NO").getMonoWeight();
      imW = EmpiricalFormula("C10H11N2").getMonoWeight();
      imF = EmpiricalFormula("C8H10N").getMonoWeight();
      imH = EmpiricalFormula("C5H8N3").getMonoWeight();
      imC = EmpiricalFormula("C2H6NS").getMonoWeight();
      imP = EmpiricalFormula("C4H8N").getMonoWeight();
      imL = EmpiricalFormula("C5H12N").getMonoWeight();
      imK1 = EmpiricalFormula("C5H13N2").getMonoWeight();
      imK2 = EmpiricalFormula("C5H10N1").getMonoWeight();
      imK3 = EmpiricalFormula("C6H13N2O").getMonoWeight();
      imQ = 101.0715;
      imE = 102.0555;
      imM = 104.0534;
   */
  float pl_im_MIC = 0;

  //
  // Scores calculated from peaks with AND without nucleotide shifts:
  //
  
  /**
       The complete TIC fraction of explained peaks (total_MIC) (excludes marker ions)
       For peptides: total_MIC = MIC + im_MIC + pc_MIC (b-,a-,y-ions, immonium ions, precursor ions)
       For XLs:      total_MIC = MIC + im_MIC + pc_MIC + pl_MIC + pl_pc_MIC + pl_im_MIC + marker_ions_sub_score 
  */
  float total_MIC = 0;

  /**
       The matched ion current in marker ions (marker_ions_score) is not considered in scoring.
  */  
  float marker_ions_score = 0;

  /**
       Coverage of peptide by prefix or suffix ions (fraction)
       For example: PEPTIDER
                    01000100 (two of eight ions observed -> 2/8)       
       Shifted and non-shifted are combined to determine coverage.
  */
  float ladder_score = 0;
  /**
       Longest sequence covered in peptide by prefix or suffix ions (fraction).
       Coverage of peptide by prefix or suffix ions (fraction)
       For example: PEPTIDER
                    01110001 (three ions form the longest sequence -> 3/8)
       Shifted and non-shifted are combined to determine coverage.
  */
  float sequence_score = 0;

  float best_localization_score = 0;
  String localization_scores = 0;
  String best_localization;
  int best_localization_position = -1; // UNKNOWN
  std::vector<PeptideHit::PeakAnnotation> fragment_annotations;

  size_t tag_unshifted = 0;
  size_t tag_shifted = 0;
  size_t tag_XLed = 0;  // tag that contains the transition from unshifted to shifted

  double explained_peak_fraction = 0;
  double matched_theo_fraction = 0;
  double wTop50 = 0;

  size_t n_theoretical_peaks = 0;

  static bool hasBetterScore(const NuXLAnnotatedHit& a, const NuXLAnnotatedHit& b)
  {
    return a.score > b.score;
  }
};
}

