// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <vector>

namespace OpenMS
{

  struct OPENMS_DLLAPI PrecursorMassComparator
  {
    bool operator() (const MSSpectrum& a, const MSSpectrum& b)
    {
      return a.getPrecursors()[0].getMZ() < b.getPrecursors()[0].getMZ();
    }

  } PrecursorMZLess;

  class OPENMS_DLLAPI SpectralMatch
  {
  public:
    /// Default constructor
    SpectralMatch();

    /// Default destructor
    ~SpectralMatch();

    /// Copy constructor
    SpectralMatch(const SpectralMatch&);

    /// Assignment operator
    SpectralMatch& operator=(const SpectralMatch&);

    double getObservedPrecursorMass() const;
    void setObservedPrecursorMass(const double&);

    double getObservedPrecursorRT() const;
    void setObservedPrecursorRT(const double&);

    double getFoundPrecursorMass() const;
    void setFoundPrecursorMass(const double&);

    Int getFoundPrecursorCharge() const;
    void setFoundPrecursorCharge(const Int&);

    double getMatchingScore() const;
    void setMatchingScore(const double&);

    Size getObservedSpectrumIndex() const;
    void setObservedSpectrumIndex(const Size&);

    Size getMatchingSpectrumIndex() const;
    void setMatchingSpectrumIndex(const Size&);

    String getObservedSpectrumNativeID() const;
    void setObservedSpectrumNativeID(const String&);

    String getPrimaryIdentifier() const;
    void setPrimaryIdentifier(const String&);

    String getSecondaryIdentifier() const;
    void setSecondaryIdentifier(const String&);

    String getCommonName() const;
    void setCommonName(const String&);

    String getSumFormula() const;
    void setSumFormula(const String&);

    String getInchiString() const;
    void setInchiString(const String&);

    String getSMILESString() const;
    void setSMILESString(const String&);

    String getPrecursorAdduct() const;
    void setPrecursorAdduct(const String&);


  private:
    double observed_precursor_mass_;
    double observed_precursor_rt_;
    double found_precursor_mass_;
    Int found_precursor_charge_;
    double matching_score_;
    Size observed_spectrum_idx_;
    Size matching_spectrum_idx_;
    String observed_spectrum_native_id_;

    // further meta information
    String primary_id_;
    String secondary_id_;
    String common_name_;
    String sum_formula_;
    String inchi_string_;
    String smiles_string_;
    String precursor_adduct_;

  };

  struct OPENMS_DLLAPI SpectralMatchScoreComparator
  {
    bool operator() (const SpectralMatch& a, const SpectralMatch& b)
    {
      return a.getMatchingScore() > b.getMatchingScore();
    }

  } SpectralMatchScoreGreater;

  class OPENMS_DLLAPI MetaboliteSpectralMatching :
  public DefaultParamHandler,
  public ProgressLogger
  {
  public:
    /// Default constructor
    MetaboliteSpectralMatching();

    /// Default destructor
    ~MetaboliteSpectralMatching() override;

    /// hyperscore computation
    static double computeHyperScore(
      double fragment_mass_error,
      bool fragment_mass_tolerance_unit_ppm,
      const MSSpectrum& exp_spectrum,
      const MSSpectrum& db_spectrum,
      double mz_lower_bound = 0.0);

    /// hyperscore computation (with output of peak annotations)
    static double computeHyperScore(
      double fragment_mass_error,
      bool fragment_mass_tolerance_unit_ppm,
      const MSSpectrum& exp_spectrum,
      const MSSpectrum& db_spectrum,
      std::vector<PeptideHit::PeakAnnotation>& annotations,
      double mz_lower_bound = 0.0);

    /// main method of MetaboliteSpectralMatching
    void run(PeakMap &, PeakMap &, MzTab &, String &);

  protected:
    void updateMembers_() override;

    // we have to use a pointer for "annotations" because mutable
    // references can't have temporary default values:
    static double computeHyperScore_(
      double fragment_mass_error,
      bool fragment_mass_tolerance_unit_ppm,
      const MSSpectrum& exp_spectrum,
      const MSSpectrum& db_spectrum,
      std::vector<PeptideHit::PeakAnnotation>* annotations = 0,
      double mz_lower_bound = 0.0);

  private:
    /// private member functions
    void exportMzTab_(const std::vector<SpectralMatch>&, MzTab&);

    double precursor_mz_error_;
    double fragment_mz_error_;
    String mz_error_unit_;
    String ion_mode_;

    String report_mode_;

    bool merge_spectra_;
  };

}

