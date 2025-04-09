// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>

using namespace std;
namespace OpenMS
{

  SpectraMerger::SpectraMerger() :
    DefaultParamHandler("SpectraMerger")
  {
    // common
    defaults_.setValue("mz_binning_width", 5.0, "minimum m/z distance for two data points (profile data) or peaks (centroided data) to be considered distinct. Closer data points or peaks will be merged.", {"advanced"});
    defaults_.setMinFloat("mz_binning_width", 0.0);

    defaults_.setValue("mz_binning_width_unit",
                       "ppm",
                       "Unit in which the distance between two data points or peaks is given.",
                       {"advanced"});
    defaults_.setValidStrings("mz_binning_width_unit", {"Da", "ppm"});

    defaults_.setValue("sort_blocks",
                       "RT_ascending",
                       "Sort blocks by <?> before merging them (useful for precursor order)",
                       {"advanced"});
    defaults_.setValidStrings("sort_blocks", {"RT_ascending", "RT_descending"});

    // Gaussian average
    defaults_.setValue("average_gaussian:spectrum_type", "automatic", "Spectrum type of the MS level to be averaged");
    defaults_.setValidStrings("average_gaussian:spectrum_type", {"profile", "centroid", "automatic"});
    defaults_
        .setValue("average_gaussian:ms_level",
                  1,
                  "If set to be 0, each MS level will be merged from 1 to max. Otherwise, average spectra of this level. All other spectra remain unchanged.");
    defaults_.setMinInt("average_gaussian:ms_level", 0);
    defaults_.setValue("average_gaussian:rt_FWHM", 5.0, "FWHM of Gauss curve in seconds to be averaged over.");
    defaults_.setMinFloat("average_gaussian:rt_FWHM", 0.0);
    defaults_.setMaxFloat("average_gaussian:rt_FWHM", 10e10);
    defaults_.setValue("average_gaussian:cutoff",
                       0.01,
                       "Intensity cutoff for Gaussian. The Gaussian RT profile decreases from 1 at its apex to 0 at infinity. Spectra for which the intensity of the Gaussian drops below the cutoff do not contribute to the average.",
                       {"advanced"});
    defaults_.setMinFloat("average_gaussian:cutoff", 0.0);
    defaults_.setMaxFloat("average_gaussian:cutoff", 1.0);
    defaults_.setValue("average_gaussian:precursor_mass_tol",
                       0.0,
                       "PPM mass tolerance for precursor mass. If set, MSn (n>2) spectra of precursor masses within the tolerance are averaged.");
    defaults_.setValue("average_gaussian:precursor_max_charge",
                       1,
                       "Possible maximum precursor ion charge. Effective only when average_gaussian:precursor_mass_tol option is active.");
    defaults_.setMinFloat("average_gaussian:precursor_mass_tol", 0.0);
    defaults_.setMinInt("average_gaussian:precursor_max_charge", 1);

    // top-hat average
    defaults_.setValue("average_tophat:spectrum_type", "automatic", "Spectrum type of the MS level to be averaged");
    defaults_.setValidStrings("average_tophat:spectrum_type", {"profile", "centroid", "automatic"});
    defaults_
        .setValue("average_tophat:ms_level",
                  1,
                  "If set to be 0, each MS level will be merged from 1 to max. Otherwise, average spectra of this level. All other spectra remain unchanged.");
    defaults_.setMinInt("average_tophat:ms_level", 0);
    defaults_.setValue("average_tophat:rt_range",
                       5.0,
                       "RT range to be averaged over, i.e. +/-(RT range)/2 from each spectrum.");
    defaults_.setMinFloat("average_tophat:rt_range", 0.0);
    defaults_.setMaxFloat("average_tophat:rt_range", 10e10);
    defaults_.setValue("average_tophat:rt_unit", "scans", "Unit for RT range.");
    defaults_.setValidStrings("average_tophat:rt_unit", {"scans", "seconds"});

    // block merging
    defaults_.setValue("block_method:ms_levels",
                       ListUtils::create<Int>("1"),
                       "Merge spectra of this level. All spectra with other MS levels remain untouched.");
    defaults_.setMinInt("block_method:ms_levels", 1);
    defaults_.setValue("block_method:rt_block_size", 5, "Maximum number of scans to be summed up.");
    defaults_.setMinInt("block_method:rt_block_size", 1);

    defaults_.setValue("block_method:rt_max_length",
                       0.0,
                       "Maximum RT size of the block in seconds (0.0 = no size restriction).");
    defaults_.setMinFloat("block_method:rt_max_length", 0.0);
    defaults_.setMaxFloat("block_method:rt_max_length", 10e10);

    // same precursor MS/MS merging
    defaults_.setValue("precursor_method:mz_tolerance",
                       10e-5,
                       "Max m/z distance of the precursor entries of two spectra to be merged in [Da].");
    defaults_.setMinFloat("precursor_method:mz_tolerance", 0);
    defaults_.setValue("precursor_method:mass_tolerance",
                       .0,
                       "Max mass distance of the precursor entries of two spectra to be merged in [Da]. Active when set to a positive value.");
    defaults_.setMinFloat("precursor_method:mass_tolerance", 0);
    defaults_.setValue("precursor_method:rt_tolerance",
                       5.0,
                       "Max RT distance of the precursor entries of two spectra to be merged in [s].");
    defaults_.setMinFloat("precursor_method:rt_tolerance", 0);

    defaultsToParam_();
  }

  SpectraMerger::SpectraMerger(const SpectraMerger & source) :
    DefaultParamHandler(source), ProgressLogger() //we probably want a new ProgressLogger when we copy
  {
  }

  SpectraMerger::~SpectraMerger() = default;

  SpectraMerger & SpectraMerger::operator=(const SpectraMerger & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

}
