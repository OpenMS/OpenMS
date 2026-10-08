// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include "MzCalibration.h"

using namespace std;

namespace OpenMS
{
  MzCalibration::MzCalibration() : mz_raw_ {}, mz_ref_ {}, no_mzml_(false)
  {
  }

  // find original m/z Value, set meta value "mz_raw" and set meta value "mz_ref"
  void MzCalibration::compute(FeatureMap& features, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum)
  {
    if (exp.empty())
    {
      no_mzml_ = true;
      OPENMS_LOG_WARN << "Metric MzCalibration received an empty mzml file. Only reporting uncalibrated mz error.\n";
    }
    else
    {
      no_mzml_ = false;
      // check for Calibration
      auto is_not_elem = [](const std::shared_ptr<const OpenMS::DataProcessing>& dp) { return (dp->getProcessingActions().count(DataProcessing::CALIBRATION) == 0); };
      auto vdp = exp[0].getDataProcessing(); // get a copy to avoid calling .begin() and .end() on two different temporaries
      if (all_of(vdp.begin(), vdp.end(), is_not_elem))
      {
        no_mzml_ = true;
        OPENMS_LOG_WARN << "Metric MzCalibration received an mzml file which did not undergo InternalCalibration. Only reporting uncalibrated mz error.\n";
      }
    }

    // set meta values for the top match of all identifications of the features and of the unassigned ones
    QCBase::annotateIdentifications(features, [&](Feature*, std::vector<QCBase::AnnotatedIdentification>& identifications) {
      for (auto& id : identifications)
      {
        addMzMetaValues_(id, exp, map_to_spectrum);
      }
    });
  }

  void MzCalibration::addMzMetaValues_(QCBase::AnnotatedIdentification& id, const PeakMap& exp, const QCBase::SpectraMap& map_to_spectrum)
  {
    if (id.top == nullptr)
    {
      return;
    }
    PeptideHit& hit = *id.top;

    mz_ref_ = hit.getSequence().getMZ(hit.getCharge());

    if (no_mzml_)
    {
      hit.setMetaValue("uncalibrated_mz_error_ppm", Math::getPPM(id.mz, mz_ref_));
    }
    else
    {
      if (id.spectrum_reference.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum reference annotated at peptide identification!");
      }

      // get spectrum from mapping and meta value
      MSSpectrum spectrum = exp[map_to_spectrum.at(id.spectrum_reference)];

      // check if spectrum fulfills all requirements
      if (spectrum.getMSLevel() == 2)
      {
        // meta value has to be there, because InternalCalibration had to be called to get here
        if (!spectrum.getPrecursors()[0].metaValueExists("mz_raw"))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected meta value 'mz_raw' at MSSpectrum, but could not find it.");
        }
        mz_raw_ = (double)spectrum.getPrecursors()[0].getMetaValue("mz_raw");
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The matching spectrum of the mzML is not an MS2 Spectrum.");
      }

      // set meta values
      hit.setMetaValue("mz_raw", mz_raw_);
      hit.setMetaValue("mz_ref", mz_ref_);
      hit.setMetaValue("uncalibrated_mz_error_ppm", Math::getPPM(mz_raw_, mz_ref_));
      hit.setMetaValue("calibrated_mz_error_ppm", Math::getPPM(id.mz, mz_ref_));
    }
  }

  // required input files
  QCBase::Status MzCalibration::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
  }

  const std::string& MzCalibration::getName() const
  {
    static const std::string& name = "MzCalibration";
    return name;
  }
} // namespace OpenMS
