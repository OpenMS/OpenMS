// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#include <OpenMS/QC/MzCalibration.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/METADATA/DataProcessing.h>

using namespace std;
 
namespace OpenMS
{
  MzCalibration::MzCalibration() :
      mz_raw_{}, mz_ref_{}, no_mzml_(false), name_("MzCalibration")
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
      auto is_not_elem = [](const boost::shared_ptr<const OpenMS::DataProcessing> &dp)
      {
        return (dp->getProcessingActions().count(DataProcessing::CALIBRATION) == 0);
      };
      auto vdp = exp[0].getDataProcessing(); // get a copy to avoid calling .begin() and .end() on two different temporaries
      if (all_of(vdp.begin(), vdp.end(), is_not_elem))
      {
        no_mzml_ = true;
        OPENMS_LOG_WARN << "Metric MzCalibration received an mzml file which did not undergo InternalCalibration. Only reporting uncalibrated mz error.\n";
      }
    }

    // set meta values for the first hit of all PeptideIdentifications of all features
    for (Feature &feature : features)
    {
      if (feature.getPeptideIdentifications().empty())
      {
        continue;
      }

      for (PeptideIdentification &peptide_ID : feature.getPeptideIdentifications())
      {
        addMzMetaValues_(peptide_ID, exp, map_to_spectrum);
      }
    }
    // set meta values for the first hit of all unasssigned PeptideIdentifications
    for (PeptideIdentification &unassigned_ID : features.getUnassignedPeptideIdentifications())
    {
      addMzMetaValues_(unassigned_ID, exp, map_to_spectrum);
    }
  }

  void MzCalibration::addMzMetaValues_(PeptideIdentification &peptide_ID,
                                       const PeakMap &exp,
                                       const QCBase::SpectraMap &map_to_spectrum)
  {
    if (peptide_ID.getHits().empty())
    {
      return;
    }

    mz_ref_ = (peptide_ID.getHits()[0].getSequence()
                                      .getMonoWeight(OpenMS::Residue::Full, peptide_ID.getHits()[0].getCharge()))
              / peptide_ID.getHits()[0].getCharge();

    if (no_mzml_)
    {
      peptide_ID.getHits()[0].setMetaValue("uncalibrated_mz_error_ppm", Math::getPPM(peptide_ID.getMZ(), mz_ref_));
    }
    else
    {
      if (!peptide_ID.metaValueExists("spectrum_reference"))
      {
        throw Exception::InvalidParameter(__FILE__,
                                          __LINE__,
                                          OPENMS_PRETTY_FUNCTION,
                                          "No spectrum reference annotated at peptide identification!");
      }

      // get spectrum from mapping and meta value
      MSSpectrum spectrum = exp[map_to_spectrum.at(peptide_ID.getMetaValue("spectrum_reference").toString())];

      // check if spectrum fulfills all requirements
      if (spectrum.getMSLevel() == 2)
      {
        // meta value has to be there, because InternalCalibration had to be called to get here
        if (!spectrum.getPrecursors()[0].metaValueExists("mz_raw"))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected meta value 'mz_raw' at MSSpectrum, but could not find it.");
        }
        mz_raw_ = spectrum.getPrecursors()[0].getMetaValue("mz_raw");
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__,
                                         __LINE__,
                                         OPENMS_PRETTY_FUNCTION,
                                         "The matching spectrum of the mzML is not a MS2 Spectrum.");
      }

      // set meta values
      peptide_ID.getHits()[0].setMetaValue("mz_raw", mz_raw_);
      peptide_ID.getHits()[0].setMetaValue("mz_ref", mz_ref_);
      peptide_ID.getHits()[0].setMetaValue("uncalibrated_mz_error_ppm", Math::getPPM(mz_raw_, mz_ref_));
      peptide_ID.getHits()[0].setMetaValue("calibrated_mz_error_ppm", Math::getPPM(peptide_ID.getMZ(), mz_ref_));
    }
  }

  // required input files
  QCBase::Status MzCalibration::requires() const
  {
    return QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
  }

  const String &MzCalibration::getName() const
  {
    return name_;
  }
}
