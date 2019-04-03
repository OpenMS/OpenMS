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
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#include <OpenMS/QC/FragmentMassError.h>
#include <include/OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

namespace OpenMS
{
  void FragmentMassError::compute(const MSExperiment& exp, FeatureMap& fmap)
  {
    //sequenz //const pep_id?
    auto lam = [](PeptideIdentification& pep_id)
    {
      if (pep_id.getHits().empty())
      {
        //Warn
        return;
      }

      //---------------------------------------------------------------------
      // CREATE THEORETICAL SPECTRUM
      //---------------------------------------------------------------------

      AASequence seq = pep_id.getHits()[0].getSequence();
      Int charge = pep_id.getHits()[0].getCharge();

      //initialize a TheoreticalSpectrumGenerator
      TheoreticalSpectrumGenerator theo_gen;

      //get current parameters (default)
      Param theo_gen_settings = theo_gen.getParameters();

      //default: b- and y-ions?
      theo_gen_settings.setValue("add_a_ions", "true");
      //theo_settings.setValue("add_b_ions", "true");
      theo_gen_settings.setValue("add_c_ions", "true");
      theo_gen_settings.setValue("add_x_ions", "true");
      //theo_settings.setValue("add_y_ions", "true");
      theo_gen_settings.setValue("add_z_ions", "true");

      //store ion types for each peak
      //theo_settings.setValue("add_metainfo", "true");

      //set changed parameters
      theo_gen.setParameters(theo_gen_settings);

      PeakSpectrum theo_spectrum;

      //generate a-, b- and y-ion spectrum of peptide seq with charge
      theo_gen.getSpectrum(theo_spectrum, seq, charge, charge);

      //-----------------------------------------------------------------------
      // GET EXPERIMENTAL SPECTRUM MATCHING TO PEPTIDEIDENTIFICTION
      //-----------------------------------------------------------------------

      //to be continued

      //-----------------------------------------------------------------------
      // COMPARE THEORETICAL AND EXPERIMENTAL SPECTRUM
      //-----------------------------------------------------------------------

      //to be continued

      //-----------------------------------------------------------------------
      // WRITE PPM ERROR IN PEPTIDEHIT
      //-----------------------------------------------------------------------

    };
  }

  const std::vector<double>& FragmentMassError::getResults() const
  {
    return result_;
  }

  QCBase::Status FragmentMassError::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }
} //namespace OpenMS