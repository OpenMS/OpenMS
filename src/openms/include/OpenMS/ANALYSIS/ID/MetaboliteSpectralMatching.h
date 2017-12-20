// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_METABOLITESPECTRALMATCHING_H
#define OPENMS_ANALYSIS_ID_METABOLITESPECTRALMATCHING_H

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
#include <algorithm>

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
    double computeHyperScore(MSSpectrum, MSSpectrum, const double&, const double&);

    /// main method of MetaboliteSpectralMatching
    void run(PeakMap &, PeakMap &, MzTab &);

  protected:
    void updateMembers_() override;

  private:
    /// private member functions
    void exportMzTab_(const std::vector<SpectralMatch>&, MzTab&);

    double precursor_mz_error_;
    double fragment_mz_error_;
    String mz_error_unit_;
    String ion_mode_;

    String report_mode_;
  };

}

#endif // OPENMS_ANALYSIS_ID_METABOLITESPECTRALMATCHING_H
