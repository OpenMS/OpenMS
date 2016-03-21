// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_PILISCROSSVALIDATION_H
#define OPENMS_ANALYSIS_ID_PILISCROSSVALIDATION_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <vector>

namespace OpenMS
{
  class PeakSpectrumCompareFunctor;
  class PILISModel;


  /** @brief Implementation of a cross validation training for the PILIS model

    This class serves as an implementation of a cross validation training for
    the PILIS model. It includes a range of parameters which can be set to
    perform a GridSearch additionally.

    @htmlinclude OpenMS_PILISCrossValidation.parameters
*/
  class OPENMS_DLLAPI PILISCrossValidation :
    public DefaultParamHandler
  {

public:

    /** @brief this struct represents a peptide spectrum pair

    */
    struct OPENMS_DLLAPI Peptide
    {
      Peptide();

      Peptide(const Peptide & rhs);

      virtual ~Peptide();

      Peptide & operator=(const Peptide & rhs);

      AASequence sequence;
      Int charge;
      RichPeakSpectrum spec;

      std::vector<PeptideHit> hits;

      bool operator<(const Peptide & peptide) const;

    };

    /** @brief This struct represents a cross validation option

    */
    struct OPENMS_DLLAPI Option
    {
      /// Type of the parameters
      enum Type
      {
        INT = 0,
        DOUBLE = 1,
        BOOL = 2,
        STRINGLIST = 3
      };

      /// Default constructor
      Option();

      /// copy constructor
      Option(const Option & rhs);

      /// detailed constructors
      Option(Type t, double min, double max, double stepsize);

      /// assignment operator
      Option & operator=(const Option & rhs);

      Type type;
      Int int_min;
      Int int_max;
      Int int_stepsize;
      double dbl_min;
      double dbl_max;
      double dbl_stepsize;
    };


    /** @name Constructors and destructors
    */
    //@{
    /// Default constructor
    PILISCrossValidation();

    /// copy constructor
    PILISCrossValidation(const PILISCrossValidation & rhs);

    /// destructor
    virtual ~PILISCrossValidation();

    /// assignment operator
    PILISCrossValidation & operator=(const PILISCrossValidation & rhs);
    //@}

    /** @name Accessors
    */
    //@{
    /// sets the options which should be used for the cross validation
    void setOptions(const Map<String, Option> & rhs)
    {
      cv_options_ = rhs;
    }

    /// sets a option to be used for the cross validation
    void setOption(const String & name, const Option & option)
    {
      cv_options_[name] = option;
    }

    /// performs a cross validation and write optimized param into PILIS_param
    void apply(Param & PILIS_param, const PILISModel & base_model, const std::vector<Peptide> & peptides);

    /// compares experimental and simulated spectra and returns a score
    double scoreHits(const std::vector<std::vector<std::vector<RichPeakSpectrum> > > & sim_spectra, const std::vector<std::vector<RichPeakSpectrum> > & exp_spectra);
    //@}

protected:

    double scoreSpectra_(const RichPeakSpectrum & spec1, const RichPeakSpectrum & spec2);

    void partition_(std::vector<std::vector<Peptide> > & parts, const std::vector<Peptide> & source);

    void generateParameters_(const Param & param, const Map<String, Option> & options, std::vector<Param> & parameters);

    Map<String, Option> cv_options_;

    void updateMembers_();

    PeakSpectrumCompareFunctor * pscf_;
  };
}

#endif
