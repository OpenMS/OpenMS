// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <iostream>
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
    struct Peptide
    {
      Peptide() :
        charge(0)
      {

      }

      Peptide(const Peptide & rhs) :
        sequence(rhs.sequence),
        charge(rhs.charge),
        spec(rhs.spec),
        hits(rhs.hits)
      {
      }

      virtual ~Peptide()
      {
      }

      Peptide & operator=(const Peptide & rhs)
      {
        if (&rhs != this)
        {
          sequence = rhs.sequence;
          charge = rhs.charge;
          spec = rhs.spec;
          hits = rhs.hits;
        }
        return *this;
      }

      AASequence sequence;
      Int charge;
      RichPeakSpectrum spec;

      std::vector<PeptideHit> hits;

      bool operator<(const Peptide & peptide) const
      {
        if (sequence < peptide.sequence)
        {
          return true;
        }
        else
        {
          if (sequence == peptide.sequence)
          {
            return charge < peptide.charge;
          }
        }
        return false;
      }

    };

    /** @brief This struct represents a cross validation option

    */
    struct Option
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
      Option() :
        type(INT),
        int_min(0),
        int_max(0),
        int_stepsize(0),
        dbl_min(0),
        dbl_max(0),
        dbl_stepsize(0)
      {
      }

      /// copy constructor
      Option(const Option & rhs) :
        type(rhs.type),
        int_min(rhs.int_min),
        int_max(rhs.int_max),
        int_stepsize(rhs.int_stepsize),
        dbl_min(rhs.dbl_min),
        dbl_max(rhs.dbl_max),
        dbl_stepsize(rhs.dbl_stepsize)
      {
      }

      /// detailed constructors
      Option(Type t, DoubleReal min, DoubleReal max, DoubleReal stepsize)
      {
        type = t;
        if (type == INT)
        {
          int_min = (Int)min;
          int_max = (Int)max;
          int_stepsize = (Int)stepsize;
        }
        else
        {
          if (type == DOUBLE)
          {
            dbl_min = min;
            dbl_max = max;
            dbl_stepsize = stepsize;
          }
          else
          {
            std::cerr << "Type: " << t << " is not known!" << std::endl;
          }
        }
      }

      /// assignment operator
      Option & operator=(const Option & rhs)
      {
        if (&rhs != this)
        {
          type = rhs.type;
          int_min = rhs.int_min;
          int_max = rhs.int_max;
          int_stepsize = rhs.int_stepsize;
          dbl_min = rhs.dbl_min;
          dbl_max = rhs.dbl_max;
          dbl_stepsize = rhs.dbl_stepsize;
        }
        return *this;
      }

      Type type;
      Int int_min;
      Int int_max;
      Int int_stepsize;
      DoubleReal dbl_min;
      DoubleReal dbl_max;
      DoubleReal dbl_stepsize;
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
    DoubleReal scoreHits(const std::vector<std::vector<std::vector<RichPeakSpectrum> > > & sim_spectra, const std::vector<std::vector<RichPeakSpectrum> > & exp_spectra);
    //@}

protected:

    DoubleReal scoreSpectra_(const RichPeakSpectrum & spec1, const RichPeakSpectrum & spec2);

    void partition_(std::vector<std::vector<Peptide> > & parts, const std::vector<Peptide> & source);

    void generateParameters_(const Param & param, const Map<String, Option> & options, std::vector<Param> & parameters);

    Map<String, Option> cv_options_;

    void updateMembers_();

    PeakSpectrumCompareFunctor * pscf_;
  };
}

#endif
