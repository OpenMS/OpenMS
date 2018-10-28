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
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost $
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/KERNEL/Peak1D.h>

// Override IsoSpec's use of mmap whenever it is available
#define ISOSPEC_GOT_SYSTEM_MMAN false
#define ISOSPEC_GOT_MMAN false

// IsoSpec doesn't (yet) support being compiled with -fvisibility=hidden
#pragma GCC visibility push(default)
#include <OpenMS/../../thirdparty/IsoSpec/IsoSpec/isoSpec++.h>
#pragma GCC visibility pop


namespace OpenMS
{

  /**
    * @brief Interface to the IsoSpec algorithm.
    * 
    * Provides an interface to the IsoSpec algorithm.
    *
    * @code
    * Łącki MK, Startek M, Valkenborg D, Gambin A.
    * IsoSpec: Hyperfast Fine Structure Calculator.
    * Anal Chem. 2017 Mar 21;89(6):3272-3277. doi: 10.1021/acs.analchem.6b01459.
    * @endcode
    *
    **/
  class OPENMS_DLLAPI IsoSpecWrapper
  {
public:

    /**
      * @brief Run the algorithm
      *
      **/
    virtual std::vector<Peak1D> run() = 0;


    virtual bool nextConf() = 0;
    virtual Peak1D getConf() = 0;
    virtual double getMass() = 0;
    virtual double getIntensity() = 0;
    virtual double getLogIntensity() = 0;

    virtual inline ~IsoSpecWrapper() {};
  };

  class OPENMS_DLLAPI IsoSpecThresholdWrapper : public IsoSpecWrapper
  {
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param threshold Intensity threshold: will only compute peaks above this threshold
      * @param absolute Whether the threshold is absolute or relative (relative to the most intense peak)
      *
      **/
  IsoSpecThresholdWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double threshold,
             bool absolute);

    /**
      * @brief Run the algorithm on a sum formula
      *
      **/
  IsoSpecThresholdWrapper(const std::string& formula, double threshold, bool absolute);

  virtual std::vector<Peak1D> run() override final;

  virtual inline bool nextConf() override final { return ITG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(ITG.mass(), ITG.prob()); };
  virtual inline double getMass() override final { return ITG.mass(); };
  virtual inline double getIntensity() override final { return ITG.prob(); };
  virtual inline double getLogIntensity() override final { return ITG.lprob(); };


protected:
  IsoSpec::IsoThresholdGenerator ITG;
  };

  class OPENMS_DLLAPI IsoSpecTotalProbWrapper : public IsoSpecWrapper
  {
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param total_prob Total coverage of probability space desired
      *
      **/
  IsoSpecTotalProbWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double total_prob);

    /**
      * @brief Run the algorithm on a sum formula
      *
      **/
  IsoSpecTotalProbWrapper(const std::string& formula, double total_prob);

  virtual std::vector<Peak1D> run() override final;

  virtual inline bool nextConf() override final { return ILG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(ILG.mass(), ILG.prob()); };
  virtual inline double getMass() override final { return ILG.mass(); };
  virtual inline double getIntensity() override final { return ILG.prob(); };
  virtual inline double getLogIntensity() override final { return ILG.lprob(); };

protected:
  IsoSpec::IsoLayeredGenerator ILG;
  };



  class OPENMS_DLLAPI IsoSpecOrderedGeneratorWrapper : public IsoSpecWrapper
  {
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param total_prob Total coverage of probability space desired
      *
      **/
  IsoSpecOrderedGeneratorWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities);

    /**
      * @brief Run the algorithm on a sum formula
      *
      **/
  IsoSpecOrderedGeneratorWrapper(const std::string& formula);

  virtual std::vector<Peak1D> run() override final
  { throw std::logic_error("There is no stop condition in OrderedGenerator - therefore it only makes sense to use it as a generator"); } ;

  virtual inline bool nextConf() override final { return IOG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(IOG.mass(), IOG.prob()); };
  virtual inline double getMass() override final { return IOG.mass(); };
  virtual inline double getIntensity() override final { return IOG.prob(); };
  virtual inline double getLogIntensity() override final { return IOG.lprob(); };

protected:
  IsoSpec::IsoOrderedGenerator IOG;
  };

}

