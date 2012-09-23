// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <iostream>

using namespace std;

namespace OpenMS
{
  const std::string EnzymaticDigestion::NamesOfEnzymes[] = {"Trypsin"};

  EnzymaticDigestion::EnzymaticDigestion() :
    missed_cleavages_(0),
    enzyme_(TRYPSIN),
    use_log_model_(false),
    log_model_threshold_(0.25),
    model_data_()
  {
    // load the cleavage model from disk (might throw exceptions)
    TextFile tf;
    tf.load(File::find("./CHEMISTRY/MissedCleavage.model"), true);
    for (Size i = 0; i < tf.size(); ++i)
    {
      if (tf[i].trim().hasPrefix("#"))
        continue;                                // skip comments
      StringList components;
      tf[i].split(' ', components);
      if (components.size() != 4)
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("split(' ',") + tf[i] + ")", String("Got ") + components.size() + " columns, expected 4!");
      BindingSite bs(components[0].toInt(), components[1].trim());
      CleavageModel cl(components[2].toDouble(), components[3].toDouble());
      model_data_[bs] = cl;
    }
  }

  SignedSize EnzymaticDigestion::getMissedCleavages() const
  {
    return missed_cleavages_;
  }

  void EnzymaticDigestion::setMissedCleavages(SignedSize missed_cleavages)
  {
    missed_cleavages_ = missed_cleavages;
  }

  EnzymaticDigestion::Enzyme EnzymaticDigestion::getEnzyme() const
  {
    return enzyme_;
  }

  void EnzymaticDigestion::setEnzyme(Enzyme enzyme)
  {
    if (enzyme < SIZE_OF_ENZYMES)
      enzyme_ = enzyme;
    else
      enzyme_ = SIZE_OF_ENZYMES;
  }

  EnzymaticDigestion::Enzyme EnzymaticDigestion::getEnzymeByName(const String & name)
  {
    if (name == "Trypsin")
      return TRYPSIN;
    else
      return SIZE_OF_ENZYMES;
  }

  bool EnzymaticDigestion::isLogModelEnabled() const
  {
    return use_log_model_;
  }

  void EnzymaticDigestion::setLogModelEnabled(bool enabled)
  {
    use_log_model_ = enabled;
  }

  DoubleReal EnzymaticDigestion::getLogThreshold() const
  {
    return log_model_threshold_;
  }

  void EnzymaticDigestion::setLogThreshold(DoubleReal threshold)
  {
    log_model_threshold_ = threshold;
  }

  void EnzymaticDigestion::nextCleavageSite_(const AASequence & protein, AASequence::ConstIterator & iterator)
  {
    switch (enzyme_)
    {
    case TRYPSIN:
      if (use_log_model_)
      {
        SignedSize pos = std::distance(AASequence::ConstIterator(protein.begin()), iterator) - 4;   // start position in sequence
        while (iterator != protein.end())
        {
          if (*iterator != 'R' && *iterator != 'K')   // wait for R or K
          {
            ++iterator;
            ++pos;
            continue;
          }
          DoubleReal score_cleave = 0, score_missed = 0;
          for (SignedSize i = 0; i < 9; ++i)
          {
            if ((pos + i >= 0) && (pos + i < (SignedSize)protein.size()))
            {
              CleavageModel p = model_data_[BindingSite(i, protein[pos + i].getOneLetterCode())]; // might not exists, but then its 0 by default
              score_cleave += p.p_cleave;
              score_missed += p.p_miss;
            }
          }
          ++iterator;
          ++pos;
          if (score_missed - score_cleave > log_model_threshold_)
            return;
        }

      }
      else   // naive digestion
      {
        while (iterator != protein.end())
        {
          //R or K at the end and not P afterwards
          if ((*iterator == 'R' || *iterator == 'K') && ((iterator + 1) == protein.end() || *(iterator + 1) != 'P'))
          {
            ++iterator;
            return;
          }
          ++iterator;
        }
      }
      break;

    default:
      return;
    }
  }

  Size EnzymaticDigestion::peptideCount(const AASequence & protein)
  {
    SignedSize count = 1;

    AASequence::ConstIterator iterator = protein.begin();
    while (nextCleavageSite_(protein, iterator), iterator != protein.end())
    {
      ++count;
    }
    if (use_log_model_)
      missed_cleavages_ = 0;                 // log model has missed cleavages build-in

    //missed cleavages
    Size sum = count;
    for (SignedSize i = 1; i < count; ++i)
    {
      if (i > missed_cleavages_)
        break;
      sum += count - i;
    }
    return sum;
  }

  void EnzymaticDigestion::digest(const AASequence & protein, std::vector<AASequence> & output)
  {
    //initialization
    SignedSize count = 1;
    output.clear();

    if (use_log_model_)
      missed_cleavages_ = 0;                 // log model has missed cleavages build-in

    //missed cleavage iterators
    vector<AASequence::ConstIterator> mc_iterators;
    if (missed_cleavages_ != 0)
      mc_iterators.push_back(protein.begin());

    AASequence::ConstIterator begin = protein.begin();
    AASequence::ConstIterator end = protein.begin();
    while (nextCleavageSite_(protein, end), end != protein.end())
    {
      ++count;
      if (missed_cleavages_ != 0)
      {
        mc_iterators.push_back(end);
      }

      output.push_back(protein.getSubsequence(begin - protein.begin(), end - begin));
      begin = end;
    }
    output.push_back(protein.getSubsequence(begin - protein.begin(), end - begin));
    if (missed_cleavages_ != 0)
    {
      mc_iterators.push_back(end);
    }

    //missed cleavages
    if (mc_iterators.size() > 2)     //there is at least one cleavage site!
    {
      //resize to number of fragments
      Size sum = count;

      for (SignedSize i = 1; i < count; ++i)
      {
        if (i > missed_cleavages_)
          break;
        sum += count - i;
      }
      output.resize(sum);

      //generate fragments with missed cleavages
      Size pos = count;
      for (SignedSize i = 1; ((i <= missed_cleavages_) && (count > i)); ++i)
      {
        vector<AASequence::ConstIterator>::const_iterator b = mc_iterators.begin();
        vector<AASequence::ConstIterator>::const_iterator e = b + (i + 1);
        while (e != mc_iterators.end())
        {
          output[pos] = AASequence(protein.getSubsequence(*b - protein.begin(), *e  - *b));
          ++b;
          ++e;
          ++pos;
        }
      }
    }
  }

}   //namespace
