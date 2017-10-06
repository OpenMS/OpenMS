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
// $Maintainer: Chris Bielow, Xiao Liang$
// $Authors: Marc Sturm, Chris Bielow, Xiao Liang $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EnzymaticDigestionLogModel.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

  EnzymaticDigestionLogModel::EnzymaticDigestionLogModel() :
    enzyme_(ProteaseDB::getInstance()->getEnzyme("Trypsin")),
    log_model_threshold_(0.25),
    model_data_()
  {
    // load the cleavage model from disk (might throw exceptions)
    TextFile tf;
    tf.load(File::find("./CHEMISTRY/MissedCleavage.model"), true);
    for (TextFile::ConstIterator it = tf.begin(); it != tf.end(); ++it)
    {
      String tmp = *it;
      if (tmp.trim().hasPrefix("#"))
      {
        continue;  // skip comments
      }
      StringList components;
      tmp.split(' ', components);
      if (components.size() != 4)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("split(' ',") + tmp + ")", String("Got ") + components.size() + " columns, expected 4!");
      }
      BindingSite_ bs(components[0].toInt(), components[1].trim());
      CleavageModel_ cl(components[2].toDouble(), components[3].toDouble());
      model_data_[bs] = cl;
    }
  }

  EnzymaticDigestionLogModel::EnzymaticDigestionLogModel(const EnzymaticDigestionLogModel& rhs) :
    enzyme_(rhs.enzyme_),
    log_model_threshold_(rhs.log_model_threshold_),
    model_data_(rhs.model_data_)
  {

  }

  /// Assignment operator
  EnzymaticDigestionLogModel& EnzymaticDigestionLogModel::operator=(const EnzymaticDigestionLogModel& rhs)
  {
    if (this != &rhs)
    {
      enzyme_ = rhs.enzyme_;
      log_model_threshold_ = rhs.log_model_threshold_;
      model_data_ = rhs.model_data_;
    }
    return *this;
  }

  void EnzymaticDigestionLogModel::setEnzyme(const String enzyme_name)
  {
    enzyme_ = ProteaseDB::getInstance()->getEnzyme(enzyme_name);
  }

  String EnzymaticDigestionLogModel::getEnzymeName() const
  {
    return enzyme_->getName();
  }

  double EnzymaticDigestionLogModel::getLogThreshold() const
  {
    return log_model_threshold_;
  }

  void EnzymaticDigestionLogModel::setLogThreshold(double threshold)
  {
    log_model_threshold_ = threshold;
  }

  bool EnzymaticDigestionLogModel::isCleavageSite_(
    const AASequence& protein, const AASequence::ConstIterator& iterator) const
  {
    if (enzyme_->getName() != "Trypsin") // no cleavage
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("EnzymaticDigestionLogModel: enzyme '") + enzyme_->getName() + " does not support logModel!");
    }
    if ((!enzyme_->getRegEx().hasSubstring(iterator->getOneLetterCode())) || *iterator == 'P') // wait for R or K
    {
      return false;
    }
    const SignedSize pos = distance(AASequence::ConstIterator(protein.begin()),
                                    iterator) - 4; // start position in sequence
    double score_cleave = 0, score_missed = 0;
    for (SignedSize i = 0; i < 9; ++i)
    {
      if ((pos + i >= 0) && (pos + i < (SignedSize) protein.size()))
      {
        BindingSite_ bs(i, protein[pos + i].getOneLetterCode());
        Map<BindingSite_, CleavageModel_>::const_iterator pos_it =
        model_data_.find(bs);
        if (pos_it != model_data_.end()) // no data for non-std. amino acids
        {
          score_cleave += pos_it->second.p_cleave;
          score_missed += pos_it->second.p_miss;
        }
      }
    }
    return score_missed - score_cleave > log_model_threshold_;
  }

  void EnzymaticDigestionLogModel::nextCleavageSite_(const AASequence& protein, AASequence::ConstIterator& iterator) const
  {
    for (; iterator != protein.end(); ++iterator)
    {
      if (isCleavageSite_(protein, iterator))
      {
        ++iterator;
        return;
      }
    }
  }

  Size EnzymaticDigestionLogModel::peptideCount(const AASequence& protein)
  {
    Size count = 0;
    for (AASequence::ConstIterator iterator = protein.begin();
         iterator != protein.end(); nextCleavageSite_(protein, iterator))
    {
        ++count;
    }
    return count;
  }

  void EnzymaticDigestionLogModel::digest(const AASequence& protein, vector<AASequence>& output) const
  {
    // initialization
    output.clear();
    AASequence::ConstIterator begin = protein.begin();
    AASequence::ConstIterator end = protein.begin();
    for (nextCleavageSite_(protein, end);
         begin != protein.end();
         begin = end, nextCleavageSite_(protein, end))
    {
      output.push_back(protein.getSubsequence(begin - protein.begin(), end - begin));
    }
  }
} //namespace

