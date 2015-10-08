// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexMassPatternList.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

  MultiplexMassPatternList::MultiplexMassPatternList(String labels, int missed_cleavages, bool knock_out, std::map<String,double> label_mass_shift) :
    labels_(labels), samples_labels_(), missed_cleavages_(missed_cleavages), knock_out_(knock_out), label_mass_shift_(label_mass_shift)
  {
    // split the labels_ string
    
    //std::vector<std::vector<String> > samples_labels;
    String temp_labels(labels_);
    std::vector<String> temp_samples;
    
    boost::replace_all(temp_labels, "[]", "no_label");
    boost::replace_all(temp_labels, "()", "no_label");
    boost::replace_all(temp_labels, "{}", "no_label");
    boost::split(temp_samples, temp_labels, boost::is_any_of("[](){}")); // any bracket allowed to separate samples
    
    for (unsigned i = 0; i < temp_samples.size(); ++i)
    {
      if (!temp_samples[i].empty())
      {
        if (temp_samples[i]=="no_label")
        {
          vector<String> temp_labels;
          temp_labels.push_back("no_label");
          samples_labels_.push_back(temp_labels);
        }
        else
        {
          vector<String> temp_labels;
          boost::split(temp_labels, temp_samples[i], boost::is_any_of(",;: ")); // various separators allowed to separate labels
          samples_labels_.push_back(temp_labels);
        }
      }
    }
    
    if (samples_labels_.empty())
    {
      vector<String> temp_labels;
      temp_labels.push_back("no_label");
      samples_labels_.push_back(temp_labels);
    }

  }

  String MultiplexMassPatternList::getLabels() const
  {
    return labels_;
  }

}
