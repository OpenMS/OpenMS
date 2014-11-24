// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/SuffixArraySeqan.h>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <stack>
#include <typeinfo>

#include <OpenMS/CHEMISTRY/ModifierRep.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Exception.h>

//using namespace seqan; // do not use since seqan defines a class set, which makes std::set ambiguous
using namespace std;


namespace OpenMS
{

  // constructor
  SuffixArraySeqan::SuffixArraySeqan(const String & st, const String & sa_file_name, const WeightWrapper::WEIGHTMODE weight_mode) :
    WeightWrapper(weight_mode),
    s_(st),
    number_of_modifications_(0),
    use_tags_(false),
    tol_(0.5)
  {
    //tol_ = 0.5;
    //use_tags_ = false;
    //number_of_modifications_ = 0;
    if (st[0] != '$')
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "String has to start with empyt string ($)", "");
    }
    if (st[st.length() - 1] != '$')
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "String has to end with $", "");
    }
    //creating array with aminoacid masses

    ResidueDB * rdb = ResidueDB::getInstance();

    char aa[] = "ARNDCEQGHILKMFPSTWYV";

    for (Size z = 0; z < 255; ++z)
    {
      masse_[z] = 0;
    }

    for (Size z = 0; z < strlen(aa); ++z)
    {
      const Residue * r = rdb->getResidue(aa[z]);
      masse_[(int)aa[z]] = this->getWeight(*r, Residue::Internal);
    }

    if (sa_file_name != "")
    {
      ifstream file;
      file.open((sa_file_name + ".txt").c_str());

      if (!file.is_open())
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "");
      }
      open(sa_file_name);
    }
    else
    {
      //create

      TIndex index(s_.c_str());
      index_ = index;
    }

    it_ = new TIter(index_);
    // TODO was: it_ = new Iter<TIndex, VSTree< TopDown< ParentLinks<Preorder> > > > (index_);

  }

  SuffixArraySeqan::SuffixArraySeqan(const SuffixArraySeqan & source) :
    SuffixArray(source),
    WeightWrapper(source),
    index_(source.index_),
    it_(source.it_),
    s_(source.s_),
    number_of_modifications_(source.number_of_modifications_),
    tags_(source.tags_),
    use_tags_(source.use_tags_),
    tol_(source.tol_)
  {
    for (Size i = 0; i < 255; ++i)
    {
      masse_[i] = source.masse_[i];
    }
  }

  bool SuffixArraySeqan::isDigestingEnd(const char, const char) const
  {
    return true;
  }

  bool SuffixArraySeqan::save(const String & file_name)
  {
    if (!seqan::save(index_, file_name.c_str()))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name + ".txt"));
    }
    return true;
  }

  bool SuffixArraySeqan::open(const String & file_name)
  {
    if (!seqan::open(index_, file_name.c_str()))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name + ".txt"));
    }

    if (!seqan::indexSupplied(index_, seqan::EsaSA()) ||
        !seqan::indexSupplied(index_, seqan::EsaLcp()) ||
        !seqan::indexSupplied(index_, seqan::EsaChildtab()))
    {
      //cout<<"creating index " << endl;

      seqan::indexRequire(index_, seqan::EsaSA());
      seqan::indexRequire(index_, seqan::EsaLcp());
      seqan::indexRequire(index_, seqan::EsaChildtab());
      seqan::save(index_, file_name.c_str());
    }

    return true;
  }

  //destructor
  SuffixArraySeqan::~SuffixArraySeqan()
  {

  }

  String SuffixArraySeqan::toString()
  {
    return "";
  }

  void SuffixArraySeqan::setTags(const vector<String> & tags)
  {
    tags_ = tags;
    use_tags_ = true;
  }

  const vector<String> & SuffixArraySeqan::getTags()
  {
    return tags_;
  }

  void SuffixArraySeqan::setUseTags(bool use_tags)
  {
    use_tags_ = use_tags;
    if (tags_.empty())
    {
      use_tags_ = false;
    }
  }

  bool SuffixArraySeqan::getUseTags()
  {
    return use_tags_;
  }

  void SuffixArraySeqan::setNumberOfModifications(Size number_of_mods)
  {
    number_of_modifications_ = number_of_mods;
  }

  Size SuffixArraySeqan::getNumberOfModifications()
  {
    return number_of_modifications_;
  }

  SignedSize SuffixArraySeqan::findFirst_(const vector<double> & spec, double & m, SignedSize start, SignedSize  end)
  {

    if (end - start <= 1)
      return (spec[start] < m - tol_) ? end : start;

    SignedSize middle = ((end - start) / 2) + start;

    if (spec[middle] < m - tol_)
    {
      return findFirst_(spec, m, middle, end);
    }
    if (spec[middle] > m + tol_)
    {
      return findFirst_(spec, m, start, middle);
    }
    while (middle >= 0 && spec[middle] >= m - tol_)
    {
      --middle;
    }
    return middle + 1;
  }

  SignedSize SuffixArraySeqan::findFirst_(const vector<double> & spec, double & m)
  {
    return findFirst_(spec, m, 0, spec.size() - 1);
  }

  // finds all occurrences of a given spectrum
  void SuffixArraySeqan::findSpec(vector<vector<pair<pair<SignedSize, SignedSize>, double> > > & candidates, const vector<double> & spec)
  {
    if (spec.empty())
    {
      return;
    }
    //check if spectrum is sorted
    //time_t t1 (time(NULL));
    for (Size i = 1; i < spec.size(); ++i)
    {
      if (spec.at(i - 1) > spec.at(i))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Spectrum has to be sorted ascendingly", "");
      }
    }
    // modification
    ModifierRep modifier;

    modifier.setNumberOfModifications(number_of_modifications_);

    Size number_of_posible_mods = modifier.getMaxModificationMasses();

    // all tags
    vector<Size> tag_indices;
    if (use_tags_)
    {
      for (Size i = 0; i < tags_.size(); ++i)
      {
        it_ = new TIter(index_);

        seqan::String<char> s(tags_.at(i).c_str());
        goDown(*it_, s);
        seqan::String<Size> occs = getOccurrences(*it_);
        for (Size j = 0; j < length(occs); ++j)
        {
          tag_indices.push_back(occs[j]);
        }
      }
      sort(tag_indices.begin(), tag_indices.end());
    }

    it_ = new TIter(index_);

    // preparing result vector
    //vector<vector<pair<pair<SignedSize, SignedSize>,double> > > res;
    for (Size i = 0; i < spec.size(); i++)
    {
      vector<pair<pair<SignedSize, SignedSize>, double> > v;
      candidates.push_back(v);
    }
    double mmax = spec.back();
    stack<double> allm;
    stack<map<double, SignedSize> > history;
    history.push(map<double, SignedSize>());

    double m = getWeight(EmpiricalFormula("H2O"));
    goNext_(*it_, m, allm, history);
    //goNextSubTree(*it_);
    SignedSize nres = 0;

    SignedSize steps4 = 0;
    //iterating over suffix array
    while (!atEnd(*it_))
    {
      SignedSize start_index_in_text = getOccurrence(*it_);
      char start_char = s_[start_index_in_text];
      char next_char = ((Size)start_index_in_text == seqan::length(s_) - 1) ? 'R' : s_[start_index_in_text + 1];

      map<double, SignedSize> modification_map(history.top());

      // because of searching only candidates generated by a specific digestion
      // enzyme we are either looking for candidates with the specific start
      // pattern or for the start indicated by the separator character
      if (start_index_in_text == 0 || start_char == '$' || isDigestingEnd(start_char, next_char))
      {
        SignedSize edge_length = length(parentEdgeLabel(*it_));
        SignedSize length_till_node = length(representative(*it_)) - edge_length;

        double subm = 0;
        double mm = 0;

        // br indicates if break was used
        bool br = false;

        for (SignedSize i = 0; i < edge_length; ++i)
        {
          // actual character is start_index_in_text + length_till_node + how many steps i walked
          char cc = s_[start_index_in_text + length_till_node + i];
          if (modification_map.size() < number_of_posible_mods)
          {
            modifier.refreshModificationList(modification_map, cc);
          }
          // for the next character we have to keep attention if we are on the
          // end of the string. Therefor the next character has to be set to a
          // amino acid so that the digesting enzyme will cut before this
          // position (i.e. for trypsin everything but P)
          char ccn = ((Size)(length_till_node + i + start_index_in_text + 1) == s_.length() - 1) ? 'R' : s_[length_till_node + i + start_index_in_text + 1];
          subm += masse_[(int)cc];
          mm = m + subm;

          // we always have to substract the mass of the start character (either $ or for trypsin K or R)
          double newm = (mm - masse_[(SignedSize)start_char]);

          // if we reached the maxmimal mass we can directly skip the sub tree
          if (newm > mmax + tol_)
          {
            allm.push(0);
            history.push(map<double, SignedSize>());
            goNextSubTree_(*it_, m, allm, history);
            br = true;
            break;
          }
          // if we are reaching a separetor character
          if ((i < 1 && length_till_node < 1) ? false : (cc == '$'))
          {
            // either we are not using tags or we have already seen one of the tags
            if (!use_tags_ ||
                binary_search(tag_indices.begin(), tag_indices.end(), Size(length_till_node + i + start_index_in_text - 2) /*, IntsInRangeLess(length_till_node+i-2,length_till_node+i+start_index_in_text-2)*/)
                )
            {
              // if the mass is in spectrum but only if the digesting enzyme
              // will not cut after the last character (because if it does the
              // candidate was already added to the result the step before)
              //
              // for every modification mass within the modification_map we
              // check if the mass + modification_mass is in the given
              // spectrum. If it is we will add it to a vector of masses and
              // adding this to result vector
              if ((!isDigestingEnd(s_[length_till_node + i + start_index_in_text - 1], cc)))
              {
                vector<double> found_masses;
                if (binary_search(spec.begin(), spec.end(), newm, FloatsWithTolLess(tol_)))
                {
                  found_masses.push_back(0);
                }
                // if the mass is in spectrum we will add the entry to all matching masses
                map<double, SignedSize>::iterator it;
                for (it = modification_map.begin(); it != modification_map.end(); ++it)
                {
                  if (binary_search(spec.begin(), spec.end(), newm + it->first, FloatsWithTolLess(tol_)))
                  {
                    found_masses.push_back(it->first);
                  }
                }
                for (Size o = 0; o < found_masses.size(); o++)
                {
                  double mass_with_mods = newm + found_masses.at(o);
                  ++steps4;
                  // getting all occurences and adding the to the specific masses
                  seqan::String<Size> occ = getOccurrences(*it_);

                  nres += length(occ);
                  for (Size k = 0; k < length(occ); ++k)
                  {
                    Size first_occ = findFirst_(spec, mass_with_mods);
                    pair<pair<SignedSize, SignedSize>, double> p(pair<SignedSize, SignedSize>(occ[k] + 1, length_till_node + i - 1), found_masses.at(o));

                    while (first_occ < spec.size() && spec.at(first_occ) <= mass_with_mods + tol_)
                    {
                      candidates.at(first_occ).push_back(p);

                      ++first_occ;
                    }
                  }
                }
              }
            }
            // because of having reached a separator we can skip the sub tree
            history.push(map<double, SignedSize>());
            allm.push(0);
            goNextSubTree_(*it_, m, allm, history);
            br = true;
            break;
          }


          // if we reached a digesting site
          // the case that i==(edge_length-1) means we are at a node. If we are
          // at a node we cannot just look at one following character but
          // instead we must look at the next of every outgoing edge and
          // deciding whether this is a digesting site
          if (i == (edge_length - 1) || isDigestingEnd(cc, ccn))
          {
            double local_newm = (mm - masse_[(int)start_char]);
            // if the mass is in the spectrum
            if (!use_tags_ ||
                binary_search(tag_indices.begin(), tag_indices.end(), Size(length_till_node + i + start_index_in_text - 2) /*,
                                IntsInRangeLess(length_till_node+i-2,length_till_node+i+start_index_in_text-2)*/)
                )
            {
              // for every modification mass within the modification_map we
              // check if the mass + modification_mass is in the given
              // spectrum. If it is we will add it to a vector of masses and
              // adding this to result vector
              vector<double> found_masses;
              if (binary_search(spec.begin(), spec.end(), local_newm, FloatsWithTolLess(tol_)))
              {
                found_masses.push_back(0);
              }
              // if the mass is in spectrum we will add the entry to all matching masses
              map<double, SignedSize>::iterator it;
              for (it = modification_map.begin(); it != modification_map.end(); ++it)
              {
                if (binary_search(spec.begin(), spec.end(), local_newm + (double)it->first, FloatsWithTolLess(tol_)))
                {
                  found_masses.push_back(it->first);
                }
              }
              for (Size o = 0; o < found_masses.size(); o++)
              {
                double mass_with_mods = local_newm + found_masses.at(o);
                //if (binary_search(spec.begin(),spec.end(),(local_newm), FloatsWithTolLess(tol_))){
                // getting all occurrences and adding the to the specific masses
                ++steps4;

                seqan::String<Size> occ = getOccurrences(*it_);
                //nres+=length(occ);
                for (Size k = 0; k < length(occ); k++)
                {
                  if (i < (edge_length - 1) || (isDigestingEnd(s_[occ[k] + length_till_node + i], s_[occ[k] + 1 + length_till_node + i])))
                  {
                    Size first_occ = findFirst_(spec, mass_with_mods);
                    pair<pair<SignedSize, SignedSize>, double> p(pair<SignedSize, SignedSize>(occ[k] + 1, length_till_node + i), found_masses.at(o));
                    ++nres;
                    while (first_occ < spec.size() && spec.at(first_occ) <= mass_with_mods + tol_)
                    {
                      candidates.at(first_occ).push_back(p);
                      ++first_occ;
                    }
                  }
                }
              }
            }
          }
        }

        // if we breaked before we are already at the next sub tree so we will not use goNext
        if (!br)
        {
          m = mm;
          //because of the on-the-fly mass update the updated mass differs from actual mass, so from time to time we can correct the actual mass
          //TODO: why would that be?! (Andreas)
          if (steps4 > 1000)
          {
            double mpart = getWeight(EmpiricalFormula("H2O"));
            seqan::String<char> seq = representative(*it_);
            for (Size w = 0; w < length(seq); ++w)
            {
              mpart += masse_[(int)seq[w]];
            }
            m = mpart;
            steps4 = 0;
          }
          history.push(map<double, SignedSize>(modification_map));
          allm.push(subm);
          goNext_(*it_, m, allm, history);
        }

      }
      else
      {
        history.push(map<double, SignedSize>());
        allm.push(0);
        goNextSubTree_(*it_, m, allm, history);
      }
    }

    //time_t t2 (time(NULL));
    //cout <<"number of hits: " << nres<<endl;
    //cout <<"used time: "<< t2-t1<<endl;
    return;
  }

  void SuffixArraySeqan::setTolerance(double t)
  {
    if (t < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Tolerance value must not be negative", (String)t);
    }
    tol_ = t;
  }

  double SuffixArraySeqan::getTolerance() const
  {
    return tol_;
  }

  void SuffixArraySeqan::printStatistic()
  {
    it_ = new TIter(index_);

    vector<pair<SignedSize, SignedSize> > out_number;
    vector<pair<SignedSize, SignedSize> > edge_length;
    vector<SignedSize> leafe_depth;
    //goNext(*it_);
    goNextSubTree_(*it_);
    parseTree_(*it_, out_number, edge_length, leafe_depth);
    for (Size i = 0; i < leafe_depth.size(); i++)
    {
      cout << leafe_depth.at(i) << ",";
    }
    cout << endl;
    for (Size i = 0; i < out_number.size(); i++)
    {
      cout << "(" << out_number.at(i).first << "," << out_number.at(i).second << ") ; ";
    }
    cout << endl;
    for (Size i = 0; i < edge_length.size(); i++)
    {
      cout << "(" << edge_length.at(i).first << "," << edge_length.at(i).second << ") ; ";
    }
    cout << endl;
  }

}
