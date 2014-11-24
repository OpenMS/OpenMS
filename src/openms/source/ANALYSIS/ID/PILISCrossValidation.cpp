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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PILISCrossValidation.h>

#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace std;

namespace OpenMS
{
  PILISCrossValidation::Peptide::Peptide() :
    charge(0)
  {
  }

  PILISCrossValidation::Peptide::Peptide(const Peptide & rhs) :
    sequence(rhs.sequence),
    charge(rhs.charge),
    spec(rhs.spec),
    hits(rhs.hits)
  {
  }

  PILISCrossValidation::Peptide::~Peptide()
  {
  }

  PILISCrossValidation::Peptide & PILISCrossValidation::Peptide::operator=(const PILISCrossValidation::Peptide & rhs)
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

  bool PILISCrossValidation::Peptide::operator<(const PILISCrossValidation::Peptide & peptide) const
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

  PILISCrossValidation::Option::Option() :
    type(INT),
    int_min(0),
    int_max(0),
    int_stepsize(0),
    dbl_min(0),
    dbl_max(0),
    dbl_stepsize(0)
  {
  }

  PILISCrossValidation::Option::Option(const PILISCrossValidation::Option & rhs) :
    type(rhs.type),
    int_min(rhs.int_min),
    int_max(rhs.int_max),
    int_stepsize(rhs.int_stepsize),
    dbl_min(rhs.dbl_min),
    dbl_max(rhs.dbl_max),
    dbl_stepsize(rhs.dbl_stepsize)
  {
  }

  PILISCrossValidation::Option::Option(Type t, double min, double max, double stepsize)
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

  PILISCrossValidation::Option & PILISCrossValidation::Option::operator=(const PILISCrossValidation::Option & rhs)
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

  PILISCrossValidation::PILISCrossValidation() :
    DefaultParamHandler("PILISCrossValidation")
  {
    defaults_.setValue("nfold", 10, "Number of partitions to use for cross validation");
    defaults_.setValue("optimization_method", "tophit_against_all_others", "Scoring method used for optimization");
    defaults_.setValidStrings("optimization_method", ListUtils::create<String>("tophit_against_all_others,only_top_hit,top_n_ions,top_n_ions_by"));
    defaults_.setValue("compare_function", "SpectrumAlignmentScore", "Spectra scoring function to use");
    defaults_.setValidStrings("compare_function", ListUtils::create<String>("SpectrumAlignmentScore,ZhangSimilarityScore"));
    defaults_.setValue("num_top_peaks", 2, "Number of highest abundant peaks to consider with top_n_ion and top_n_ions_by optimization_methods");
    defaults_.setValue("min_intensity", 0.30, "Min relative intensity of highest abundant peaks to consider in top_n_ions_by");
    defaults_.setValue("fragment_mass_tolerance", 0.5, "Fragment mass tolerance, mainly used in compare function.");
    defaults_.setValue("normalize_to_TIC", "true", "Whether the spectra should be normalized to TIC before training, to max of one otherwise.");
    defaults_.setValidStrings("normalize_to_TIC", ListUtils::create<String>("true,false"));
    defaultsToParam_();
  }

  PILISCrossValidation::PILISCrossValidation(const PILISCrossValidation & rhs) :
    DefaultParamHandler(rhs),
    cv_options_(rhs.cv_options_)
  {
  }

  PILISCrossValidation::~PILISCrossValidation()
  {
  }

  PILISCrossValidation & PILISCrossValidation::operator=(const PILISCrossValidation & rhs)
  {
    if (this != &rhs)
    {
      cv_options_ = rhs.cv_options_;
    }
    return *this;
  }

  void PILISCrossValidation::partition_(vector<vector<Peptide> > & parts, const vector<Peptide> & source)
  {
    Size nfold(param_.getValue("nfold"));
    parts = vector<vector<Peptide> >(nfold);

    Size count = 0;
    set<SignedSize> used_numbers;
    Size part = 0;
    while (count != source.size())
    {
      // get new random variable [0, source.size()]
      Size r = (Size)(double(rand()) / double(RAND_MAX) * (source.size()));

      if (used_numbers.find(r) == used_numbers.end())
      {
        ++count;
        used_numbers.insert(r);
        parts[part++].push_back(source[r]);

        if (part == nfold)
        {
          part = 0;
        }
      }
    }

    return;
  }

// generate parameter permutations
  void PILISCrossValidation::generateParameters_(const Param & param, const Map<String, Option> & options, vector<Param> & parameters)
  {
    if (options.empty())
    {
      parameters.push_back(param);
      return;
    }
    for (Map<String, Option>::ConstIterator it = options.begin(); it != options.end(); ++it)
    {
      Map<String, Option> new_options = options;
      new_options.erase(new_options.find(it->first));
      Param new_param = param;

      if (it->second.type == Option::DOUBLE)
      {
        double dbl_min(it->second.dbl_min), dbl_max(it->second.dbl_max);
        if (dbl_min > dbl_max)
        {
          cerr << "PILISCrossValidation: " << it->first << " min-value > max-value! (" << dbl_min << ", " << dbl_max << ")" << endl;
        }
        for (double value = dbl_min; value <= dbl_max; value += it->second.dbl_stepsize)
        {
          new_param.setValue(it->first, value);
          generateParameters_(new_param, new_options, parameters);
        }
        continue;
      }
      if (it->second.type == Option::INT)
      {
        Int int_min(it->second.int_min), int_max(it->second.int_max);
        if (int_min > int_max)
        {
          cerr << "PILISCrossValidation: " << it->first << " min-value > max-value! (" << int_min << ", " << int_max << ")" << endl;
        }
        for (Int value = int_min; value <= int_max; value += it->second.int_stepsize)
        {
          new_param.setValue(it->first, value);
          generateParameters_(new_param, new_options, parameters);
        }
      }
    }
    return;
  }

  void PILISCrossValidation::apply(Param & PILIS_param, const PILISModel & base_model, const vector<Peptide> & peptides)
  {
    vector<vector<Peptide> > partitions;
    partition_(partitions, peptides);

    SpectrumAlignmentScore sa;
    Param sa_param(sa.getParameters());
    sa_param.setValue("tolerance", (double)PILIS_param.getValue("fragment_mass_tolerance"));
    sa_param.setValue("use_linear_factor", "true");
    sa.setParameters(sa_param);

    Param model_param = base_model.getParameters();
    vector<Param> all_parameters;
    generateParameters_(model_param, cv_options_, all_parameters);

    // gnaaah
    vector<Param> all_new_parameters;
    for (vector<Param>::const_iterator it1 = all_parameters.begin(); it1 != all_parameters.end(); ++it1)
    {
      bool found(false);
      for (vector<Param>::const_iterator it2 = all_new_parameters.begin(); it2 != all_new_parameters.end(); ++it2)
      {
        if (*it1 == *it2)
        {
          found = true;
        }
      }
      if (!found)
      {
        all_new_parameters.push_back(*it1);
      }
    }
    all_parameters = all_new_parameters;
    //cerr << all_parameters.size() << " parameters generated" << endl;
    for (Size i = 0; i != all_parameters.size(); ++i)
    {
      for (Param::ParamIterator it = all_parameters[i].begin(); it != all_parameters[i].end(); ++it)
      {
        //cerr << i+1 << " " << it->name << " " << it->value.toString() << endl;
      }
    }

    Normalizer to_one_normalizer, tic_normalizer;
    Param normalizer_param(to_one_normalizer.getParameters());
    normalizer_param.setValue("method", "to_one");
    to_one_normalizer.setParameters(normalizer_param);
    normalizer_param.setValue("method", "to_TIC");
    tic_normalizer.setParameters(normalizer_param);

    bool normalize_to_TIC(param_.getValue("normalize_to_TIC").toBool());

    // iterate over all parameter setting permutations
    vector<double> scores;
    for (vector<Param>::const_iterator pait = all_parameters.begin(); pait != all_parameters.end(); ++pait)
    {
      cerr << "Param #" << pait - all_parameters.begin() + 1 << "/" << all_parameters.size() << endl;
      PILISModel model = base_model;
      model.setParameters(*pait);
      vector<vector<RichPeakSpectrum> > exp_spectra;
      vector<vector<vector<RichPeakSpectrum> > > sim_spectra;

      // perform cross validation
      for (Size i = 0; i != partitions.size(); ++i)
      {
        cerr << "Partition #" << i + 1 << endl;
        for (Size j = 0; j != partitions.size(); ++j)
        {
          for (vector<Peptide>::const_iterator it = partitions[j].begin(); it != partitions[j].end(); ++it)
          {
            if (i != j)
            {
              RichPeakSpectrum spec = it->spec;
              spec.setMSLevel(2);

              if (normalize_to_TIC)
              {
                tic_normalizer.filterSpectrum(spec);
              }
              else
              {
                to_one_normalizer.filterSpectrum(spec);
              }
              //cerr << "training with " << it->sequence << " " << it->charge << "...";
              model.train(it->spec, it->sequence, it->charge);
              //cerr << "ended" << endl;
            }
          }
        }

        vector<vector<RichPeakSpectrum> > sim_spectra_part;
        vector<RichPeakSpectrum> exp_spectra_part;
        // evaluate this setting
        for (vector<Peptide>::const_iterator it = partitions[i].begin(); it != partitions[i].end(); ++it)
        {
          RichPeakSpectrum exp_spec = it->spec;
          to_one_normalizer.filterSpectrum(exp_spec);
          exp_spec.sortByPosition();
          exp_spectra_part.push_back(exp_spec);

          //cerr << "Evalulating...(#peptides=" << it->hits.size() << ") ";
          // for all hits
          vector<RichPeakSpectrum> sim_spectra_hits;
          for (vector<PeptideHit>::const_iterator pit = it->hits.begin(); pit != it->hits.end(); ++pit)
          {
            RichPeakSpectrum sim_spec;
            //cerr << "simulate " << pit->getSequence() << " " << pit->getCharge() << "..";
            model.getSpectrum(sim_spec, pit->getSequence(), pit->getCharge());

            Precursor prec;
            prec.setCharge(pit->getCharge());
            prec.setPosition((pit->getSequence().getMonoWeight() + (double)pit->getCharge()) / (double)pit->getCharge());
            sim_spec.getPrecursors().push_back(prec);

            //sim_spec.getPrecursorPeak().setCharge(pit->getCharge());
            //sim_spec.getPrecursorPeak().setPosition((pit->getSequence().getMonoWeight() + (double)pit->getCharge())/(double)pit->getCharge());

            to_one_normalizer.filterSpectrum(sim_spec);
            sim_spec.sortByPosition();
            sim_spectra_hits.push_back(sim_spec);
          }
          //cerr << "ended" << endl;
          sim_spectra_part.push_back(sim_spectra_hits);
        }
        sim_spectra.push_back(sim_spectra_part);
        exp_spectra.push_back(exp_spectra_part);
      }

      scores.push_back(scoreHits(sim_spectra, exp_spectra));
    }

    Size best_param_pos(0);
    double max_score(0);
    Size pos(0);
    for (vector<double>::const_iterator it = scores.begin(); it != scores.end(); ++it, ++pos)
    {
      if (*it > max_score)
      {
        max_score = *it;
        best_param_pos = pos;
      }
    }


    cerr << "Best parameters (score=" << max_score << ")" << endl;
    for (Param::ParamIterator it = all_parameters[best_param_pos].begin(); it != all_parameters[best_param_pos].end(); ++it)
    {
      cerr << " " << it->name << " " << it->value.toString() << endl;
    }
    PILIS_param = all_parameters[best_param_pos];

    return;
  }

  double PILISCrossValidation::scoreHits(const vector<vector<vector<RichPeakSpectrum> > > & sim_spectra, const vector<vector<RichPeakSpectrum> > & exp_spectra)
  {
    String optimization_method = param_.getValue("optimization_method");

    if (optimization_method == "tophit_against_all_others")
    {
      // consider all against all
      vector<double> top_scores, non_top_scores;
      for (Size i = 0; i != sim_spectra.size(); ++i)
      {
        for (Size j = 0; j != sim_spectra[i].size(); ++j)
        {
          for (Size k = 0; k != sim_spectra[i][j].size(); ++k)
          {
            double score = scoreSpectra_(sim_spectra[i][j][k], exp_spectra[i][j]);
            if (k == 0)
            {
              top_scores.push_back(score);
            }
            else
            {
              non_top_scores.push_back(score);
            }
          }
        }
      }

      double sum = 0;
      for (vector<double>::const_iterator it1 = top_scores.begin(); it1 != top_scores.end(); ++it1)
      {
        for (vector<double>::const_iterator it2 = non_top_scores.begin(); it2 != non_top_scores.end(); ++it2)
        {
          sum += *it1 - *it2;
        }
      }
      double score = sum / (double)(top_scores.size() * non_top_scores.size());
      cerr << "Avg. score-diff for param: " << score << endl;
      return score;
    }

    if (optimization_method == "only_top_hit")
    {
      // consider only top-scores
      vector<double> top_scores;
      for (Size i = 0; i != sim_spectra.size(); ++i)
      {
        for (Size j = 0; j != sim_spectra[i].size(); ++j)
        {
          if (!sim_spectra[i][j].empty())
          {
            double score = scoreSpectra_(sim_spectra[i][j][0], exp_spectra[i][j]);
            top_scores.push_back(score);
          }
        }
      }
      double sum(0);
      for (vector<double>::const_iterator it = top_scores.begin(); it != top_scores.end(); ++it)
      {
        sum += *it;
      }
      double score(sum / (double)top_scores.size());
      cerr << "Avg. score for param: " << score << endl;
      return score;
    }

    if (optimization_method == "top_n_ions")
    {
      MRMFragmentSelection fragment_selection;
      Size num_correct_topn(0), num_all(0);
      for (Size i = 0; i != sim_spectra.size(); ++i)
      {
        for (Size j = 0; j != sim_spectra[i].size(); ++j)
        {
          ++num_all;
          if (!sim_spectra[i][j].empty())
          {
            vector<RichPeak1D> exp_highest_peak, sim_highest_peak;
            fragment_selection.selectFragments(exp_highest_peak, exp_spectra[i][j]);
            fragment_selection.selectFragments(sim_highest_peak, sim_spectra[i][j][0]);

            double fragment_mass_tolerance((double)param_.getValue("fragment_mass_tolerance"));
            bool has_topn(false);
            for (Size ii = 0; ii < exp_highest_peak.size(); ++ii)
            {
              for (Size jj = ii; jj < sim_highest_peak.size(); ++jj)
              {
                if (fabs(exp_highest_peak[ii].getMZ() - sim_highest_peak[jj].getMZ()) < fragment_mass_tolerance)
                {
                  ++num_correct_topn;
                  has_topn = true;
                  break;
                }
              }
              if (has_topn)
              {
                break;
              }
            }
          }
        }
      }
      double score((double)num_correct_topn / (double)num_all);
      cerr << "Avg. score in top " << param_.getValue("num_top_peaks") << ": " << score << endl;
      return score;
    }

    if (optimization_method == "top_n_ions_by")
    {
      MRMFragmentSelection fragment_selection;
      double fragment_mass_tolerance((double)param_.getValue("fragment_mass_tolerance"));
      double min_intensity((double)param_.getValue("min_intensity"));

      TheoreticalSpectrumGenerator tsg;
      Param tsg_param(tsg.getParameters());
      tsg_param.setValue("add_metainfo", "true");
      tsg.setParameters(tsg_param);

      SpectrumAlignment aligner;
      Param aligner_param(aligner.getParameters());
      aligner_param.setValue("tolerance", fragment_mass_tolerance);
      aligner.setParameters(aligner_param);

      Size num_correct_topn(0), num_all(0);
      for (Size i = 0; i != sim_spectra.size(); ++i)
      {
        for (Size j = 0; j != sim_spectra[i].size(); ++j)
        {
          ++num_all;
          if (!sim_spectra[i][j].empty())
          {
            vector<RichPeak1D> sim_highest_peak;
            //cerr << "Peptide: " << exp_spectra[i][j].getPeptideIdentifications().begin()->getHits().begin()->getSequence() << " "
            //<< exp_spectra[i][j].getPeptideIdentifications().begin()->getHits().begin()->getCharge() << endl;
            cerr << "Sim-highest peaks: (sim-spectrum size=" << sim_spectra[i][j][0].size() << ") ";
            /*
            for (Size p = 0; p != sim_spectra[i][j][0].size(); ++p)
            {
                cerr << sim_spectra[i][j][0][p].getMZ() << " " << sim_spectra[i][j][0][p].getIntensity() << " " << sim_spectra[i][j][0][p].getMetaValue("IonName") << endl;
            }*/
            fragment_selection.selectFragments(sim_highest_peak, sim_spectra[i][j][0]);
            cerr << sim_highest_peak.size();
            cerr << endl;

            // normalize the exp_spectrum to the highest peaks which could be picked
            RichPeakSpectrum exp_spec(exp_spectra[i][j]);
            RichPeakSpectrum theo_spec;
            const AASequence & peptide(exp_spec.getPeptideIdentifications().begin()->getHits().begin()->getSequence());
            tsg.addPeaks(theo_spec, peptide, Residue::BIon, 1);
            //tsg.addPeaks(theo_spec, peptide, Residue::BIon, 2); // TODO check which charge states are allowed
            tsg.addPeaks(theo_spec, peptide, Residue::YIon, 1);
            //tsg.addPeaks(theo_spec, peptide, Residue::YIon, 2);
            theo_spec.sortByPosition();
            vector<pair<Size, Size> > alignment;
            aligner.getSpectrumAlignment(alignment, exp_spec, theo_spec);
            vector<RichPeak1D> exp_highest_peak;

            for (vector<pair<Size, Size> >::const_iterator it = alignment.begin(); it != alignment.end(); ++it)
            {
              exp_spec[it->first].setMetaValue("IonName", (String)theo_spec[it->second].getMetaValue("IonName"));
              exp_highest_peak.push_back(exp_spec[it->first]);
            }

            //cerr << "Exp-highest peaks: ";
            fragment_selection.selectFragments(exp_highest_peak, exp_spec);
            //cerr << endl;

            // normalize the exp spectrum to max possible intensity to be selected
            double max_exp_int(0);
            for (vector<RichPeak1D>::const_iterator it = exp_highest_peak.begin(); it != exp_highest_peak.end(); ++it)
            {
              if (it->getIntensity() > max_exp_int)
              {
                max_exp_int = it->getIntensity();
              }
            }

            for (vector<RichPeak1D>::iterator it = exp_highest_peak.begin(); it != exp_highest_peak.end(); ++it)
            {
              it->setIntensity(it->getIntensity() / max_exp_int);
            }

            bool has_topn(false);
            for (Size ii = 0; ii < exp_highest_peak.size(); ++ii)
            {
              for (Size jj = 0; jj < sim_highest_peak.size(); ++jj)
              {
                //cerr << ii << " " << jj << " " <<  exp_highest_peak[ii].getMZ() << " " << sim_highest_peak[jj].getMZ() << " " << exp_highest_peak[ii].getIntensity() << endl;
                if (fabs(exp_highest_peak[ii].getMZ() - sim_highest_peak[jj].getMZ()) < fragment_mass_tolerance
                   && exp_highest_peak[ii].getIntensity() >= min_intensity)
                {
                  cerr << "Found: " << exp_spec.getPeptideIdentifications().begin()->getHits().begin()->getSequence()
                       << " " << exp_spec.getPeptideIdentifications().begin()->getHits().begin()->getCharge()
                       << " exp m/z=" << exp_highest_peak[ii].getMZ()  << ", sim m/z=" << sim_highest_peak[jj].getMZ()
                       << " exp int=" << exp_highest_peak[ii].getIntensity() <<  ", sim m/z=" << sim_highest_peak[jj].getIntensity() << endl;
                  ++num_correct_topn;
                  has_topn = true;
                  break;
                }
              }
              if (has_topn)
              {
                break;
              }
            }
          }
        }
      }
      double score((double)num_correct_topn / (double)num_all);
      cerr << "Avg. score in top " << param_.getValue("num_top_peaks") << " with a least " << min_intensity * 100.0 << "% intensity: " << score << endl;
      return score;
    }

    cerr << "PILISCrossValidation: unknown optimization_method: " << optimization_method << endl;
    return 0;
  }

  double PILISCrossValidation::scoreSpectra_(const RichPeakSpectrum & spec1, const RichPeakSpectrum & spec2)
  {
    PeakSpectrum s1, s2;
    for (RichPeakSpectrum::ConstIterator piit = spec1.begin(); piit != spec1.end(); ++piit)
    {
      Peak1D p;
      p.setPosition(piit->getPosition()[0]);
      p.setIntensity(piit->getIntensity());
      s1.push_back(p);
    }
    for (RichPeakSpectrum::ConstIterator piit = spec2.begin(); piit != spec2.end(); ++piit)
    {
      Peak1D p;
      p.setPosition(piit->getPosition()[0]);
      p.setIntensity(piit->getIntensity());
      s2.push_back(p);
    }
    return (*pscf_)(s1, s2);
  }

  void PILISCrossValidation::updateMembers_()
  {
    pscf_ = Factory<PeakSpectrumCompareFunctor>::create(param_.getValue("compare_function"));
    Param compare_param(pscf_->getParameters());
    if (compare_param.exists("tolerance"))
    {
      compare_param.setValue("tolerance", (double)param_.getValue("fragment_mass_tolerance"));
      pscf_->setParameters(compare_param);
    }
    return;
  }

}
