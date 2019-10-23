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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>

//#define DAC_DEBUG
//#define ESTIMATE_PRECURSOR_DEBUG

//#define ETD_SPECTRUM_DEBUG

//#define WRITE_SCORED_SPEC
//#define REDUCE_PERMUTS_DEBUG

//#define SPIKE_IN

using namespace std;

namespace OpenMS
{
  CompNovoIdentification::CompNovoIdentification() :
    CompNovoIdentificationBase()
  {
  }

  CompNovoIdentification::CompNovoIdentification(const CompNovoIdentification & rhs) :
    CompNovoIdentificationBase(rhs)
  {
  }

  CompNovoIdentification & CompNovoIdentification::operator=(const CompNovoIdentification & rhs)
  {
    if (this != &rhs)
    {
      CompNovoIdentificationBase::operator=(rhs);
    }
    return *this;
  }

  CompNovoIdentification::~CompNovoIdentification()
  {
  }

  void CompNovoIdentification::getIdentifications(vector<PeptideIdentification> & pep_ids, const PeakMap & exp)
  {
    Size count(1);
    for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it, ++count)
    {
      //cerr << count << "/" << exp.size() / 2 << endl;

      PeptideIdentification id;
      PeakSpectrum CID_spec(*it);
      double cid_rt(it->getRT());
      double cid_mz(0);
      if (!it->getPrecursors().empty())
      {
        cid_mz = it->getPrecursors().begin()->getMZ();
      }

      if (it->getPrecursors().empty() || cid_mz == 0)
      {
        cerr << "CompNovoIdentification: Spectrum id=\"" << it->getNativeID() << "\" at RT=" << cid_rt << " does not have valid precursor information." << endl;
        continue;
      }
      id.setRT(cid_rt);
      id.setMZ(cid_mz);

      if ((it + 1) != exp.end() && !(it + 1)->getPrecursors().empty())
      {
        double etd_rt = (it + 1)->getRT();
        double etd_mz = (it + 1)->getPrecursors().begin()->getMZ();

        if (fabs(etd_rt - cid_rt) < 10 &&         // RT distance is not too large
            fabs(etd_mz - cid_mz) < 0.01)             // same precursor used
        {
          PeakSpectrum ETD_spec(*(++it));

          subspec_to_sequences_.clear();
          permute_cache_.clear();

          getIdentification(id, CID_spec, ETD_spec);
          //cerr << "size_of id=" << id.getHits().size() << endl;
          pep_ids.push_back(id);
        }
      }
    }
  }

  void CompNovoIdentification::getIdentification(PeptideIdentification & id, const PeakSpectrum & CID_spec, const PeakSpectrum & ETD_spec)
  {
    PeakSpectrum new_CID_spec(CID_spec), new_ETD_spec(ETD_spec);
    windowMower_(new_CID_spec, 0.3, 1);
    windowMower_(new_ETD_spec, 0.3, 1);

    Param zhang_param;
    zhang_param = zhang_.getParameters();
    zhang_param.setValue("tolerance", fragment_mass_tolerance_);
    zhang_param.setValue("use_gaussian_factor", "true");
    zhang_param.setValue("use_linear_factor", "false");
    zhang_.setParameters(zhang_param);


    Normalizer normalizer;
    Param n_param(normalizer.getParameters());
    n_param.setValue("method", "to_one");
    normalizer.setParameters(n_param);

    normalizer.filterSpectrum(new_CID_spec);
    normalizer.filterSpectrum(new_ETD_spec);

    Size charge(0);
    double precursor_weight = 0;

    bool estimate_precursor_mz = param_.getValue("estimate_precursor_mz").toBool();

    if (estimate_precursor_mz)
    {
      precursor_weight = estimatePrecursorWeight_(new_ETD_spec, charge);
    }

    if (precursor_weight == 0 || charge == 0)
    {
      if (CID_spec.getPrecursors().empty())
      {
        cerr << "No precursors found, skipping identification." << endl;
        return;
      }

      if (CID_spec.getPrecursors().begin()->getCharge() != 0)
      {
        charge = CID_spec.getPrecursors().begin()->getCharge();
      }
      else
      {
        cerr << "No charge annotated with precursor, estimating as 2+" << endl;
        charge = 2;
      }
      precursor_weight = CID_spec.getPrecursors().begin()->getMZ() * (double)charge - (double)(charge - 1) * Constants::PROTON_MASS_U;
    }


    if (precursor_weight > 2000.0)
    {
      cerr << "Weight of precursor has been estimated to exceed 2000.0 Da which is the current limit: " << precursor_weight << endl;
      return;
    }

    // now delete all peaks that are right of the estimated precursor weight
    Size peak_counter(0);
    for (PeakSpectrum::ConstIterator it = new_CID_spec.begin(); it != new_CID_spec.end(); ++it, ++peak_counter)
    {
      if (it->getPosition()[0] > precursor_weight)
      {
        break;
      }
    }
    if (peak_counter < new_CID_spec.size())
    {
      new_CID_spec.resize(peak_counter);
    }
    peak_counter = 0;
    for (PeakSpectrum::ConstIterator it = new_ETD_spec.begin(); it != new_ETD_spec.end(); ++it, ++peak_counter)
    {
      if (it->getPosition()[0] > precursor_weight)
      {
        break;
      }
    }
    if (peak_counter < new_ETD_spec.size())
    {
      new_ETD_spec.resize(peak_counter);
    }

    double precursor_mass_tolerance((double)param_.getValue("precursor_mass_tolerance"));

    // delete the precursor peaks from the ETD spec
    PeakSpectrum ETD_copy;
    for (PeakSpectrum::ConstIterator it = new_ETD_spec.begin(); it != new_ETD_spec.end(); ++it, ++peak_counter)
    {
      double pre_pos((precursor_weight + 1.0 * Constants::PROTON_MASS_U) / precursor_mass_tolerance);
      if (fabs(it->getPosition()[0] - pre_pos) > precursor_mass_tolerance)
      {
        ETD_copy.push_back(*it);
      }
    }

    new_ETD_spec = ETD_copy;

    static double oxonium_mass = EmpiricalFormula("H2O+").getMonoWeight();

    Peak1D p;
    p.setIntensity(1.0f);
    p.setPosition(oxonium_mass);

    new_CID_spec.push_back(p);

    p.setPosition(precursor_weight);
    new_CID_spec.push_back(p);

    //cerr << "Estimated precursor weight: " << precursor_weight << ", from file: " << CID_spec.getPrecursors().begin()->getMZ() << endl;

    if (charge == 3)
    {
      // add complement to spectrum
      for (PeakSpectrum::ConstIterator it1 = CID_spec.begin(); it1 != CID_spec.end(); ++it1)
      {
        // get m/z of complement
        double mz_comp = precursor_weight - it1->getPosition()[0] + Constants::PROTON_MASS_U;

        // search if peaks are available that have similar m/z values
        Size count(0);
        bool found(false);
        for (PeakSpectrum::ConstIterator it2 = CID_spec.begin(); it2 != CID_spec.end(); ++it2, ++count)
        {
          if (fabs(mz_comp - it2->getPosition()[0]) < fragment_mass_tolerance_)
          {
            // add peak intensity to corresponding peak in new_CID_spec
            new_CID_spec[count].setIntensity(new_CID_spec[count].getIntensity());
          }
        }
        if (!found)
        {
          // infer this peak
          Peak1D lp;
          lp.setIntensity(it1->getIntensity());
          lp.setPosition(mz_comp);
          new_CID_spec.push_back(lp);
        }
      }
      // add putative b/y-ions from c/z ions

    }

    new_CID_spec.sortByPosition();
    new_ETD_spec.sortByPosition();

    CompNovoIonScoring ion_scoring;
    Param ion_scoring_param(ion_scoring.getParameters());
    ion_scoring_param.setValue("fragment_mass_tolerance", fragment_mass_tolerance_);
    ion_scoring_param.setValue("decomp_weights_precision", decomp_weights_precision_);
    ion_scoring_param.setValue("double_charged_iso_threshold", (double)param_.getValue("double_charged_iso_threshold"));
    ion_scoring_param.setValue("max_isotope_to_score", param_.getValue("max_isotope_to_score"));
    ion_scoring_param.setValue("max_isotope", max_isotope_);
    ion_scoring.setParameters(ion_scoring_param);

    Map<double, IonScore> ion_scores;
    ion_scoring.scoreSpectra(ion_scores, new_CID_spec, new_ETD_spec, precursor_weight, charge);

    new_CID_spec.sortByPosition();
    new_ETD_spec.sortByPosition();

    // check the distances between the ions and reset the max_decomp_weight_ if necessary
    /*
    double max_distance(0);
    for (PeakSpectrum::ConstIterator it1 = new_CID_spec.begin(); it1 != new_CID_spec.end(); ++it1)
    {
        PeakSpectrum::ConstIterator it2 = it1;
        ++it2;

        if (it2 != new_CID_spec.end())
        {
            double dist = it2->getPosition()[0] - it1->getPosition()[0];
            if (dist > max_distance)
            {
                max_distance = dist;
            }
        }
    }

    if (max_distance > max_decomp_weight_)
    {
        param_.setValue("max_decomp_weight", max_distance + 10);
        max_decomp_weight_ = max_distance + 10;
        cerr << "new max_decomp_weight_=" << max_decomp_weight_ << endl;
    }*/


    /*
    cerr << "Size of ion_scores " << ion_scores.size() << endl;
    for (Map<double, IonScore>::const_iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
    {
        cerr << it->first << " " << it->second.score << endl;
    }*/

#ifdef WRITE_SCORED_SPEC
    PeakSpectrum filtered_spec(new_CID_spec);
    filtered_spec.clear();
    for (Map<double, CompNovoIonScoring::IonScore>::const_iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
    {
      Peak1D p;
      p.setIntensity(it->second.score);
      p.setPosition(it->first);
      filtered_spec.push_back(p);
    }
    DTAFile().store("spec_scored.dta", filtered_spec);

    PeakSpectrum test_etd_spec, test_cid_spec;
#define TEST_PEPTIDE "KEELLLPEWILQR"
    getETDSpectrum_(test_etd_spec, TEST_PEPTIDE, 2, 0.0, 0.0);
    getCIDSpectrum_(test_cid_spec, TEST_PEPTIDE, 2, 0.0, 0.0);
    DTAFile().store(TEST_PEPTIDE + String("_etd_test.dta"), test_etd_spec);
    DTAFile().store(TEST_PEPTIDE + String("_cid_test.dta"), test_cid_spec);
#endif

    set<String> sequences;
    getDecompositionsDAC_(sequences, 0, new_CID_spec.size() - 1, precursor_weight, new_CID_spec, new_ETD_spec, ion_scores);

#ifdef SPIKE_IN
    sequences.insert("AFCVDGEGR");
    sequences.insert("APEFAAPWPDFVPR");
    sequences.insert("AVKQFEESQGR");
    sequences.insert("CCTESLVNR");
    sequences.insert("DAFLGSFLYEYSR");
    sequences.insert("DAIPENLPPLTADFAEDK");
    sequences.insert("DDNKVEDIWSFLSK");
    sequences.insert("DDPHACYSTVFDK");
    sequences.insert("DEYELLCLDGSR");
    sequences.insert("DGAESYKELSVLLPNR");
    sequences.insert("DGASCWCVDADGR");
    sequences.insert("DLFIPTCLETGEFAR");
    sequences.insert("DTHKSEIAHR");
    sequences.insert("DVCKNYQEAK");
    sequences.insert("EACFAVEGPK");
    sequences.insert("ECCHGDLLECADDR");
    sequences.insert("EFLGDKFYTVISSLK");
    sequences.insert("EFTPVLQADFQK");
    sequences.insert("ELFLDSGIFQPMLQGR");
    sequences.insert("ETYGDMADCCEK");
    sequences.insert("EVGCPSSSVQEMVSCLR");
    sequences.insert("EYEATLEECCAK");
    sequences.insert("FADLIQSGTFQLHLDSK");
    sequences.insert("FFSASCVPGATIEQK");
    sequences.insert("FLANVSTVLTSK");
    sequences.insert("FLSGSDYAIR");
    sequences.insert("FTASCPPSIK");
    sequences.insert("GAIEWEGIESGSVEQAVAK");
    sequences.insert("GDVAFIQHSTVEENTGGK");
    sequences.insert("GEPPSCAEDQSCPSER");
    sequences.insert("GEYVPTSLTAR");
    sequences.insert("GQEFTITGQKR");
    sequences.insert("GTFAALSELHCDK");
    sequences.insert("HLVDEPQNLIK");
    sequences.insert("HQDCLVTTLQTQPGAVR");
    sequences.insert("HTTVNENAPDQK");
    sequences.insert("ILDCGSPDTEVR");
    sequences.insert("KCPSPCQLQAER");
    sequences.insert("KGTEFTVNDLQGK");
    sequences.insert("KQTALVELLK");
    sequences.insert("KVPQVSTPTLVEVSR");
    sequences.insert("LALQFTTNAKR");
    sequences.insert("LCVLHEKTPVSEK");
    sequences.insert("LFTFHADICTLPDTEK");
    sequences.insert("LGEYGFQNALIVR");
    sequences.insert("LHVDPENFK");
    sequences.insert("LKECCDKPLLEK");
    sequences.insert("LKHLVDEPQNLIK");
    sequences.insert("LKPDPNTLCDEFK");
    sequences.insert("LLGNVLVVVLAR");
    sequences.insert("LLVVYPWTQR");
    sequences.insert("LRVDPVNFK");
    sequences.insert("LTDEELAFPPLSPSR");
    sequences.insert("LVNELTEFAK");
    sequences.insert("MFLSFPTTK");
    sequences.insert("MPCTEDYLSLILNR");
    sequences.insert("NAPYSGYSGAFHCLK");
    sequences.insert("NECFLSHKDDSPDLPK");
    sequences.insert("NEPNKVPACPGSCEEVK");
    sequences.insert("NLQMDDFELLCTDGR");
    sequences.insert("QAGVQAEPSPK");
    sequences.insert("RAPEFAAPWPDFVPR");
    sequences.insert("RHPEYAVSVLLR");
    sequences.insert("RPCFSALTPDETYVPK");
    sequences.insert("RSLLLAPEEGPVSQR");
    sequences.insert("SAFPPEPLLCSVQR");
    sequences.insert("SAGWNIPIGTLLHR");
    sequences.insert("SCWCVDEAGQK");
    sequences.insert("SGNPNYPHEFSR");
    sequences.insert("SHCIAEVEK");
    sequences.insert("SISSGFFECER");
    sequences.insert("SKYLASASTMDHAR");
    sequences.insert("SLHTLFGDELCK");
    sequences.insert("SLLLAPEEGPVSQR");
    sequences.insert("SPPQCSPDGAFRPVQCK");
    sequences.insert("SREGDPLAVYLK");
    sequences.insert("SRQIPQCPTSCER");
    sequences.insert("TAGTPVSIPVCDDSSVK");
    sequences.insert("TCVADESHAGCEK");
    sequences.insert("TQFGCLEGFGR");
    sequences.insert("TVMENFVAFVDK");
    sequences.insert("TYFPHFDLSHGSAQVK");
    sequences.insert("TYMLAFDVNDEK");
    sequences.insert("VDEVGGEALGR");
    sequences.insert("VDLLIGSSQDDGLINR");
    sequences.insert("VEDIWSFLSK");
    sequences.insert("VGGHAAEYGAEALER");
    sequences.insert("VGTRCCTKPESER");
    sequences.insert("VKVDEVGGEALGR");
    sequences.insert("VKVDLLIGSSQDDGLINR");
    sequences.insert("VLDSFSNGMK");
    sequences.insert("VLSAADKGNVK");
    sequences.insert("VPQVSTPTLVEVSR");
    sequences.insert("VTKCCTESLVNR");
    sequences.insert("VVAASDASQDALGCVK");
    sequences.insert("VVAGVANALAHR");
    sequences.insert("YICDNQDTISSK");
    sequences.insert("YLASASTMDHAR");
    sequences.insert("YNGVFQECCQAEDK");
#endif

    SpectrumAlignmentScore spectra_zhang;
    spectra_zhang.setParameters(zhang_param);

    vector<PeptideHit> hits;
    Size missed_cleavages = param_.getValue("missed_cleavages");
    for (set<String>::const_iterator it = sequences.begin(); it != sequences.end(); ++it)
    {

      Size num_missed = countMissedCleavagesTryptic_(*it);
      if (missed_cleavages < num_missed)
      {
        //cerr << "Two many missed cleavages: " << *it << ", found " << num_missed << ", allowed " << missed_cleavages << endl;
        continue;
      }
      PeakSpectrum ETD_sim_spec, CID_sim_spec;
      getETDSpectrum_(ETD_sim_spec, *it, charge);
      getCIDSpectrum_(CID_sim_spec, *it, charge);

      //normalizer.filterSpectrum(ETD_sim_spec);
      //normalizer.filterSpectrum(CID_sim_spec);

      double cid_score = zhang_(CID_sim_spec, CID_spec);
      double etd_score = zhang_(ETD_sim_spec, ETD_spec);

      PeptideHit hit;
      hit.setScore(cid_score + etd_score);

      hit.setSequence(getModifiedAASequence_(*it));
      hit.setCharge((Int)charge);   //TODO unify charge interface: int or size?
      hits.push_back(hit);
      //cerr << getModifiedAASequence_(*it) << " " << cid_score << " " << etd_score << " " << cid_score + etd_score << endl;
    }

    // rescore the top hits
    id.setHits(hits);
    id.assignRanks();

    hits = id.getHits();

    SpectrumAlignmentScore alignment_score;
    Param align_param(alignment_score.getParameters());
    align_param.setValue("tolerance", fragment_mass_tolerance_);
    align_param.setValue("use_linear_factor", "true");
    alignment_score.setParameters(align_param);

    /*
    for (vector<PeptideHit>::iterator it = hits.begin(); it != hits.end(); ++it)
    {
        cerr << "Pre: " << it->getRank() << " " << it->getSequence() << " " << it->getScore() << " " << endl;
    }
    */

    Size number_of_prescoring_hits = param_.getValue("number_of_prescoring_hits");
    if (hits.size() > number_of_prescoring_hits)
    {
      hits.resize(number_of_prescoring_hits);
    }

    for (vector<PeptideHit>::iterator it = hits.begin(); it != hits.end(); ++it)
    {
      PeakSpectrum ETD_sim_spec, CID_sim_spec;
      String mod_string = getModifiedStringFromAASequence_(it->getSequence());
      getETDSpectrum_(ETD_sim_spec, mod_string, charge);
      getCIDSpectrum_(CID_sim_spec, mod_string, charge);

      normalizer.filterSpectrum(ETD_sim_spec);
      normalizer.filterSpectrum(CID_sim_spec);

      double cid_score = alignment_score(CID_sim_spec, CID_spec);
      double etd_score = alignment_score(ETD_sim_spec, ETD_spec);

      /*
      cerr << "Final: " << it->getSequence() << " " << cid_score << " " << etd_score << " " << 2 * cid_score + etd_score << endl;
      */

      it->setScore(cid_score + etd_score);
    }

    id.setHits(hits);
    id.assignRanks();
    hits = id.getHits();

    /*
    for (vector<PeptideHit>::iterator it = hits.begin(); it != hits.end(); ++it)
    {
        cerr << "Fin: " << it->getRank() << " " << it->getSequence() << " " << it->getScore() << " " << endl;
    }
    */

    Size number_of_hits = param_.getValue("number_of_hits");
    if (id.getHits().size() > number_of_hits)
    {
      hits.resize(number_of_hits);
    }

    id.setHits(hits);
    id.assignRanks();
  }

  void CompNovoIdentification::reducePermuts_(set<String> & permuts, const PeakSpectrum & CID_spec, const PeakSpectrum & ETD_spec, double prefix, double suffix)
  {
    if (permuts.size() < max_subscore_number_)
    {
      return;
    }

    vector<Permut> score_permuts;

    Size i(0);
    score_permuts.resize(permuts.size(), Permut(permuts.begin(), 0.0));
    for (set<String>::const_iterator it = permuts.begin(); it != permuts.end(); ++it, ++i)
    {
#ifdef REDUCE_PERMUTS_DEBUG
      if (i % 1000 == 0)
      {
        cerr << (double)i / permuts.size() * 100 << "%" << endl;
      }
#endif

      PeakSpectrum ETD_sim_spec, CID_sim_spec;
      getETDSpectrum_(ETD_sim_spec, *it, 1, prefix, suffix);
      getCIDSpectrumLight_(CID_sim_spec, *it, prefix, suffix);
      //getCIDSpectrum_(CID_sim_spec, *it, 1, prefix, suffix);

      double cid_score = zhang_(CID_sim_spec, CID_spec);
      //double cid_score = compareSpectra_(CID_spec, CID_sim_spec);
      /*if (isnan(cid_score))
          {
              cid_score = 0;
          }*/

      double etd_score = zhang_(ETD_sim_spec, ETD_spec);
      /*if (isnan(etd_score))
      {
          etd_score = 0;
      }*/

      //double etd_score = compareSpectra_(ETD_spec, ETD_sim_spec);
      double score = cid_score + etd_score;

      score /= it->size();

      if (boost::math::isnan(score))
      {
        score = 0;
      }

#ifdef REDUCE_PERMUTS_DEBUG
      cerr << "Subscoring: " << *it << " " << cid_score << " " << etd_score << " " << score << " (CID=";
/*      for (PeakSpectrum::ConstIterator pit = CID_sim_spec.begin(); pit != CID_sim_spec.end(); ++pit)
        {
        cerr << pit->getPosition()[0] << "|" << pit->getIntensity() << "; ";
        }*/
      cerr << endl;
#endif

      Permut new_permut(it, score);
      score_permuts[i].setScore(score);
      score_permuts[i].setPermut(it);

      //cerr << "permut=" << *it << ", score=" << score << endl;
    }

    //cerr << "Size of score_permuts: " << score_permuts.size() << endl;

    sort(score_permuts.begin(), score_permuts.end(), Internal::PermutScoreComparator);

    set<String> new_permuts;
    Size count(0);
    for (vector<Permut>::const_iterator it = score_permuts.begin(); it != score_permuts.end() && count < max_subscore_number_; ++it, ++count)
    {
      new_permuts.insert(*it->getPermut());
#ifdef REDUCE_PERMUTS_DEBUG
      cerr << "Subscore winner: " << *it->getPermut() << " " << it->getScore() << endl;
#endif
    }

    permuts = new_permuts;
    return;
  }

// divide and conquer algorithm of the sequencing
  void CompNovoIdentification::getDecompositionsDAC_(set<String> & sequences, Size left, Size right, double peptide_weight, const PeakSpectrum & CID_spec, const PeakSpectrum & ETD_spec, Map<double, CompNovoIonScoring::IonScore> & ion_scores)
  {
    static double oxonium_mass = EmpiricalFormula("H2O+").getMonoWeight();
    double offset_suffix(CID_spec[left].getPosition()[0] - oxonium_mass);
    double offset_prefix(peptide_weight - CID_spec[right].getPosition()[0]);

#ifdef DAC_DEBUG
    static int depth_(0);
    ++depth_;
    String tabs_(depth_, '\t');
    cerr << tabs_ << "void getDecompositionsDAC(sequences[" << sequences.size() << "], " << left << ", " << right << ") ";
    cerr << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " diff=";
#endif

    double diff = CID_spec[right].getPosition()[0] - CID_spec[left].getPosition()[0];

#ifdef DAC_DEBUG
    cerr << diff << endl;
    cerr << "offset_prefix=" << offset_prefix << ", offset_suffix=" << offset_suffix << endl;
#endif

    if (subspec_to_sequences_.has(left) && subspec_to_sequences_[left].has(right))
    {
      sequences = subspec_to_sequences_[left][right];

#ifdef DAC_DEBUG
      depth_--;
      cerr << tabs_ << "from cache DAC: " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << sequences.size() << " " << left << " " << right << endl;
#endif
      return;
    }

    // no further solutions possible?
    if (diff < min_aa_weight_)
    {
#ifdef DAC_DEBUG
      depth_--;
#endif
      return;
    }

    // no further division needed?
    if (diff <= max_decomp_weight_)
    {
      vector<MassDecomposition> decomps;
      getDecompositions_(decomps, diff);
      //filterDecomps_(decomps);

#ifdef DAC_DEBUG
      cerr << tabs_ << "Found " << decomps.size() << " decomps" << endl;
      cerr << tabs_ << "Permuting...";
#endif

      //static Map<String, set<String> > permute_cache;
      for (vector<MassDecomposition>::const_iterator it = decomps.begin(); it != decomps.end(); ++it)
      {
#ifdef DAC_DEBUG
        cerr << it->toString() << endl;
#endif

        String exp_string = it->toExpandedString();
        if (!permute_cache_.has(exp_string))
        {
          permute_("", exp_string, sequences);
          permute_cache_[exp_string] = sequences;
        }
        else
        {
          sequences = permute_cache_[exp_string];
        }
      }

#ifdef DAC_DEBUG
      cerr << tabs_ << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << peptide_weight << endl;
      if (sequences.size() > max_subscore_number_)
      {
        cerr << tabs_ << "Reducing #sequences from " << sequences.size() << " to " << max_subscore_number_ << "(prefix=" << offset_prefix  << ", suffix=" << offset_suffix << ")...";
      }
#endif

      // C-terminus
      if (offset_suffix <= fragment_mass_tolerance_)
      {
        filterPermuts_(sequences);
      }

      // reduce the sequences
      reducePermuts_(sequences, CID_spec, ETD_spec, offset_prefix, offset_suffix);
#ifdef DAC_DEBUG
      cerr << "Writing to cache " << left << " " << right << endl;
#endif
      subspec_to_sequences_[left][right] = sequences;

#ifdef DAC_DEBUG
      cerr << "ended" << endl;
      cerr << tabs_ << "DAC: " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << sequences.size() << endl;
      depth_--;
#endif

      return;
    }

    // select suitable pivot peaks
    vector<Size> pivots;

    if (offset_suffix < fragment_mass_tolerance_ && offset_prefix < fragment_mass_tolerance_)
    {
      selectPivotIons_(pivots, left, right, ion_scores, CID_spec, peptide_weight, true);
    }
    else
    {
      selectPivotIons_(pivots, left, right, ion_scores, CID_spec, peptide_weight, false);
    }

    // run divide step
#ifdef DAC_DEBUG
    cerr << tabs_ << "Selected " << pivots.size() << " pivot ions: ";
    for (vector<Size>::const_iterator it = pivots.begin(); it != pivots.end(); ++it)
    {
      cerr << *it << "(" << CID_spec[*it].getPosition()[0] << ") ";
    }
    cerr << endl;
#endif

    for (vector<Size>::const_iterator it = pivots.begin(); it != pivots.end(); ++it)
    {
      set<String> seq1, seq2, new_sequences;

      // the smaller the 'gap' the greater the chance of not finding anything
      // so we we compute the smaller gap first
      double diff1(CID_spec[*it].getPosition()[0] - CID_spec[left].getPosition()[0]);
      double diff2(CID_spec[right].getPosition()[0] - CID_spec[*it].getPosition()[0]);

      if (diff1 < diff2)
      {
        getDecompositionsDAC_(seq1, left, *it, peptide_weight, CID_spec, ETD_spec, ion_scores);
        if (seq1.empty())
        {
#ifdef DAC_DEBUG
          cerr << tabs_ << "first call produced 0 candidates (" << diff1 << ")" << endl;
#endif
          continue;
        }

        getDecompositionsDAC_(seq2, *it, right, peptide_weight, CID_spec, ETD_spec, ion_scores);
      }
      else
      {
        getDecompositionsDAC_(seq2, *it, right, peptide_weight, CID_spec, ETD_spec, ion_scores);
        if (seq2.empty())
        {
#ifdef DAC_DEBUG
          cerr << tabs_ << "second call produced 0 candidates (" << diff2 << ")" << endl;
#endif
          continue;
        }

        getDecompositionsDAC_(seq1, left, *it, peptide_weight, CID_spec, ETD_spec, ion_scores);
      }

#ifdef DAC_DEBUG
      cerr << tabs_ << "Found " << seq1.size() << " solutions (1) " << diff1 << endl;
      cerr << tabs_ << "Found " << seq2.size() << " solutions (2) " << diff2 << endl;
      cerr << tabs_ << "inserting " << seq1.size() * seq2.size()  << " sequences" << endl;
#endif

      // C-terminus
      if (offset_suffix <= fragment_mass_tolerance_)
      {
        filterPermuts_(seq1);
      }

      // test if we found enough sequence candidates
      if (seq1.empty() || seq2.empty())
      {
        continue;
      }

      for (set<String>::const_iterator it1 = seq1.begin(); it1 != seq1.end(); ++it1)
      {
        for (set<String>::const_iterator it2 = seq2.begin(); it2 != seq2.end(); ++it2)
        {
          new_sequences.insert(*it2 + *it1);
        }
      }

      if (seq1.size() * seq2.size() > max_subscore_number_)
      {
#ifdef DAC_DEBUG
        cerr << tabs_ << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << peptide_weight << endl;
        cerr << tabs_ << "Reducing #sequences from " << new_sequences.size() << " to " << max_subscore_number_ << "(prefix=" << offset_prefix  << ", suffix=" << offset_suffix << ")...";
#endif
        if (offset_prefix > fragment_mass_tolerance_ || offset_suffix > fragment_mass_tolerance_)
        {
          reducePermuts_(new_sequences, CID_spec, ETD_spec, offset_prefix, offset_suffix);
        }

#ifdef DAC_DEBUG
        for (set<String>::const_iterator it1 = new_sequences.begin(); it1 != new_sequences.end(); ++it1)
        {
          cerr << tabs_ << *it1 << endl;
        }
        cerr << endl;
#endif
      }

      for (set<String>::const_iterator sit = new_sequences.begin(); sit != new_sequences.end(); ++sit)
      {
        sequences.insert(*sit);
      }
    }
#ifdef DAC_DEBUG
    cerr << tabs_ << "Found sequences for " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << endl;
    for (set<String>::const_iterator sit = sequences.begin(); sit != sequences.end(); ++sit)
    {
      cerr << tabs_ << *sit << endl;
    }
#endif

    // reduce the permuts once again to reduce complexity
    if (offset_prefix > fragment_mass_tolerance_ || offset_suffix > fragment_mass_tolerance_)
    {
      reducePermuts_(sequences, CID_spec, ETD_spec, offset_prefix, offset_suffix);
    }

#ifdef DAC_DEBUG
    cerr << "Writing to cache " << left << " " << right << endl;
#endif

    subspec_to_sequences_[left][right] = sequences;

#ifdef DAC_DEBUG
    depth_--;
    cerr << tabs_ << "DAC: " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << sequences.size() << endl;
#endif
    return;

  }

  double CompNovoIdentification::estimatePrecursorWeight_(const PeakSpectrum & ETD_spec, Size & charge)
  {
    CompNovoIonScoring ion_scoring;
    double precursor_mass_tolerance((double)param_.getValue("precursor_mass_tolerance"));
    double peptide_weight(0);
    // for each possible charge state just get all possible precursor peaks
    double precursor_mz(ETD_spec.getPrecursors().begin()->getMZ());
    Map<Size, Map<Size, vector<Peak1D> > > peaks;
    Map<Size, Map<Size, vector<double> > > correlations;
    for (PeakSpectrum::ConstIterator it = ETD_spec.begin(); it != ETD_spec.end(); ++it)
    {
      for (Size prec_z = 1; prec_z <= 3; ++prec_z)
      {
        for (Size peptide_z = 2; peptide_z <= 3; ++peptide_z)
        {
          if (prec_z > peptide_z)
          {
            continue;
          }
          double pre_mz = ((double)(precursor_mz * peptide_z) - (double)(peptide_z - prec_z) * Constants::PROTON_MASS_U) / (double)prec_z;
          if (fabs(it->getMZ() * (double)prec_z - pre_mz * (double)prec_z) < precursor_mass_tolerance /* / (double) prec_z*/)
          {
            peaks[peptide_z][prec_z].push_back(*it);
            correlations[peptide_z][prec_z].push_back(ion_scoring.scoreIsotopes(ETD_spec, it, prec_z) /* *  it->getIntensity()*/);
            //cerr << "Assumed charge=" << peptide_z << ", precursor peak z=" << prec_z << ", experimental m/z=" << precursor_mz << ", diff=" << fabs(it->getMZ() - pre_mz) << " m/z=" << it->getMZ() << ", int=" << it->getIntensity() << " " << correlations[peptide_z][prec_z].back() <<  endl;
          }
        }
      }
    }

    // for each possible peptide charge state
    Map<Size, double> correlation_sums;
    Map<Size, Map<Size, pair<double, double> > > best_corr_ints;
    for (Map<Size, Map<Size, vector<double> > >::ConstIterator it1 = correlations.begin(); it1 != correlations.end(); ++it1)
    {
      double correlation_sum(0);
      // search for the best correlation
      for (Map<Size, vector<double> >::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        double best_correlation(0);
        Size best_pos(0);
        Size i = 0;
        for (vector<double>::const_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3, ++i)
        {
          if (best_correlation == 0 || *it3 > (best_correlation * 1.25))           // must be really better!
          {
            best_correlation = *it3;
            best_pos = i;
          }
        }
        best_corr_ints[it1->first][it2->first] = make_pair(best_correlation, peaks[it1->first][it2->first][best_pos].getMZ());
        correlation_sum += best_correlation;
      }

      correlation_sums[it1->first] = correlation_sum;
    }

    double best_correlation = 0;
    Size best_charge = 0;
    for (Map<Size, double>::ConstIterator it = correlation_sums.begin(); it != correlation_sums.end(); ++it)
    {
      //cerr << "Correlations z=" << it->first << ", corr=" << it->second << endl;
      for (Map<Size, pair<double, double> >::ConstIterator mit = best_corr_ints[it->first].begin(); mit != best_corr_ints[it->first].end(); ++mit)
      {
        //cerr << "CorrelationIntensity: z=" << mit->first << ", corr=" << mit->second.first << ", m/z=" << mit->second.second << " [M+H]=" << (mit->second.second * (double)mit->first) - ((double)mit->first - 1) * Constants::NEUTRON_MASS_U  << endl;
      }
      if (best_correlation < it->second)
      {
        best_correlation = it->second;
        best_charge = it->first;
      }
    }

    //cerr << "Best correlation=" << best_correlation << " best_charge=" << best_charge << endl;
    charge = best_charge;

    // check whether charge one is available
    if (best_corr_ints[best_charge].has(1))
    {
      peptide_weight = best_corr_ints[best_charge][1].second;
    }
    else
    {
      // find other best charge state
      best_correlation = 0;
      double best_corr_mz = 0;
      Size best_corr_z = 0;
      for (Map<Size, pair<double, double> >::ConstIterator it = best_corr_ints[best_charge].begin(); it != best_corr_ints[best_charge].end(); ++it)
      {
        if (it->second.first > best_correlation)
        {
          best_correlation = it->second.first;
          best_corr_mz = it->second.second;
          best_corr_z = it->first;
        }
      }
      peptide_weight = best_corr_mz * (double)best_corr_z - (double)(best_corr_z - 1) * Constants::PROTON_MASS_U;
      //cerr << "BestCorr: " << best_correlation << " " << best_corr_mz << " " << best_corr_z << " " << peptide_weight << endl;
    }

    return peptide_weight;
  }

/*
    double CompNovoIdentification::estimatePrecursorWeight_(const PeakSpectrum& ETD_spec, Size& charge)
    {
        CompNovoIonScoring ion_scoring;
        double precursor_mass_tolerance((double)param_.getValue("precursor_mass_tolerance"));
        double precursor_weight(0.0);

        // first we assume the charge is 3;
        double pre_weight_3_z2 = (ETD_spec.getPrecursors().begin()->getMZ() * 3.0 - 1.0 * Constants::PROTON_MASS_U) / 2.0;

        // now get the charge 2 peaks
        vector<double> precursor_ints_3_z2, iso_scores_3_z2;
        vector<Peak1D> precursor_peaks_3_z2;
      for (PeakSpectrum::ConstIterator it = ETD_spec.begin(); it != ETD_spec.end(); ++it)
    {
      if (fabs(it->getPosition()[0] - pre_weight_3_z2) < precursor_mass_tolerance)
      {
          double iso_score = ion_scoring.scoreIsotopes(ETD_spec, it, 2);
          precursor_ints_3_z2.push_back(it->getIntensity());
          precursor_peaks_3_z2.push_back(*it);
                    iso_scores_3_z2.push_back(iso_score);
#ifdef ESTIMATE_PRECURSOR_DEBUG
                    cerr << "Pre: " << it->getPosition()[0] << " 2 " << iso_score << " int=" << it->getIntensity() << endl;
#endif
      }
    }

        double pre_weight_2_z1 = ETD_spec.getPrecursors().begin()->getMZ() * 2.0 - 1.0 * Constants::PROTON_MASS_U;
        vector<double> precursor_ints_2_z1, iso_scores_2_z1;
        vector<Peak1D> precursor_peaks_2_z1;
        for (PeakSpectrum::ConstIterator it = ETD_spec.begin(); it != ETD_spec.end(); ++it)
        {
            if (fabs(it->getPosition()[0] - pre_weight_2_z1) < precursor_mass_tolerance)
            {
                double iso_score = ion_scoring.scoreIsotopes(ETD_spec, it, 1);
                precursor_ints_2_z1.push_back(it->getIntensity());
                precursor_peaks_2_z1.push_back(*it);
                iso_scores_2_z1.push_back(iso_score);
#ifdef ESTIMATE_PRECURSOR_DEBUG
                cerr << "Pre: " << it->getPosition()[0] << " 1 " << iso_score << " int=" << it->getIntensity() << endl;
#endif
            }
        }

        double pre_weight_2_z2 = ETD_spec.getPrecursors().begin()->getMZ();
        vector<double> precursor_ints_2_z2, iso_scores_2_z2;
        vector<Peak1D> precursor_peaks_2_z2;
        for (PeakSpectrum::ConstIterator it = ETD_spec.begin(); it != ETD_spec.end(); ++it)
        {
            if (fabs(it->getPosition()[0] - pre_weight_2_z2) < precursor_mass_tolerance)
            {
                double iso_score = ion_scoring.scoreIsotopes(ETD_spec, it, 2);
                precursor_ints_2_z2.push_back(it->getIntensity());
                precursor_peaks_2_z2.push_back(*it);
                iso_scores_2_z2.push_back(iso_score);
#ifdef ESTIMATE_PRECURSOR_DEBUG
                cerr << "Pre: " << it->getPosition()[0] << " 2 " << iso_score << " int=" << it->getIntensity() << endl;
#endif
            }
        }

        // now decide which charge variant is more likely
        double max_element_z2(0), max_element_z3(0);
        if (iso_scores_3_z2.size() > 0)
        {
            max_element_z3 = *max_element(iso_scores_3_z2.begin(), iso_scores_3_z2.end());

            if (max_element_z3 < 0)
            {
                // isotope scoring was not successful, only decide on the intensities, however scale to prefer ions which clearly have good isopattern
                max_element_z3 = *max_element(precursor_ints_3_z2.begin(), precursor_ints_3_z2.end()) / 100;
            }
            else
            {
                max_element_z3 *= *max_element(precursor_ints_3_z2.begin(), precursor_ints_3_z2.end());
            }
        }

        if (iso_scores_2_z1.size() > 0)
        {
            max_element_z2 = *max_element(iso_scores_2_z1.begin(), iso_scores_2_z1.end());
            if (max_element_z2 < 0)
            {
                // isotope scoring was not successful, only decide on the intensities, however scale to prefer ions which clearly have good isopattern
                max_element_z2 = *max_element(precursor_ints_2_z1.begin(), precursor_ints_2_z1.end()) / 100;
            }
            else
            {
                max_element_z2 *= *max_element(precursor_ints_2_z1.begin(), precursor_ints_2_z1.end());
            }
        }

        double max_element_z2_2(0);
        if (iso_scores_2_z2.size() > 0)
        {
            max_element_z2_2 = *max_element(iso_scores_2_z2.begin(), iso_scores_2_z2.end());
            if (max_element_z2_2 < 0)
            {
                max_element_z2_2 = *max_element(precursor_ints_2_z2.begin(), precursor_ints_2_z2.end());
            }
            else
            {
                max_element_z2_2 *= *max_element(precursor_ints_2_z2.begin(), precursor_ints_2_z2.end());
            }
        }

        cerr << "max_elements: " << max_element_z3 << " " << max_element_z2 << endl;

        if (max_element_z3 > max_element_z2 && max_element_z3 > max_element_z2_2)
        {
            charge = 3;
        }
        else
        {
            charge = 2;
        }


        // now get the exact precursor weight of the charge variant

        if (charge == 3 && precursor_peaks_3_z2.size() > 0)
        {
            precursor_weight = precursor_peaks_3_z2.begin()->getMZ() * 2.0 - 1.0 * Constants::PROTON_MASS_U - 1.0 * Constants::NEUTRON_MASS_U;
        }
        else
        {
            if (precursor_ints_2_z1.size() > 1)
            {
                precursor_weight = precursor_peaks_2_z1.begin()->getMZ();
                double pre_int_first = *precursor_ints_2_z1.begin();
                double pre_int_second = *(++precursor_ints_2_z1.begin());

                if (pre_int_second / pre_int_first > 10)
                {
                    precursor_weight = (++precursor_peaks_2_z1.begin())->getMZ();
                }

                if (iso_scores_2_z1.size() > 1)
                {
                    cerr << "Corr. of: " << precursor_peaks_2_z1[0].getMZ() << " " << iso_scores_2_z1[0] << endl;
                    cerr << "Corr. of: " << precursor_peaks_2_z1[1].getMZ() << " " << iso_scores_2_z1[1] << endl;
                    double iso_diff = iso_scores_2_z1[1] - iso_scores_2_z1[0];
                    cerr << "Corr. diff: " << iso_diff << endl;

                    if ((iso_scores_2_z1[0] < 0 && iso_scores_2_z1[1] > 0) ||
                            (iso_scores_2_z1[0] > 0 && iso_scores_2_z1[1] > 0 &&
                             iso_diff > 0.2
                            ))
                    {
                        precursor_weight = (++precursor_peaks_2_z1.begin())->getMZ();
                    }

                    if (iso_scores_2_z1.size() > 2 &&
                            iso_scores_2_z1[2] > 0)
                    {
                        iso_diff = iso_scores_2_z1[2] - iso_scores_2_z1[1];
                        cerr << "Corr. of: " << precursor_peaks_2_z1[2].getMZ() << " " << iso_scores_2_z1[2] << endl;
                        cerr << "Corr. diff: " << iso_diff << endl;
                        if (iso_diff > 0.5)
                        {
                            precursor_weight = (++(++precursor_peaks_2_z1.begin()))->getMZ();
                        }
                    }
                }
            }
        }

        return precursor_weight;
    }
    */
}
