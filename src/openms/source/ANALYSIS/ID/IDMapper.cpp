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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace std;

namespace OpenMS
{

  IDMapper::IDMapper() :
    DefaultParamHandler("IDMapper"),
    rt_tolerance_(5.0),
    mz_tolerance_(20),
    measure_(MEASURE_PPM),
    ignore_charge_(false)
  {
    defaults_.setValue("rt_tolerance", rt_tolerance_, "RT tolerance (in seconds) for the matching");
    defaults_.setMinFloat("rt_tolerance", 0);
    defaults_.setValue("mz_tolerance", mz_tolerance_, "m/z tolerance (in ppm or Da) for the matching");
    defaults_.setMinFloat("mz_tolerance", 0);
    defaults_.setValue("mz_measure", "ppm", "unit of 'mz_tolerance' (ppm or Da)");
    defaults_.setValidStrings("mz_measure", ListUtils::create<String>("ppm,Da"));
    defaults_.setValue("mz_reference", "precursor", "source of m/z values for peptide identifications");
    defaults_.setValidStrings("mz_reference", ListUtils::create<String>("precursor,peptide"));

    defaults_.setValue("ignore_charge", "false", "For feature/consensus maps: Assign an ID independently of whether its charge state matches that of the (consensus) feature.");
    defaults_.setValidStrings("ignore_charge", ListUtils::create<String>("true,false"));

    defaultsToParam_();
  }

  IDMapper::IDMapper(const IDMapper & cp) :
    DefaultParamHandler(cp),
    rt_tolerance_(cp.rt_tolerance_),
    mz_tolerance_(cp.mz_tolerance_),
    measure_(cp.measure_),
    ignore_charge_(cp.ignore_charge_)
  {
    updateMembers_();
  }

  IDMapper & IDMapper::operator=(const IDMapper & rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    rt_tolerance_ = rhs.rt_tolerance_;
    mz_tolerance_ = rhs.mz_tolerance_;
    measure_ = rhs.measure_;
    ignore_charge_ = rhs.ignore_charge_;
    updateMembers_();

    return *this;
  }

  void IDMapper::updateMembers_()
  {
    rt_tolerance_ = param_.getValue("rt_tolerance");
    mz_tolerance_ = param_.getValue("mz_tolerance");
    measure_ = param_.getValue("mz_measure") == "ppm" ? MEASURE_PPM : MEASURE_DA;
    ignore_charge_ = param_.getValue("ignore_charge") == "true";
  }

  void IDMapper::annotate(ConsensusMap & map, const std::vector<PeptideIdentification> & ids, const std::vector<ProteinIdentification> & protein_ids, bool measure_from_subelements)
  {
    // validate "RT" and "MZ" metavalues exist
    checkHits_(ids);

    //append protein identifications to Map
    map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

    //keep track of assigned/unassigned peptide identifications
    std::map<Size, Size> assigned;

    // store which peptides fit which feature (and avoid double entries)
    // consensusMap -> {peptide_index}
    std::vector<std::set<size_t> > mapping(map.size());

    DoubleList mz_values;
    double rt_pep;
    IntList charges;

    //iterate over the peptide IDs
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (ids[i].getHits().empty())
        continue;

      getIDDetails_(ids[i], rt_pep, mz_values, charges);

      //iterate over the features
      for (Size cm_index = 0; cm_index < map.size(); ++cm_index)
      {
        // if set to TRUE, we leave the i_mz-loop as we added the whole ID with all hits
        bool was_added = false;       // was current pep-m/z matched?!

        // iterate over m/z values of pepIds
        for (Size i_mz = 0; i_mz < mz_values.size(); ++i_mz)
        {
          double mz_pep = mz_values[i_mz];

          // charge states to use for checking:
          IntList current_charges;
          if (!ignore_charge_)
          {
            // if "mz_ref." is "precursor", we have only one m/z value to check,
            // but still one charge state per peptide hit that could match:
            if (mz_values.size() == 1)
            {
              current_charges = charges;
            }
            else
            {
              current_charges.push_back(charges[i_mz]);
            }
            current_charges.push_back(0);             // "not specified" always matches
          }

          //check if we compare distance from centroid or subelements
          if (!measure_from_subelements)
          {
            if (isMatch_(rt_pep - map[cm_index].getRT(), mz_pep, map[cm_index].getMZ()) && (ignore_charge_ || ListUtils::contains(current_charges, map[cm_index].getCharge())))
            {
              was_added = true;
              map[cm_index].getPeptideIdentifications().push_back(ids[i]);
              ++assigned[i];
            }
          }
          else
          {
            for (ConsensusFeature::HandleSetType::const_iterator it_handle = map[cm_index].getFeatures().begin();
                 it_handle != map[cm_index].getFeatures().end();
                 ++it_handle)
            {
              if (isMatch_(rt_pep - it_handle->getRT(), mz_pep, it_handle->getMZ())  && (ignore_charge_ || ListUtils::contains(current_charges, it_handle->getCharge())))
              {
                was_added = true;
                if (mapping[cm_index].count(i) == 0)
                {
                  map[cm_index].getPeptideIdentifications().push_back(ids[i]);
                  ++assigned[i];
                  mapping[cm_index].insert(i);
                }
                break;                 // we added this peptide already.. no need to check other handles
              }
            }
            // continue to here
          }

          if (was_added)
            break;

        }         // m/z values to check

        // break to here

      }       // features
    }     // Identifications


    Size matches_none(0);
    Size matches_single(0);
    Size matches_multi(0);

    //append unassigned peptide identifications
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (assigned[i] == 0)
      {
        map.getUnassignedPeptideIdentifications().push_back(ids[i]);
        ++matches_none;
      }
      else if (assigned[i] == 1)
      {
        ++matches_single;
      }
      else if (assigned[i] > 1)
      {
        ++matches_multi;
      }
    }

    //some statistics output
    LOG_INFO << "Unassigned peptides: " << matches_none << "\n"
             << "Peptides assigned to exactly one feature: "
             << matches_single << "\n"
             << "Peptides assigned to multiple features: "
             << matches_multi << std::endl;

  }

  double IDMapper::getAbsoluteMZTolerance_(const double mz) const
  {
    if (measure_ == MEASURE_PPM)
    {
      return mz * mz_tolerance_ / 1e6;
    }
    else if (measure_ == MEASURE_DA)
    {
      return mz_tolerance_;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDMapper::getAbsoluteTolerance_(): illegal internal state of measure_!", String(measure_));
  }

  bool IDMapper::isMatch_(const double rt_distance, const double mz_theoretical, const double mz_observed) const
  {
    if (measure_ == MEASURE_PPM)
    {
      return (fabs(rt_distance) <= rt_tolerance_) && (fabs((mz_theoretical * 1e6 - mz_observed * 1e6) / mz_theoretical) <= mz_tolerance_);
    }
    else if (measure_ == MEASURE_DA)
    {
      return (fabs(rt_distance) <= rt_tolerance_) && (fabs(mz_theoretical - mz_observed) <= mz_tolerance_);
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDMapper::getAbsoluteTolerance_(): illegal internal state of measure_!", String(measure_));
  }

  void IDMapper::checkHits_(const std::vector<PeptideIdentification> & ids) const
  {
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (!ids[i].metaValueExists("RT"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDMapper: meta data value 'RT' missing for peptide identification!");
      }
      if (!ids[i].metaValueExists("MZ"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDMapper: meta data value 'MZ' missing for peptide identification!");
      }
    }
  }

  void IDMapper::getIDDetails_(const PeptideIdentification & id, double & rt_pep, DoubleList & mz_values, IntList & charges, bool use_avg_mass) const
  {
    mz_values.clear();
    charges.clear();

    rt_pep = id.getMetaValue("RT");

    // collect m/z values of pepId
    if (param_.getValue("mz_reference") == "precursor") // use precursor m/z of pepId
    {
      mz_values.push_back(id.getMetaValue("MZ"));
    }

    for (vector<PeptideHit>::const_iterator hit_it = id.getHits().begin();
         hit_it != id.getHits().end(); ++hit_it)
    {
      Int charge = hit_it->getCharge();
      charges.push_back(charge);

      if (param_.getValue("mz_reference") == "peptide") // use mass of each pepHit (assuming H+ adducts)
      {
        double mass = use_avg_mass ?
                          hit_it->getSequence().getAverageWeight(Residue::Full, charge) :
                          hit_it->getSequence().getMonoWeight(Residue::Full, charge);

        mz_values.push_back( mass / (double) charge);
      }
    }
  }

  void IDMapper::increaseBoundingBox_(DBoundingBox<2> & box)
  {
    DPosition<2> sub_min(rt_tolerance_,
                         getAbsoluteMZTolerance_(box.minPosition().getY())),
    add_max(rt_tolerance_, getAbsoluteMZTolerance_(box.maxPosition().getY()));

    box.setMin(box.minPosition() - sub_min);
    box.setMax(box.maxPosition() + add_max);
  }

  bool IDMapper::checkMassType_(const vector<DataProcessing> & processing) const
  {
    bool use_avg_mass = false;
    String before;
    for (vector<DataProcessing>::const_iterator proc_it = processing.begin();
         proc_it != processing.end(); ++proc_it)
    {
      if (proc_it->getSoftware().getName() == "FeatureFinder")
      {
        String reported_mz = proc_it->
                             getMetaValue("parameter: algorithm:feature:reported_mz");
        if (reported_mz.empty())
          continue;                                // parameter info not available
        if (!before.empty() && (reported_mz != before))
        {
          LOG_WARN << "The m/z values reported for features in the input seem to be of different types (e.g. monoisotopic/average). They will all be compared against monoisotopic peptide masses, but the mapping results may not be meaningful in the end." << endl;
          return false;
        }
        if (reported_mz == "average")
        {
          use_avg_mass = true;
        }
        else if (reported_mz == "maximum")
        {
          LOG_WARN << "For features, m/z values from the highest mass traces are reported. This type of m/z value is not available for peptides, so the comparison has to be done using average peptide masses." << endl;
          use_avg_mass = true;
        }
        before = reported_mz;
      }
    }
    return use_avg_mass;
  }

} // namespace OpenMS
