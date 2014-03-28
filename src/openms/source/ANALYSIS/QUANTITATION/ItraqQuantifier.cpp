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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>

// default isotope correction (via A^-1)
// NNLS isotope correction
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/Utils/MatrixUtils.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include <limits>



//#define ITRAQ_DEBUG 1

namespace OpenMS
{
  ItraqQuantifier::ItraqQuantifier() :
    DefaultParamHandler("ItraqQuantifier"),
    itraq_type_(FOURPLEX),
    isotope_corrections_()
  {
    initIsotopeCorrections_();
    setDefaultParams_();
  }

  ItraqQuantifier::ItraqQuantifier(Int itraq_type) :
    DefaultParamHandler("ItraqQuantifier"),
    itraq_type_(itraq_type),
    isotope_corrections_()
  {
    initIsotopeCorrections_();
    setDefaultParams_();
  }

  ItraqQuantifier::ItraqQuantifier(Int itraq_type, const Param& param) :
    DefaultParamHandler("ItraqQuantifier"),
    itraq_type_(itraq_type),
    isotope_corrections_()
  {
    initIsotopeCorrections_();
    setDefaultParams_();
    setParameters(param);
    updateMembers_();
  }

  ItraqQuantifier::ItraqQuantifier(const ItraqQuantifier& cp) :
    DefaultParamHandler(cp),
    ItraqConstants(cp),
    itraq_type_(cp.itraq_type_),
    channel_map_(cp.channel_map_),
    isotope_corrections_(cp.isotope_corrections_)
  {
  }

  ItraqQuantifier& ItraqQuantifier::operator=(const ItraqQuantifier& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    ItraqConstants::operator=(rhs);
    itraq_type_ = rhs.itraq_type_;
    channel_map_ = rhs.channel_map_;
    isotope_corrections_ = rhs.isotope_corrections_;

    return *this;
  }
  
  bool ItraqQuantifier::isIdentityCorrectionMatrix_(const Matrix<double>& channel_frequency) const
  {
    // check if we have an identity matrix
    bool isIdentity = true;
    for (Size i = 0; i < channel_frequency.cols(); ++i)
    {
      if (channel_frequency.getValue(i,i) != 1.0)
      {
        isIdentity = false;
        break;
      }
    }
    return isIdentity;
  }

  void ItraqQuantifier::run(const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out)
  {
    // new stats
    stats_ = ItraqQuantifierStats();
    stats_.channel_count = CHANNEL_COUNT[itraq_type_];
    if (consensus_map_in.empty())
    {
      LOG_WARN << "Warning: Empty iTRAQ container. No quantitative information available!" << std::endl;
      return;
    }

    reconstructChannelInfo_(consensus_map_in);
    consensus_map_out = consensus_map_in;

    // first do isotope correction
    if (String(param_.getValue("isotope_correction")) == "true")
    {
      // translate isotope_corrections_ to a channel_frequency matrix
      Matrix<double> channel_frequency = ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_);

      // if it is an identity matrix, performing isotope correction makes no sense
      if (isIdentityCorrectionMatrix_(channel_frequency))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier: The given isotope correction matrix is an identity matrix leading to no correction. Please provide a valid isotope_correction matrix as it was provided with the iTRAQ/TMT kit!");
      }
      
#ifdef ITRAQ_DEBUG
      std::cout << "channel_frequency matrix: \n" << channel_frequency << "\n" << std::endl;
#endif

      // ISOTOPE CORRECTION: this solves the system naively via matrix inversion
      EigenMatrixXdPtr m ( convertOpenMSMatrix2EigenMatrixXd( channel_frequency ) );
      Eigen::FullPivLU<Eigen::MatrixXd> ludecomp (*m);
      Eigen::VectorXd b;
      b.resize(CHANNEL_COUNT[itraq_type_]);
      b.setZero();
      std::vector<double> x (CHANNEL_COUNT[itraq_type_], 0);


      if(!ludecomp.isInvertible())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; the Matrix is not invertible!");
      }

      LOG_INFO << "SOLVING isotope correction via NNLS\n";

      Matrix<double> m_b(CHANNEL_COUNT[itraq_type_], 1);
      Matrix<double> m_x(CHANNEL_COUNT[itraq_type_], 1);

      // correct all consensus elements
      for (size_t i = 0; i < consensus_map_out.size(); ++i)
      {
#ifdef ITRAQ_DEBUG
        std::cout << "\nMAP element  #### " << i << " #### \n" << std::endl;
#endif

        consensus_map_out[i].clear(); // delete only the consensus handles
        // fill b vector
        for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_in[i].getFeatures().begin();
             it_elements != consensus_map_in[i].getFeatures().end();
             ++it_elements)
        {


          //find channel_id of current element
          Int index = Int(consensus_map_in.getFileDescriptions()[it_elements->getMapIndex()].getMetaValue("channel_id"));
#ifdef ITRAQ_DEBUG
          std::cout << "  map_index " << it_elements->getMapIndex() << "-> id " << index << " with intensity " << it_elements->getIntensity() << "\n" << std::endl;
#endif

          // this is deprecated, but serves as quality measurement
          b(index) = it_elements->getIntensity();
          m_b(index, 0) = it_elements->getIntensity();
        }

        // solve
        Eigen::MatrixXd x = ludecomp.solve( b );
        // check if a solutioon exists
        if (! ((*m) * x).isApprox(b))
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; Cannot multiply!");
        }
        Int status = NonNegativeLeastSquaresSolver::solve(channel_frequency, m_b, m_x);
        if (status != NonNegativeLeastSquaresSolver::SOLVED)
        {
          throw Exception::FailedAPICall(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier: Failed to find least-squares fit!");
        }

        Size s_negative(0);
        Size s_different_count(0); // happens when naive solution is negative in other channels
        DoubleReal s_different_intensity(0);
        // ISOTOPE CORRECTION: compare solutions of Matrix inversion vs. NNLS
        for (Size index = 0; index < (Size)CHANNEL_COUNT[itraq_type_]; ++index)
        {
          if (x(index) < 0.0)
          {
            ++s_negative;
          }
          else if (std::fabs(m_x(index, 0) - x(index)) > 0.000001)
          {
            ++s_different_count;
            s_different_intensity += std::fabs(m_x(index, 0) - x(index));
          }
        }

        if (s_negative == 0 && s_different_count > 0) // solutions are inconsistent, despite being positive! This should not happen!
        {
          throw Exception::Postcondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Isotope correction values of alternative method differ!");
        }

        // update global stats
        stats_.iso_number_reporter_negative += s_negative;
        stats_.iso_number_reporter_different += s_different_count;
        stats_.iso_solution_different_intensity += s_different_intensity;

        // write back the values to the map
        Peak2D::IntensityType cf_intensity(0);
        for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_in[i].begin();
             it_elements != consensus_map_in[i].end();
             ++it_elements)
        {
          FeatureHandle handle = *it_elements;
          //find channel_id of current element
          Int index = Int(consensus_map_out.getFileDescriptions()[it_elements->getMapIndex()].getMetaValue("channel_id"));

          handle.setIntensity(Peak2D::IntensityType(m_x(index, 0)));

          consensus_map_out[i].insert(handle);

          cf_intensity += handle.getIntensity(); // sum up all channels for CF

#ifdef ITRAQ_DEBUG
          std::cout <<  it_elements->getIntensity() << " -> " << handle.getIntensity() << std::endl;
#endif
        }
        consensus_map_out[i].setIntensity(cf_intensity); // set overall intensity of CF (sum of all channels)

        if (s_negative > 0)
        {
          ++stats_.iso_number_ms2_negative;
          stats_.iso_total_intensity_negative += cf_intensity;
        }
      }

    } // ! isotope_correction
    else
    {
      LOG_WARN << "Warning: Due to deactivated isotope-correction labeling statistics will be based on raw intensities, which might give too optimistic results." << std::endl;
    }

    stats_.number_ms2_total = consensus_map_out.size();
    // ------------------------------
    // Labeling efficiency statistics
    // ------------------------------
    std::map<Size, Size> empty_channel;

    for (size_t i = 0; i < consensus_map_out.size(); ++i)
    {
      // is whole scan empty?!
      if (consensus_map_out[i].getIntensity() == 0)
        ++stats_.number_ms2_empty;

      // look at single reporters
      for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_out[i].begin();
           it_elements != consensus_map_out[i].end();
           ++it_elements)
      {
        if (it_elements->getIntensity() == 0)
        {
          Int ch_index = consensus_map_out.getFileDescriptions()[it_elements->getMapIndex()].getMetaValue("channel_name");
          ++empty_channel[ch_index];
        }
      }
    }
    LOG_INFO << "iTRAQ: skipped " << stats_.number_ms2_empty << " of " << consensus_map_out.size() << " selected scans due to lack of iTRAQ information:\n";
    consensus_map_out.setMetaValue("itraq:scans_noquant", stats_.number_ms2_empty);
    consensus_map_out.setMetaValue("itraq:scans_total", consensus_map_out.size());

    stats_.empty_channels = empty_channel;

    LOG_INFO << "iTRAQ: channels with signal\n";
    for (std::map<Size, Size>::const_iterator it_m = empty_channel.begin(); it_m != empty_channel.end(); ++it_m)
    {
      LOG_INFO << "      channel " << it_m->first << ": " << (consensus_map_out.size() - it_m->second) << " / " <<  consensus_map_out.size() << " (" << ((consensus_map_out.size() - it_m->second) * 100 / consensus_map_out.size()) << "%)\n";
      consensus_map_out.setMetaValue(String("itraq:quantifyable_ch") + it_m->first, (consensus_map_out.size() - it_m->second));
    }

    // ****************************
    // ** find reference channel **
    // ****************************
    Int reference_channel = Int(param_.getValue("channel_reference"));
    if (itraq_type_ == ItraqConstants::FOURPLEX && (reference_channel < 114 || reference_channel > 117))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier:Invalid entry in Param 'channel_reference'; Valid channels for 4plex are 114-117!");
    }
    else if (itraq_type_ == ItraqConstants::EIGHTPLEX && (reference_channel < 113 || reference_channel > 121))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier:Invalid entry in Param 'channel_reference'; Valid channels for 8plex are 113-121!");
    }
    else if (itraq_type_ == ItraqConstants::TMT_SIXPLEX && (reference_channel < 126 || reference_channel > 131))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier:Invalid entry in Param 'channel_reference'; Valid channels for TMT-6plex are 126-131!");
    }

#ifdef ITRAQ_DEBUG
    std::cout << "reference_channel is: " << reference_channel  << std::endl;
#endif

    // determine reference channel as vector index
    Map<Size, Size> map_to_vectorindex;
    Size ref_mapid = 0;
    Size index = 0;
    for (ConsensusMap::FileDescriptions::const_iterator file_it = consensus_map_out.getFileDescriptions().begin();
         file_it != consensus_map_out.getFileDescriptions().end();
         ++file_it)
    {
      if ((Int) file_it->second.getMetaValue("channel_name") == reference_channel)
      {
        ref_mapid = file_it->first;
#ifdef ITRAQ_DEBUG
        std::cout << "reference_map_id is: " << ref_mapid <<  std::endl;
#endif
      }
      map_to_vectorindex[file_it->first] = index;
      ++index;
    }

    // ** NORMALIZATION ** //

    // normalize median of channel-to-reference ratio to 1
    if (String(param_.getValue("do_normalization")) == "true")
    {
      if (channel_map_.has(reference_channel))
      {
        std::vector<std::vector<double> > peptide_ratios;
        // this is a control (the normalization factors should be about the same)
        std::vector<std::vector<double> > peptide_intensities;

        // build mapping of map_index to ratio_array_index
        peptide_ratios.resize(channel_map_.size());
        peptide_intensities.resize(channel_map_.size());

        //build up ratios for each peptide of non-reference channels
        ConsensusFeature::HandleSetType::iterator ref_it;
        Peak2D::IntensityType ref_intensity;
        for (size_t i = 0; i < consensus_map_out.size(); ++i)
        {
          // find reference index (this is inefficient to do every time,
          // but the most robust against anyone who tries to change the internals of ConsensusFeature):
          ref_it = consensus_map_out[i].end();
          for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
               it_elements != consensus_map_out[i].end();
               ++it_elements)
          {
            if ((Int) consensus_map_out.getFileDescriptions()[it_elements->getMapIndex()].getMetaValue("channel_name") == reference_channel)
            {
              ref_it = it_elements;
              break;
            }
          }

          // reference channel not found in this ConsensusFeature
          if (ref_it == consensus_map_out[i].end())
          {
            LOG_ERROR << "ItraqQuantifier::run() WARNING: ConsensusFeature " << i << " does not have a reference channel! Skipping" << std::endl;
            continue;
          }

          ref_intensity = ref_it->getIntensity();

          // now collect the ratios and intensities
          for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
               it_elements != consensus_map_out[i].end();
               ++it_elements)
          {
            if (ref_intensity == 0) //avoid nan's and inf's
            {
              if (it_elements->getIntensity() == 0) // 0/0 will give 'nan'
              {
                //so leave it out completely (there is no information to be gained)
              }
              else                  // x/0 is 'inf' but std::sort() has problems with that
              {
                peptide_ratios[map_to_vectorindex[it_elements->getMapIndex()]].push_back(std::numeric_limits<double>::max());
              }
            }
            else             // everything seems fine
            {
              peptide_ratios[map_to_vectorindex[it_elements->getMapIndex()]].push_back(it_elements->getIntensity() / ref_intensity);
            }

            // control
            peptide_intensities[map_to_vectorindex[it_elements->getMapIndex()]].push_back(it_elements->getIntensity());
          }
        } // ! collect ratios

        double max_deviation_from_control = 0;
        // find MEDIAN of ratios for each channel (store as 0th element in sorted vector)
        for (Map<Size, Size>::const_iterator it_map = map_to_vectorindex.begin(); it_map != map_to_vectorindex.end(); ++it_map)
        {
          // sort vector (partial_sort might improve performance here)
          std::sort(peptide_ratios[it_map->second].begin(), peptide_ratios[it_map->second].end());
          // save median as first element
          peptide_ratios[it_map->second][0] = peptide_ratios[it_map->second][peptide_ratios[it_map->second].size() / 2];

          // sort control (intensities)
          std::sort(peptide_intensities[it_map->second].begin(), peptide_intensities[it_map->second].end());
          // find MEDIAN of control-method (intensities) for each channel
          peptide_intensities[it_map->second][0] = peptide_intensities[it_map->second][peptide_intensities[it_map->second].size() / 2] /
                                                   peptide_intensities[ref_mapid][peptide_intensities[ref_mapid].size() / 2];
          //#ifdef ITRAQ_DEBUG
          LOG_INFO << "iTRAQ-normalize:  map-id " << (it_map->first) << " has factor " << (peptide_ratios[it_map->second][0]) << " (control: " << (peptide_intensities[it_map->second][0]) << ")" << std::endl;
          //#endif
          double dev = (peptide_ratios[it_map->second][0] - peptide_intensities[it_map->second][0]) / peptide_ratios[it_map->second][0];
          if (fabs(max_deviation_from_control) < fabs(dev))
          {
            max_deviation_from_control = dev;
          }
        }

        LOG_INFO << "iTRAQ-normalization: max ratio deviation of alternative method is " << (max_deviation_from_control * 100) << "%\n";

#ifdef ITRAQ_DEBUG
        std::cout << "debug OUTPUT\n";
        for (Size i = 1; i < peptide_ratios[0].size(); ++i)
        {
          if (i == peptide_intensities[0].size() / 2)
          {
            std::cout << "++++++++++ median: \n";
          }
          for (Size j = 0; j < peptide_ratios.size(); ++j)
          {
            std::cout << peptide_ratios[j][i] << " ";
          }
          std::cout << " -- int -- ";
          for (Size j = 0; j < peptide_intensities.size(); ++j)
          {
            std::cout << peptide_intensities[j][i] << " ";
          }
          if (i == peptide_intensities[0].size() / 2)
          {
            std::cout << "\n----------- median: ";
          }
          std::cout << "\n";
        }
#endif

        // adjust intensity ratios
        for (size_t i = 0; i < consensus_map_out.size(); ++i)
        {
          // find reference index (this is inefficient to do every time,
          // but the most robust against anyone who tries to change the internals of ConsensusFeature):
          ref_it = consensus_map_out[i].end();
          for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
               it_elements != consensus_map_out[i].end();
               ++it_elements)
          {
            if ((Int) consensus_map_out.getFileDescriptions()[it_elements->getMapIndex()].getMetaValue("channel_name") == reference_channel)
            {
              ref_it = it_elements;
              break;
            }
          }

          // reference channel not found in this ConsensusFeature
          if (ref_it == consensus_map_out[i].end())
          {
            continue;
          }

          ref_intensity = ref_it->getIntensity();

          // now adjust the ratios
          ConsensusFeature cf = consensus_map_out[i];
          cf.clear(); // delete its handles
          for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
               it_elements != consensus_map_out[i].end();
               ++it_elements)
          {
            FeatureHandle hd = *it_elements;
            if (it_elements == ref_it)
            {
              hd.setIntensity(1);
            }
            else          // divide current intensity by normalization factor (which was stored at position 0)
            {
              hd.setIntensity(hd.getIntensity() / peptide_ratios[map_to_vectorindex[it_elements->getMapIndex()]][0]);
            }
            cf.insert(hd);
          }
          // replace consensusFeature with updated intensity
          consensus_map_out[i] = cf;
        } // ! adjust ratios

      } // ! ref_channel valid
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier::run() Parameter 'channel_reference' does not name a valid channel!");
      }
    } // !do_normalization


    // ** PEPTIDE PROTEIN MAPPING ** //

    consensus_map_out.setExperimentType("itraq");

    return;
  }

  ItraqQuantifier::ItraqQuantifierStats ItraqQuantifier::getStats() const
  {
    return stats_;
  }

  void ItraqQuantifier::setDefaultParams_()
  {
    // choose default and documentation depending on itraq/tmt, since we can provide no stable default for TMT
    defaults_.setValue("isotope_correction",
                       (itraq_type_ == TMT_SIXPLEX
                        ? "false"
                        : "true"),
                       (itraq_type_ == TMT_SIXPLEX
                        ? "Enable isotope correction (highly recommended). Note that you need to provide a correction matrix (see isotope_correction:tmt-6plex otherwise the tool will fail."
                        : "Enable isotope correction (highly recommended)."),
                       ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("isotope_correction", ListUtils::create<String>("true,false"));

    defaults_.setValue("do_normalization", "false", "Normalize channels? Done by using the Median of Ratios (every channel / Reference). Also the ratio of medians (from any channel and reference) is provided as control measure!", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("do_normalization", ListUtils::create<String>("true,false"));

    if (itraq_type_ == TMT_SIXPLEX)
    {
      defaults_.setValue("isotope_correction:tmt-6plex",
                         ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::TMT_SIXPLEX, isotope_corrections_),
                         "Override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '126:0/0.3/4/0' , '128:0.1/0.3/3/0.2'.",
                         ListUtils::create<String>("advanced"));
    }
    else
    {
      defaults_.setValue("isotope_correction:4plex",
                         ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::FOURPLEX, isotope_corrections_),
                         "Override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2'.",
                         ListUtils::create<String>("advanced"));
      defaults_.setValue("isotope_correction:8plex",
                         ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::EIGHTPLEX, isotope_corrections_),
                         "Override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2'.",
                         ListUtils::create<String>("advanced"));
    }

    defaults_.setSectionDescription("isotope_correction",
                                    (itraq_type_ == TMT_SIXPLEX
                                     ? "Isotope correction matrices for tmt-6plex."
                                     : "Isotope correction matrices for 4plex and 8plex. Only one of them will be used (depending on iTRAQ mode)."));


    // for 4 & 8 plex. Max value is again checked during runtime
    defaults_.setValue("channel_reference",
                       (itraq_type_ != TMT_SIXPLEX
                        ? 114
                        : 126),
                       (itraq_type_ != TMT_SIXPLEX
                        ? "Number of the reference channel (114-117 for 4plex)."
                        : "Number of the reference channel (126-131)."));
    if (itraq_type_ == TMT_SIXPLEX)
    {
      defaults_.setMinInt("channel_reference", 126);
      defaults_.setMaxInt("channel_reference", 131);
    }
    else if (itraq_type_ == FOURPLEX)
    {
      defaults_.setMinInt("channel_reference", 114);
      defaults_.setMaxInt("channel_reference", 117);
    }
    else // EIGHTPLEX
    {
      defaults_.setMinInt("channel_reference", 113);
      defaults_.setMaxInt("channel_reference", 121);
    }

    defaultsToParam_();
  }

  void ItraqQuantifier::updateMembers_()
  {
    StringList channels;
    // update isotope_corrections_ Matrix with custom values
    if (itraq_type_ == ItraqConstants::FOURPLEX)
    {
      channels = param_.getValue("isotope_correction:4plex");
    }
    else if (itraq_type_ == ItraqConstants::EIGHTPLEX)
    {
      channels = param_.getValue("isotope_correction:8plex");
    }
    else if (itraq_type_ == ItraqConstants::TMT_SIXPLEX)
    {
      channels = param_.getValue("isotope_correction:tmt-6plex");
    }

    if (channels.size() > 0)
    {
      ItraqConstants::updateIsotopeMatrixFromStringList(itraq_type_, channels, isotope_corrections_);
    }
  }

  /// initialize
  void ItraqQuantifier::initIsotopeCorrections_()
  {
    isotope_corrections_.resize(3);
    isotope_corrections_[0].setMatrix<4, 4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
    isotope_corrections_[1].setMatrix<8, 4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
    isotope_corrections_[2].setMatrix<6, 4>(ItraqConstants::ISOTOPECORRECTIONS_TMT_SIXPLEX);
  }

  /// extract channel information (active channels, names, etc) from ConsensusMap
  void ItraqQuantifier::reconstructChannelInfo_(const ConsensusMap& consensus_map)
  {
    channel_map_.clear();

    for (ConsensusMap::FileDescriptions::const_iterator file_it = consensus_map.getFileDescriptions().begin();
         file_it != consensus_map.getFileDescriptions().end();
         ++file_it)
    {
      if (file_it->second.metaValueExists("channel_name"))
      {
        ChannelInfo info;
        // fill info
        info.name = file_it->second.getMetaValue("channel_name");
        info.id = file_it->second.getMetaValue("channel_id");
        info.description = file_it->second.getMetaValue("channel_description");
        info.center = file_it->second.getMetaValue("channel_center");
        info.active = (String(file_it->second.getMetaValue("channel_active")) == "true" ? true : false);
        channel_map_[info.name] = info;
#ifdef ITRAQ_DEBUG
        std::cout << " setting info.name " << (info.name) << " and id " << (info.id) << std::endl;
#endif
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier::reconstructChannelInfo_ The ConsensusMap provided is missing MetaInfo from ItraqChannelExtractor!");
      }
    }
  }

  std::ostream& operator<<(std::ostream& os, const ItraqQuantifier::ItraqQuantifierStats& stats)
  {
    os << "name\tvalue\t(value in %)\n";
    os << "# channels\t" << stats.channel_count << "\tNA\n";
    os << "# spectra total\t" << stats.number_ms2_total << "\tNA\n";
    os << "# spectra negative\t" << stats.iso_number_reporter_negative << "\tNA\n";
    os << "# negative reporter intensity\t" << stats.iso_number_reporter_negative << "\tNA\n";
    os << "# alternative positive reporter intensity\t" << stats.iso_number_reporter_different << "\tNA\n";
    os << "total intensity (affected spectra)\t" << stats.iso_total_intensity_negative << "\tNA\n";
    os << "total intensity difference (affected spectra)\t" << stats.iso_solution_different_intensity << "\t" << (stats.iso_solution_different_intensity * 100 / stats.iso_total_intensity_negative) << "\n";

    for (std::map<Size, Size>::const_iterator it_m = stats.empty_channels.begin(); it_m != stats.empty_channels.end(); ++it_m)
    {
      os << "labeling_efficiency_channel_" << it_m->first << "\t" << (stats.number_ms2_total - it_m->second) << "\t" << ((stats.number_ms2_total - it_m->second) * 100 / stats.number_ms2_total) << "\n";
    }

    return os;
  }

}
