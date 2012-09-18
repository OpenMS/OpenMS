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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>


namespace OpenMS
{

  // number of channels for iTRAQ types. (make sure it corresponds to enum ITRAQ_TYPES)
  const Int ItraqConstants::CHANNEL_COUNT[3] = {4, 8, 6};

  const Int ItraqConstants::CHANNELS_FOURPLEX[4][1] = {{114}, {115}, {116}, {117}};
  const Int ItraqConstants::CHANNELS_EIGHTPLEX[8][1] = {{113}, {114}, {115}, {116}, {117}, {118}, {119}, {121}};
  const Int ItraqConstants::CHANNELS_TMT_SIXPLEX[6][1] = {{126}, {127}, {128}, {129}, {130}, {131}};

  // currently from http://www.matrixscience.com/help/quant_config_help.html
  // (@ 117; +2 the value is 0.1 not 0.0 as confirmed by ABSciex)
  const double ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX[4][4] =
  {
    {0.0, 1.0, 5.9, 0.2},   //114
    {0.0, 2.0, 5.6, 0.1},
    {0.0, 3.0, 4.5, 0.1},
    {0.1, 4.0, 3.5, 0.1}    //117
  };

  //taken from Applied Biosystems Website
  // http://www.absciex.com/Documents/Support/AB_SCIEX_Question_and_Answer.xls
  const double ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX[8][4] =
  {
    {0.00, 0.00, 6.89, 0.22},   //113
    {0.00, 0.94, 5.90, 0.16},
    {0.00, 1.88, 4.90, 0.10},
    {0.00, 2.82, 3.90, 0.07},
    {0.06, 3.77, 2.99, 0.00},
    {0.09, 4.71, 1.88, 0.00},
    {0.14, 5.66, 0.87, 0.00},
    {0.27, 7.44, 0.18, 0.00}    //121
  };

  //taken from ThermoFischer Scientific
  // http://www.piercenet.com/coapdfs/CofA-90064-SPECS.pdf
  const double ItraqConstants::ISOTOPECORRECTIONS_TMT_SIXPLEX[6][4] =
  {
    {0.00, 0.00, 0.00, 0.00},      //126
    {0.00, 0.00, 0.00, 0.00},      //127
    {0.00, 0.00, 0.00, 0.00},      //128
    {0.00, 0.00, 0.00, 0.00},      //129
    {0.00, 0.00, 0.00, 0.00},      //130
    {0.00, 0.00, 0.00, 0.00},      //131
  };

  StringList ItraqConstants::getIsotopeMatrixAsStringList(const int itraq_type, const IsotopeMatrices & isotope_corrections)
  {
    OPENMS_PRECONDITION(itraq_type < SIZE_OF_ITRAQ_TYPES && itraq_type >= 0, "Error while trying to access invalid isotope correction matrix.");

    StringList isotopes;
    std::vector<Matrix<Int> > channel_names(3);
    channel_names[0].setMatrix<4, 1>(CHANNELS_FOURPLEX);
    channel_names[1].setMatrix<8, 1>(CHANNELS_EIGHTPLEX);
    channel_names[2].setMatrix<6, 1>(CHANNELS_TMT_SIXPLEX);
    for (Int i = 0; i < CHANNEL_COUNT[itraq_type]; ++i)
    {
      String line = String(channel_names[itraq_type].getValue(i, 0)) + ":";
      for (Int j = 0; j < 3; ++j)
      {
        line += String(isotope_corrections[itraq_type].getValue(i, j)) + "/";
      }
      line += String(isotope_corrections[itraq_type].getValue(i, 3));
      isotopes.push_back(line);
    }

    return isotopes;
  }

  void ItraqConstants::updateIsotopeMatrixFromStringList(const int itraq_type, const StringList & channels, IsotopeMatrices & isotope_corrections)
  {
    // TODO: make generic .. why do we need to initialize all matrices, we are only interested in itraq_type
    isotope_corrections.resize(SIZE_OF_ITRAQ_TYPES);
    isotope_corrections[0].setMatrix<4, 4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
    isotope_corrections[1].setMatrix<8, 4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
    isotope_corrections[2].setMatrix<6, 4>(ItraqConstants::ISOTOPECORRECTIONS_TMT_SIXPLEX);

    // split the channels key:name pairs apart
    for (StringList::const_iterator it = channels.begin(); it != channels.end(); ++it)
    {
      StringList result;
      it->split(':', result);
      if (result.size() != 2)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; expected one ':', got this: '" + (*it) + "'");
      }
      result[0] = result[0].trim();       // hold channel name
      result[1] = result[1].trim();       // holds 4 values

      Int channel = result[0].toInt();
      Int line = 0;
      if (itraq_type == FOURPLEX)
        line = channel - 114;
      else if (itraq_type == EIGHTPLEX)
        line = channel - 113;
      else // TODO: what do we need as offset here
        line = channel - 126;

      if ((itraq_type == FOURPLEX && (line < 0 || line > 3))
         ||
          ((itraq_type == EIGHTPLEX && (line < 0 || line > 8)) || channel == 120)
         ||
          (itraq_type == TMT_SIXPLEX && (line < 0 || line > 5)))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; channel-name is not valid for ")
                                          + String(
                                            itraq_type == FOURPLEX ? "4plex" : (itraq_type == EIGHTPLEX ? "8plex" : "TMT-6plex")
                                            )
                                          + String(": '")
                                          + result[0]
                                          + String("'"));
      }

      // if set to 121 we still want to change line 7 of the matrix
      if (line == 8 && itraq_type == EIGHTPLEX)
        line = 7;

      StringList corrections;
      result[1].split('/', corrections);
      if (corrections.size() != 4)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; expected four correction values separated by '&', got this: '" + result[1] + "'");
      }

      // overwrite line in Matrix with custom values
      isotope_corrections[itraq_type].setValue(line, 0, corrections[0].toDouble());
      isotope_corrections[itraq_type].setValue(line, 1, corrections[1].toDouble());
      isotope_corrections[itraq_type].setValue(line, 2, corrections[2].toDouble());
      isotope_corrections[itraq_type].setValue(line, 3, corrections[3].toDouble());

#ifdef ITRAQ_DEBUG
      std::cout << "Channel " << channel << " has values " << corrections << std::endl;
#endif
    }
  }

  void ItraqConstants::initChannelMap(const int itraq_type, ChannelMapType & map)
  {
    static Map<Int, DoubleReal> reporter_mass_exact;
    if (reporter_mass_exact.size() == 0 && (itraq_type == EIGHTPLEX || itraq_type == FOURPLEX)) // exact monoisotopic reporter ion masses (taken from AB Sciex)
    {
      reporter_mass_exact[113] = 113.1078;
      reporter_mass_exact[114] = 114.1112;
      reporter_mass_exact[115] = 115.1082;
      reporter_mass_exact[116] = 116.1116;
      reporter_mass_exact[117] = 117.1149;
      reporter_mass_exact[118] = 118.1120;
      reporter_mass_exact[119] = 119.1153;
      reporter_mass_exact[121] = 121.1220;
    }
    else // TMT(CID): exact masses taken from Thermo Scientific
    {
      reporter_mass_exact[126] = 126.127725;
      reporter_mass_exact[127] = 127.124760;
      reporter_mass_exact[128] = 128.134433;
      reporter_mass_exact[129] = 129.131468;
      reporter_mass_exact[130] = 130.141141;
      reporter_mass_exact[131] = 131.138176;
    }
    /// valid names for 4 and 8plex, i.e. 114,115,116,117 for 4plex
    std::vector<Matrix<Int> > channel_names;
    channel_names.resize(3);
    channel_names[0].setMatrix<4, 1>(CHANNELS_FOURPLEX);
    channel_names[1].setMatrix<8, 1>(CHANNELS_EIGHTPLEX);
    channel_names[2].setMatrix<6, 1>(CHANNELS_TMT_SIXPLEX);

    map.clear();
    for (Size i = 0; i < channel_names[itraq_type].rows(); ++i)
    {
      ChannelInfo info;
      info.description = "";
      info.name = channel_names[itraq_type].getValue(i, 0);
      info.id = (Int)i;
      if (!reporter_mass_exact.has(info.name))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unexpected reporter name during initialization.", String(info.name));
      }
      info.center = reporter_mass_exact[info.name];
      info.active = false;
      map[info.name] = info;
    }

#ifdef ITRAQ_DEBUG
    std::cout << "INIT: map has " << map.size() << " entries!" << std::endl;
#endif
  }

  void ItraqConstants::updateChannelMap(const StringList & active_channels, ChannelMapType & map)
  {
    // split the channels key:name pairs apart
    for (StringList::const_iterator it = active_channels.begin(); it != active_channels.end(); ++it)
    {
      StringList result;
      it->split(':', result);
      if (result.size() != 2)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqConstants: Invalid entry in Param 'channel_active'; expected one semicolon ('" + (*it) + "')");
      }
      result[0] = result[0].trim();
      result[1] = result[1].trim();
      if (result[0] == String::EMPTY || result[1] == String::EMPTY)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqConstants: Invalid entry in Param 'channel_active'; key or value is empty ('" + (*it) + "')");
      }
      Int channel = result[0].toInt();
      if (map.find(channel) == map.end())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqConstants: Invalid entry in Param 'channel_active'; channel is not valid ('" + String(channel) + "')");
      }
      // update name (description) of channel
      map[channel].description = result[1];
      map[channel].active = true;

#ifdef ITRAQ_DEBUG
      std::cout << "Channel " << channel << " has description " << channel_map_[channel].description << " and center " << channel_map_[channel].center << std::endl;
#endif
    }
  }

  Matrix<double> ItraqConstants::translateIsotopeMatrix(const int & itraq_type, const IsotopeMatrices & isotope_corrections)
  {
    // translate isotope_corrections to a channel_frequency matrix

    /*
      take special care of 8plex case, as ch-121 has a gap just before, thus matrix values need to be adapted
    */

    Matrix<double> channel_frequency(CHANNEL_COUNT[itraq_type], CHANNEL_COUNT[itraq_type]);
    for (Int i = 0; i < CHANNEL_COUNT[itraq_type]; ++i)
    {
      for (Int j = 0; j < CHANNEL_COUNT[itraq_type]; ++j)
      {
        // diagonal (should be close to 1 = 100%)
        if (i == j)
        {
          double val = 1.0;

          // subtract all isotope deviations of row i
          for (Int col_idx = 0; col_idx < 4; ++col_idx)
          {
            val += -isotope_corrections[itraq_type].getValue(i, col_idx) / 100;
          }
          channel_frequency.setValue(i, j, val);
        }
        else      // from mass i to mass j (directly copy the deviation)
        {
          if (i != 7 && j != 7)
          {
            if (j < i && i <= j + 2)      // -2, -1 cases of row 'i'
            {
              channel_frequency.setValue(j, i, isotope_corrections[itraq_type].getValue(i, j - i + 2) / 100);
            }
            else if (i < j && j <= i + 2)     // +1, +2 cases of row 'i'
            {
              channel_frequency.setValue(j, i, isotope_corrections[itraq_type].getValue(i, j - i + 1) / 100);
            }
          }
          else // special case of ch-121 for 8plex
          { // make everything more 'extreme' by 1 index
            if (i == 7 && j == 6) // -2 case of ch-121
            {
              channel_frequency.setValue(j, i, isotope_corrections[itraq_type].getValue(i, 0) / 100);
            }
            else if (i == 6 && j == 7) // +2 case of ch-121
            {
              channel_frequency.setValue(j, i, isotope_corrections[itraq_type].getValue(i, 3) / 100);
            }
          }
        }
      }
    }

    return channel_frequency;
  }

} // !namespace
