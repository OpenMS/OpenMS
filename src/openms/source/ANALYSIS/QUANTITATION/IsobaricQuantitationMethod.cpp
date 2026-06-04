// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <array>

namespace OpenMS
{
  IsobaricQuantitationMethod::~IsobaricQuantitationMethod() = default;

  IsobaricQuantitationMethod::IsobaricQuantitationMethod() :
    DefaultParamHandler("IsobaricQuantitationMethod")
  {
  }

  Matrix<double> IsobaricQuantitationMethod::stringListToIsotopeCorrectionMatrix_(const StringList& stringlist) const
  {
    // check the string list
    if (stringlist.size() != getNumberOfChannels())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("IsobaricQuantitationMethod: Invalid string representation of the isotope correction matrix. Expected ") + getNumberOfChannels() + " entries but got " + stringlist.size() + ".");
    }

    // compute frequency matrix based on the deviation matrix
    Matrix<double> channel_frequency(getNumberOfChannels(), getNumberOfChannels(), 0.0);

    // channel index
    Size contributing_channel = 0;

    // fill row-wise
    for (const auto& l : stringlist)
    {
      StringList corrections;
      l.split('/', corrections);

      auto number_of_columns = getChannelInformation()[contributing_channel].affected_channels.size();
      if (corrections.size() != number_of_columns )
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Corrections for channel ID #") + contributing_channel + " must contain " + number_of_columns + " values, but has " + corrections.size() + "!", String(corrections.size()));
      }

      // overwrite line in Matrix with custom values
      Size affected_channel_idx = 0;
      double self_contribution = 100.0;
      double correction;
      Int target_channel;
      for (auto& c : corrections)
      {
        c = c.trim().toUpper();
        if (c != "NA" && c != "-1" && c != "0.0")
        {
          target_channel = getChannelInformation()[contributing_channel].affected_channels[affected_channel_idx];
          try
          {
            correction = c.toDouble();
          }
          catch (Exception::ConversionError& e)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Correction entry #") + affected_channel_idx + " in channel ID " + contributing_channel + " must be one of na/NA/-1 or a floating point number representation!", c);
          }

          if (correction < 0.0 || correction > 100.0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Correction entry #") + affected_channel_idx + " in channel ID " + contributing_channel + " must be a percentage between 0 and 100!", c);
          }
          
          if (target_channel >= 0 && Size(target_channel) < getNumberOfChannels())
          {
            channel_frequency(target_channel, contributing_channel) = correction / 100.0;
          }
          self_contribution -= correction; // count reduced self-contribution even if it does not affect another channel
        }
        affected_channel_idx++;
      }
      // set reduced self contribution
      channel_frequency(contributing_channel, contributing_channel) = self_contribution / 100.0;
      // increment channel index
      ++contributing_channel;
    }
    return channel_frequency;
  }

  Matrix<double> IsobaricQuantitationMethod::stringListToIsotopeCorrectionMatrixSplit_(const std::vector<String>& nondeuterated,
                                                                                       const std::vector<String>& deuterated) const
  {
    // Positions in the 14-column layout that the eight non-deuterated values map onto.
    //  8-col order: <-2C13>/<-N15-C13>/<-C13>/<-N15>/<+N15>/<+C13>/<+N15+C13>/<+2C13>
    // 14-col order: <-C13-H2>/<-2C13>/<-N15-H2>/<-C13-N15>/<-H2>/<-C13>/<-N15>/<+N15>/<+C13>/<+H2>/<+N15+C13>/<+N15+H2>/<+2C13>/<+C13+H2>
    // The six remaining positions are the 2H (deuterium) offsets, which are 'NA' for non-deuterated reagents.
    static const std::array<Size, 8> non_deut_to_full = {{1, 3, 5, 6, 7, 8, 10, 12}};
    const Size full_columns = 14;

    const IsobaricChannelList& channels = getChannelInformation();

    // count how many deuterated / non-deuterated channels this method expects
    Size n_deut = 0;
    for (const auto& ch : channels)
    {
      if (ch.name.find('D') != std::string::npos) { ++n_deut; }
    }
    const Size n_nondeut = getNumberOfChannels() - n_deut;

    if (nondeuterated.size() != n_nondeut)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("IsobaricQuantitationMethod: Expected ") + n_nondeut + " non-deuterated correction matrix entries but got " + nondeuterated.size() + ".");
    }
    if (deuterated.size() != n_deut)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("IsobaricQuantitationMethod: Expected ") + n_deut + " deuterated correction matrix entries but got " + deuterated.size() + ".");
    }

    // interleave both inputs into one 14-column row per channel, in channel order
    std::vector<String> full;
    full.reserve(getNumberOfChannels());
    Size i_nondeut = 0;
    Size i_deut = 0;
    for (const auto& ch : channels)
    {
      if (ch.name.find('D') != std::string::npos)
      {
        // deuterated channel: already provided in the 14-column layout
        full.push_back(deuterated[i_deut++]);
        continue;
      }

      // non-deuterated channel: expand the 8-column row into the 14-column layout
      std::vector<String> values;
      nondeuterated[i_nondeut++].split('/', values);
      if (values.size() != non_deut_to_full.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Non-deuterated correction entry for channel '") + ch.name + "' must contain " + non_deut_to_full.size() + " values, but has " + values.size() + "!", String(values.size()));
      }
      std::vector<String> row(full_columns, "NA");
      for (Size k = 0; k < non_deut_to_full.size(); ++k)
      {
        row[non_deut_to_full[k]] = values[k].trim();
      }
      String joined;
      for (Size k = 0; k < full_columns; ++k)
      {
        if (k != 0) { joined += "/"; }
        joined += row[k];
      }
      full.push_back(joined);
    }

    return stringListToIsotopeCorrectionMatrix_(full);
  }

} // namespace
