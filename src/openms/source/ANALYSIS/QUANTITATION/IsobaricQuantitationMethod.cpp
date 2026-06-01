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

} // namespace
