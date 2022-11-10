// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
            channel_frequency.setValue(target_channel, contributing_channel, correction / 100.0);
          }
          self_contribution -= correction; // count reduced self-contribution even if it does not affect another channel
        }
        affected_channel_idx++;
      }
      // set reduced self contribution
      channel_frequency.setValue(contributing_channel, contributing_channel, self_contribution / 100.0);
      // increment channel index
      ++contributing_channel;
    }
    return channel_frequency;
  }

} // namespace
