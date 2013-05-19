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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

namespace OpenMS
{
  IsobaricQuantitationMethod::~IsobaricQuantitationMethod()
  {
  }

  IsobaricQuantitationMethod::IsobaricQuantitationMethod() :
    DefaultParamHandler("IsobaricQuantitationMethod")
  {
  }

  Matrix<double> IsobaricQuantitationMethod::stringListToIsotopCorrectionMatrix_(const StringList& stringlist) const
  {
    // check the string list
    if (stringlist.size() != getNumberOfChannels())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("IsobaricQuantitationMethod: Invalid string representation of the isotope correction matrix. Expected ") + getNumberOfChannels() + " entries but got " + stringlist.size() + ".");
    }

    // create matrix to fill from stringlist
    Matrix<double> isotope_correction_matrix(getNumberOfChannels(), 4);

    // channel index
    Size channel_index = 0;

    // fill row-wise
    for (StringList::ConstIterator it = stringlist.begin(); it != stringlist.end(); ++it)
    {
      StringList corrections;
      it->split('/', corrections);
      if (corrections.size() != 4)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IsobaricQuantitationMethod: Invalid entry in string representation of the isotope correction matrx. Expected four correction values separated by '/', got: '" + *it + "'");
      }

      // overwrite line in Matrix with custom values
      isotope_correction_matrix.setValue(channel_index, 0, corrections[0].toDouble());
      isotope_correction_matrix.setValue(channel_index, 1, corrections[1].toDouble());
      isotope_correction_matrix.setValue(channel_index, 2, corrections[2].toDouble());
      isotope_correction_matrix.setValue(channel_index, 3, corrections[3].toDouble());

      // increment channel index
      ++channel_index;
    }

    // compute frequency matrix based on the deviation matrix
    Matrix<double> channel_frequency(getNumberOfChannels(), getNumberOfChannels());

    // matrix element i, j == what contributes channel i to the intensity of channel j
    for (Size contributing_channel = 0; contributing_channel < getNumberOfChannels(); ++contributing_channel)
    {
      for (Size target_channel = 0; target_channel < getNumberOfChannels(); ++target_channel)
      {
        // as the isotope_correction_matrix encodes what channel i contributes to theoretical -2/-1/+1/+2 channels
        // hence we check if the target_channel corresponds to the contributing_channel -2/-1/+1/+2
        // if yes, we assign the corresponding contribution
        // contribution to itself is handled separately
        if ((getChannelInformation()[contributing_channel].name - 2) == getChannelInformation()[target_channel].name)
        {
          channel_frequency.setValue(target_channel, contributing_channel, isotope_correction_matrix.getValue(contributing_channel, 0) / 100);
        }
        else if ((getChannelInformation()[contributing_channel].name - 1) == getChannelInformation()[target_channel].name)
        {
          channel_frequency.setValue(target_channel, contributing_channel, isotope_correction_matrix.getValue(contributing_channel, 1) / 100);
        }
        else if ((getChannelInformation()[contributing_channel].name + 1) == getChannelInformation()[target_channel].name)
        {
          channel_frequency.setValue(target_channel, contributing_channel, isotope_correction_matrix.getValue(contributing_channel, 2) / 100);
        }
        else if ((getChannelInformation()[contributing_channel].name + 2) == getChannelInformation()[target_channel].name)
        {
          channel_frequency.setValue(target_channel, contributing_channel, isotope_correction_matrix.getValue(contributing_channel, 3) / 100);
        }
        else if (target_channel == contributing_channel)
        {
          double self_contribution = 100.0;
          for (Size column_idx = 0; column_idx < 4; ++column_idx)
          {
            self_contribution -= isotope_correction_matrix.getValue(contributing_channel, column_idx);
          }
          channel_frequency.setValue(contributing_channel, contributing_channel, (self_contribution / 100));
        }
      }
    }

    return channel_frequency;
  }

} // namespace
