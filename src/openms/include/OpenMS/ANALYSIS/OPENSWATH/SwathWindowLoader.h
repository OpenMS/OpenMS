// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <vector>

namespace OpenMS
{

  /**
   * @brief Class to read a file describing the Swath Windows
   *
   * The file must of be tab delimited and of the following format:
   *    window_lower window_upper
   *    400 425
   *    425 450
   *    ...
   *
   * Note that the first line is a header and will be skipped.
   *
   */
  class OPENMS_DLLAPI SwathWindowLoader
  {

  public:

    /**
     * @brief Annotate a Swath map using a Swath window file specifying the individual windows
     *
     * @note It is assumed that the files in the swath_maps vector are in the
     * same order as the windows in the provided file (usually from lowest to
     * highest).
     *
     * @param filename The filename of the tab delimited file
     * @param swath_maps The list of SWATH maps (assumed to be in the same order as in the file)
     * @param doSort Sort the windows after reading in ascending order
     * @throw Exception::IllegalArgument if the maps in the file and the provided input maps to not match
     *
     */
    static void annotateSwathMapsFromFile(const std::string & filename,
                                          std::vector< OpenSwath::SwathMap >& swath_maps,
                                          bool doSort);

    /**
     * @brief Reading a tab delimited file specifying the SWATH windows
     *
     * The file must of be tab delimited and of the following format:
     *    window_lower window_upper
     *    400 425
     *    425 450
     *    ...
     *
     * Note that the first line is a header and will be skipped.
     *
     * @param filename The filename of the tab delimited file
     * @param swath_prec_lower_ The output vector for the window start
     * @param swath_prec_upper_ The output vector for the window end
     *
     */
    static void readSwathWindows(const std::string & filename,
                                 std::vector<double> & swath_prec_lower_,
                                 std::vector<double> & swath_prec_upper_ );
  };
}

