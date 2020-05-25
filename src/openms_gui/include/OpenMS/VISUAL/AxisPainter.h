// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <vector>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QPaintEvent>
#include <QPainter>


namespace OpenMS
{
  /**
    @brief Draws a coordinate axis.
    It has only static methods, that's why the constructor is private.
    @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI AxisPainter
  {
public:
    /// Typedef for the grid vector
    typedef std::vector<std::vector<double> > GridVector;

    /// Where the axis is placed
    enum Alignment
    {
      TOP,
      BOTTOM,
      LEFT,
      RIGHT
    };

    /// Draws an axis
    static void paint(QPainter * painter, QPaintEvent * e, const double & min, const double & max, const GridVector & grid,
                      const Int width, const Int height, const Alignment alignment, const UInt margin,
                      const bool show_legend, const String legend, const bool shorten_number,
                      const bool is_log, const bool is_inverse_orientation);
private:
    /// Constructor: only static methods
    AxisPainter();

    /// sets @p short_num to a shortened string representation ("123.4 k/M/G") of @p number
    static void getShortenedNumber_(QString & short_num, double number);

    /// Scale axis values to correct value (i.e. reverse log, unit conversion)
    static double scale_(double x, bool is_log);
  };
}
