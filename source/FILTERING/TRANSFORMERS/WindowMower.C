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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>

using namespace std;
namespace OpenMS
{
  WindowMower::WindowMower() :
    DefaultParamHandler("WindowMower")
  {
    defaults_.setValue("windowsize", 50.0, "The size of the sliding window along the m/z axis.");
    defaults_.setValue("peakcount", 2, "The number of peaks that should be kept.");
    defaults_.setValue("movetype", "slide", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    defaults_.setValidStrings("movetype", StringList::create("slide,jump"));
    defaultsToParam_();
  }

  WindowMower::~WindowMower()
  {
  }

  WindowMower::WindowMower(const WindowMower & source) :
    DefaultParamHandler(source)
  {
  }

  WindowMower & WindowMower::operator=(const WindowMower & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void WindowMower::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    bool sliding = (String)param_.getValue("movetype") == "slide" ? true : false;

    if (sliding)
    {
      filterPeakSpectrumForTopNInSlidingWindow(spectrum);
    } else
    {
      filterPeakSpectrumForTopNInJumpingWindow(spectrum);
    }
  }

  void WindowMower::filterPeakMap(PeakMap & exp)
  {
    bool sliding = (String)param_.getValue("movetype") == "slide" ? true : false;
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      if (sliding)
      {
        filterPeakSpectrumForTopNInSlidingWindow(*it);
      } else
      {
        filterPeakSpectrumForTopNInJumpingWindow(*it);
      }
    }
  }

}
