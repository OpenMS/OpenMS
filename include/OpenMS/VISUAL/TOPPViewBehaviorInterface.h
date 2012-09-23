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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPVIEWBEHAVIORINTERFACE_H
#define OPENMS_VISUAL_TOPPVIEWBEHAVIORINTERFACE_H

#include <QtCore/QObject>
#include <vector>

namespace OpenMS
{
  /** @brief Interface class to model different behaviors of TOPPView.

   @ingroup TOPPView
   */
  class TOPPViewBehaviorInterface :
    public QObject
  {
    Q_OBJECT

public:
    /// Destructor
    virtual ~TOPPViewBehaviorInterface() {}

public slots:
    /// Behavior for showSpectraumAs1D
    virtual void showSpectrumAs1D(int index) = 0;

    /// Behavior for activate1DSpectrum
    virtual void activate1DSpectrum(int index) = 0;
    virtual void activate1DSpectrum(std::vector<int, std::allocator<int> > indices) = 0;

    /// Behavior for deactivate1DSpectrum
    virtual void deactivate1DSpectrum(int index) = 0;

    /// Slot for behavior activation
    virtual void activateBehavior() = 0;

    /// Slot for behavior deactivation
    virtual void deactivateBehavior() = 0;
  };
}

#endif // OPENMS_VISUAL_TOPPVIEWBEHAVIORINTERFACE_H
