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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/EnhancedTabBarWidgetInterface.h>

#include <OpenMS/VISUAL/EnhancedTabBar.h>

namespace OpenMS
{

  EnhancedTabBarWidgetInterface::EnhancedTabBarWidgetInterface()
  {
    /// every new window gets a new ID automatically
    static Int window_counter_ = getFirstWindowID();
    window_id_ = ++window_counter_;
  }

  EnhancedTabBarWidgetInterface::~EnhancedTabBarWidgetInterface()
  { // remove ourselves from the tab bar
    if (parent_) parent_->removeId(window_id_);
  }

  void EnhancedTabBarWidgetInterface::addToTabBar(EnhancedTabBar* const parent, const String& caption, const bool make_active_tab)
  {
    parent_ = parent;
    parent_->addTab(caption.toQString(), window_id_);
    if (make_active_tab) parent_->setCurrentId(window_id_);
  }

  Int EnhancedTabBarWidgetInterface::getWindowId()
  {
    return window_id_;
  }

  /*static*/ Int EnhancedTabBarWidgetInterface::getFirstWindowID()
  {
    return 1234;
  }

}

