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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_APPLICATIONS_MISC_QAPPLICATIONTOPP_H
#define OPENMS_VISUAL_APPLICATIONS_MISC_QAPPLICATIONTOPP_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//Qt
#include <QtGui/QApplication>

namespace OpenMS
{
  /**
    @brief Extension to the QApplication for running TOPPs GUI tools.

    Basically re-implements notify of QApplication to prevent ungraceful exit.
  */
  class OPENMS_GUI_DLLAPI QApplicationTOPP :
    public QApplication
  {

    Q_OBJECT

public:
    /// Constructor (no NOT remove the "&" from argc, since Qt will segfault on some platforms otherwise!)
    QApplicationTOPP(int& argc, char** argv);

    /// Destructor
    ~QApplicationTOPP() override;

    /*
      @brief: Catch exceptions in Qt GUI applications, preventing ungraceful exit

      Re-implementing QApplication::notify() to catch exception thrown in event
      handlers (which is most likely OpenMS code).
    */
    bool notify(QObject* rec, QEvent* ev) override;

    /*
      Reimplemented from QApplication, to handle QEvent::FileOpen to enable handling of odoc event on MacOSX
    */
    bool event(QEvent*) override;

signals:
    void fileOpen(QString file);

  };

}

#endif
