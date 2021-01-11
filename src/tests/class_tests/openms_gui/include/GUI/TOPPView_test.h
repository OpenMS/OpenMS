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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef TEST_GUI_TESTTOPPVIEW_H
#define TEST_GUI_TESTTOPPVIEW_H


#include <QtTest/QtTest>
#include <QtGui>
#include <QQueue>


#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  class String;
/**

  @todo write a proper GUI base class for the scheduler below (Chris)

*/
class TestTOPPView: public QObject
{
  Q_OBJECT
 
  /**
    @brief Store information on timed keyboard input events

    Store the time offset, the key sequence and the expected window title.

  */
  struct ScheduleInfo
  {
    ScheduleInfo() :
      delay(0)
    {}

    ScheduleInfo(QString p_keys, QString p_title, int p_delay) :
      keys(p_keys),
      title(p_title),
      delay(p_delay)
    {}

    QString keys; //< key sequence
    QString title;//< expected window title
    int delay;    //< delay in ms when event is fired off
  };

  public slots:

  /**
    @brief Slot that tries to process the current event queue for modal dialogs until its empty.

    Slot that tries to process the current event queue for modal dialogs until its empty.
    The function will repeatedly invoke itself until the queue is empty, to allow other
    incoming events (e.g. loading a file) to be processed in between two scheduled dialogs.
  */
  void simulateClick_();

  private slots:
    void testGui();

  private:
    /**
      @brief Schedule a keyboard input using a QTimer signal to direct input to a modal window.

      Modal windows have their own event queue and once launched will halt the
      execution of the test script until closed. This implies one cannot simply direct keyboard
      input to them. To do that we pre-schedule the input in the main event loop
      using a timer. The keyboard input sequence @p key_sequence and a subsequent 'return' key press
      is then issued when the time out occurs, given that the current window has the correct
      @p title! Otherwise the event is rescheduled until the title is correct.
      The @p delay then the timer pops is relative to the last timed events successfull completion.
    */
    void scheduleModalWidget_(const QString& key_sequence, const QString& title, const int delay=100);

    /**
      @brief Waits until the scheduled event queue is emtpy

      Waits until the scheduled event queue is emtpy.
    */
    void waitForModalWidget(const int max_wait, const String& line);

    /// event queue for modal/popup dialogs
    QQueue<ScheduleInfo> modal_key_sequence_;
};

}

#endif

