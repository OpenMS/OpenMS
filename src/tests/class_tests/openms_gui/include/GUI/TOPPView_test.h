// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

