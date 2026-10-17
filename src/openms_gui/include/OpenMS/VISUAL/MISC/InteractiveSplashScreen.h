// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtWidgets/QSplashScreen>

#include <memory>

// forward declarations (Qt classes live in the global namespace)
class QKeyEvent;
class QMouseEvent;
class QPixmap;

namespace OpenMS
{
  /**
    @brief A QSplashScreen which stays responsive and can be dismissed early by the user.

    Applications typically want to show a splash screen for a minimum duration (so users can
    actually read it).

    This class instead offers showFor(), which runs a nested, interruptible event loop. The
    splash stays responsive for at most the requested time and closes immediately once the user
    clicks it or presses a key.

    Usage:
    @code
    InteractiveSplashScreen splash(pixmap);
    splash.show();
    // load data, optionally calling splash.showMessage(...)
    // ...
    splash.showFor(3.0); // keep it up for <= 3s since construction, but let the user dismiss it early
    @endcode
  */
  class OPENMS_GUI_DLLAPI InteractiveSplashScreen : public QSplashScreen
  {
  public:
    /// C'tor from a pixmap (as for QSplashScreen). Starts an internal timer used by showFor().
    explicit InteractiveSplashScreen(const QPixmap& pixmap);

    /// D'tor
    ~InteractiveSplashScreen() override;

    /**
      @brief Keep the splash visible and responsive for at most @p max_seconds, then close it.

      Returns early (closing the splash) as soon as the user clicks the splash or presses a key.
      The elapsed time is counted from construction, so any work done between constructing the
      splash and calling this method counts towards @p max_seconds.
    */
    void showFor(double max_seconds);

  protected:
    void mousePressEvent(QMouseEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
} // namespace OpenMS
