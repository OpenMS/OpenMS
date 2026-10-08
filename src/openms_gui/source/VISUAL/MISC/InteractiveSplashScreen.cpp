// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/MISC/InteractiveSplashScreen.h>

#include <OpenMS/SYSTEM/StopWatch.h>

#include <QtCore/QEventLoop>
#include <QtCore/QTimer>
#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>

namespace OpenMS
{
  /// Hidden state of InteractiveSplashScreen (PImpl), so the header stays free of these includes.
  struct InteractiveSplashScreen::Impl
  {
    StopWatch timer;             ///< measures time since construction
    QEventLoop* loop {nullptr};  ///< the nested loop running inside showFor(), or nullptr
    bool dismissed {false};      ///< set once the user has clicked / pressed a key

    /// mark as dismissed and break out of the nested event loop (if one is running)
    void dismiss()
    {
      dismissed = true;
      if (loop != nullptr)
      {
        loop->quit();
        loop = nullptr;
      }
    }
  };

  InteractiveSplashScreen::InteractiveSplashScreen(const QPixmap& pixmap)
    : QSplashScreen(pixmap), impl_(std::make_unique<Impl>())
  {
    // allow the splash to receive key events, so that any key press can dismiss it
    setFocusPolicy(Qt::StrongFocus);
    impl_->timer.start();
  }

  InteractiveSplashScreen::~InteractiveSplashScreen() = default;

  void InteractiveSplashScreen::showFor(double max_seconds)
  {
    // best effort to grab keyboard focus, so key presses are delivered to the splash
    raise();
    activateWindow();
    setFocus();

    const int remaining_ms = static_cast<int>((max_seconds - impl_->timer.getClockTime()) * 1000);
    if (!impl_->dismissed && remaining_ms > 0)
    {
      QEventLoop loop;
      impl_->loop = &loop;
      // quit once the minimum display time has elapsed ...
      QTimer::singleShot(remaining_ms, [&]() { impl_->dismiss(); });
      // ... or earlier, when the user dismisses the splash (see Impl::dismiss()).
      loop.exec();
      impl_->loop = nullptr; // 'loop' is about to be destroyed; clear the pointer so a late dismiss() never dereferences it
    }
    close();
  }

  void InteractiveSplashScreen::mousePressEvent(QMouseEvent* e)
  {
    impl_->dismiss();
    e->accept();
  }

  void InteractiveSplashScreen::keyPressEvent(QKeyEvent* e)
  {
    impl_->dismiss();
    e->accept();
  }
} // namespace OpenMS
