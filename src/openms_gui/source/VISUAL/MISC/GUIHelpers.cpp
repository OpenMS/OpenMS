// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDesktopServices>
#include <QGuiApplication>
#include <QMessageBox>
#include <QPainter>
#include <QPoint>
#include <QPointF>
#include <QProcess>
#include <QRectF>
#include <QString>
#include <QStringList>
#include <QtWidgets/QFileDialog>
#include <QUrl>

namespace OpenMS
{

  void GUIHelpers::openFolder(const QString& folder)
  {
#if defined(__APPLE__)
    QProcess p;
    p.setProcessChannelMode(QProcess::ForwardedChannels);
    p.start("/usr/bin/open", QStringList() << folder);
    if (!p.waitForStarted())
    {
      // execution failed
      QMessageBox::warning(0, "Open Folder Error", "The folder '" + folder + "' could not be opened!");
      OPENMS_LOG_ERROR << "Failed to open folder '" << folder.toStdString() << "'" << std::endl;
      OPENMS_LOG_ERROR << p.errorString().toStdString() << std::endl;
    }
#else
    if (!QDir(folder).exists() || (!QDesktopServices::openUrl(QUrl("file:///" + folder, QUrl::TolerantMode))))
    {
      QMessageBox::warning(nullptr, "Open Folder Error", "The folder '" + folder + "' could not be opened!");
    }
#endif
  }

  OPENMS_GUI_DLLAPI QString GUIHelpers::getSaveFilename(QWidget* parent, const QString& caption, const QString& dir, const FileTypeList& supported_file_types, bool add_all_filter,
                                                        const FileTypes::Type fallback_extension)
  {
    QString selected_filter;
    QString file_name = QFileDialog::getSaveFileName(parent, caption, dir, supported_file_types.toFileDialogFilter(FilterLayout::ONE_BY_ONE, add_all_filter).toQString(), &selected_filter);
    if (file_name.isEmpty())
    {
      return file_name;
    }
    // check whether a file type suffix has been given, or fall back to @p fallback_extension (if 'all filter' was used)
    file_name = FileHandler::swapExtension(file_name, supported_file_types.fromFileDialogFilter(selected_filter, fallback_extension)).toQString();
    return file_name;
  }



  bool GUIHelpers::startTOPPView(QStringList args)
  {
    QString app_path;
#if defined(__APPLE__)
    // check if we can find the TOPPView.app
    app_path = (File::getExecutablePath() + "../../../TOPPView.app").toQString();

    if (File::exists(app_path))
    {
      // we found the app
      QStringList app_args;
      app_args.append("-a");
      app_args.append(app_path);
      app_args.append("--args");
      app_args.append(args);
      args = app_args;
      app_path = "/usr/bin/open";
    }
    else
    { // we could not find the app, try it the Linux way
      app_path = (File::findSiblingTOPPExecutable("TOPPView")).toQString();
    }
#else
    // LINUX+WIN
    app_path = (File::findSiblingTOPPExecutable("TOPPView")).toQString();
#endif

    if (!QProcess::startDetached(app_path, args))
    {
      // execution failed
      OPENMS_LOG_ERROR << "Could not start '" << app_path.toStdString() << "'. Please see above for error messages." << std::endl;
  #if defined(__APPLE__)
      OPENMS_LOG_ERROR << "Please check if TOPPAS and TOPPView are located in the same directory" << std::endl;
  #endif
      return false;
    }
    return true;
  }

  void GUIHelpers::openURL(const QString& target)
  {
    QUrl url_target;

    // add protocol handler if none is given
    if (!(target.startsWith("http://") || target.startsWith("https://")))
    {
      // we expect all unqualified urls to be file urls
      try
      {
        String local_url = File::findDoc(target);
        url_target = QUrl::fromLocalFile(local_url.toQString());
      }
      catch (Exception::FileNotFound&)
      {
        // we fall back to the web url
        url_target = QUrl(QString("http://www.openms.de/current_doxygen/%1").arg(target), QUrl::TolerantMode);
      }
    }
    else
    {
      url_target = QUrl(target, QUrl::TolerantMode);
    }

    if (!QDesktopServices::openUrl(url_target))
    {
      QMessageBox::warning(nullptr,
                           QObject::tr("Error"),
                           QObject::tr("Unable to open\n") + target + 
                           QObject::tr("\n\nPossible reason: security settings or misconfigured Operating System"));
    }
  }

  void GUIHelpers::drawText(QPainter& painter, const QStringList& text, const QPoint& where, const QColor& col_fg, const QColor& col_bg, const QFont& f)
  {
    painter.save();

    // font
    painter.setFont(f);

    int line_spacing;
    QRectF dim = getTextDimension(text, painter.font(), line_spacing);

    // draw background for text
    if (col_bg.isValid()) painter.fillRect(where.x(), where.y(), dim.width(), dim.height(), col_bg);

    // draw text
    if (col_fg.isValid()) painter.setPen(col_fg);

    for (int i = 0; i < text.size(); ++i)
    {
      painter.drawText(where.x() + 1, where.y() + (i + 1) * line_spacing, text[i]);
    }
    painter.restore();
  }

  QRectF GUIHelpers::getTextDimension(const QStringList& text, const QFont& f, int& line_spacing)
  {
    // determine width and height of the box we need
    QFontMetrics metrics(f);
    line_spacing = metrics.lineSpacing();
    int height = 6 + text.size() * line_spacing;
    int width = 4;
    for (int i = 0; i < text.size(); ++i)
    {
      width = std::max(width, 4 + metrics.horizontalAdvance(text[i]));
    }
    return QRectF(0, 0, width, height);
  }

  QPointF GUIHelpers::nearestPoint(const QPointF& origin, const QList<QPointF>& list)
  {
    if (list.empty())
    {
      return QPointF();
    }
    QPointF nearest = list.first();
    qreal min_distance = (std::numeric_limits<qreal>::max)();

    for (QList<QPointF>::const_iterator it = list.begin(); it != list.end(); ++it)
    {
      qreal sqr_distance = (it->x() - origin.x()) * (it->x() - origin.x()) + (it->y() - origin.y()) * (it->y() - origin.y());
      if (sqr_distance < min_distance)
      {
        min_distance = sqr_distance;
        nearest = *it;
      }
    }

    return nearest;
  }

  QPointF GUIHelpers::intersectionPoint(const QRectF& rect, const QPointF& p)
  {
    if (rect.contains(p))
      return rect.center();

    QPointF delta = rect.center() - p;
    qreal slope;
    if (delta.x() == 0)
    {
      slope = std::numeric_limits<qreal>::infinity();
    }
    else
    {
      slope = delta.y() / delta.x();
    }

    qreal y_1 = p.y() + slope * (rect.left() - p.x());
    qreal y_2 = p.y() + slope * (rect.right() - p.x());

    slope = 1.0 / slope;

    qreal x_3 = p.x() + slope * (rect.top() - p.y());
    qreal x_4 = p.x() + slope * (rect.bottom() - p.y());

    QList<QPointF> point_list;
    if (y_1 <= rect.bottom() && y_1 >= rect.top())
    {
      point_list.push_back(QPointF(rect.left(), y_1));
    }
    if (y_2 <= rect.bottom() && y_2 >= rect.top())
    {
      point_list.push_back(QPointF(rect.right(), y_2));
    }
    if (x_3 <= rect.right() && x_3 >= rect.left())
    {
      point_list.push_back(QPointF(x_3, rect.top()));
    }
    if (x_4 <= rect.right() && x_4 >= rect.left())
    {
      point_list.push_back(QPointF(x_4, rect.bottom()));
    }

    return nearestPoint(p, point_list);
  }

  StringList GUIHelpers::convert(const QStringList& in)
  {
    StringList out;
    for (const auto& s : in) out.push_back(s);
    return out;
  }

  QStringList GUIHelpers::convert(const StringList& in)
  {
    QStringList out;
    for (const auto& s : in) out.push_back(s.toQString());
    return out;
  }

  GUIHelpers::GUILock::GUILock(QWidget* gui)
    : locked_widget_(gui)
  {
    lock();
  }

  GUIHelpers::GUILock::~GUILock()
  {
    unlock();
  }

  void GUIHelpers::GUILock::lock()
  {
    if (currently_locked_) return;
    if (locked_widget_ == nullptr) return;

    was_enabled_ = locked_widget_->isEnabled();
    locked_widget_->setEnabled(false);
    QGuiApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    currently_locked_ = true;
  }

  void GUIHelpers::GUILock::unlock()
  {
    if (!currently_locked_) return;
    if (locked_widget_ == nullptr) return;

    locked_widget_->setEnabled(was_enabled_);
    QGuiApplication::restoreOverrideCursor(); 
    currently_locked_ = false;
  }
  
  GUIHelpers::OverlapDetector::OverlapDetector(int levels)
  {
    if (levels <= 0)
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, levels);
    rows_.resize(levels, 0);
  }

  size_t GUIHelpers::OverlapDetector::placeItem(double x_start, double x_end)
  {
    if (x_start < 0)
      OPENMS_LOG_WARN << "Warning: x coordinates should be positive!\n";
    if (x_start > x_end)
      OPENMS_LOG_WARN << "Warning: x-end is larger than x-start!\n";

    size_t best_index = 0;
    double best_distance = std::numeric_limits<double>::max();
    for (size_t i = 0; i < rows_.size(); ++i)
    {
      if (rows_[i] < x_start)
      {                   // easy win; row[i] does not overlap; take it
        rows_[i] = x_end; // update space for next call
        return i;
      }
      // x_start is smaller than row's end...
      if ((rows_[i] - x_start) < best_distance)
      {
        best_distance = rows_[i] - x_start;
        best_index = i;
      }
    }

    rows_[best_index] = x_end; // update space for next call
    return best_index;
  }

} //namespace OpenMS

