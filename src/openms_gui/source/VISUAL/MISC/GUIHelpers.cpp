// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDesktopServices>
#include <QDir>
#include <QGuiApplication>
#include <QMessageBox>
#include <QPainter>
#include <QPoint>
#include <QProcess>
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

  OPENMS_GUI_DLLAPI QString GUIHelpers::getSaveFilename(QWidget* parent, const QString& caption, const QString& dir, FileTypeList supported_file_types, bool add_all_filter,
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

  void GUIHelpers::startTOPPView(const QStringList& args)
  {
    QProcess p;
    p.setProcessChannelMode(QProcess::ForwardedChannels);
#if defined(__APPLE__)
    // check if we can find the TOPPView.app
    QString app_path = (File::getExecutablePath() + "../../../TOPPView.app").toQString();

    if (File::exists(app_path))
    {
      // we found the app
      QStringList app_args;
      app_args.append("-a");
      app_args.append(app_path);
      app_args.append("--args");
      app_args.append(args);
      p.start("/usr/bin/open", app_args);
    }
    else
    {
      // we could not find the app, try it the Linux way
      QString toppview_executable = (File::findSiblingTOPPExecutable("TOPPView")).toQString();
      p.start(toppview_executable, args);
    }
#else
    // LINUX+WIN
    QString toppview_executable = (File::findSiblingTOPPExecutable("TOPPView")).toQString();
    p.start(toppview_executable, args);
#endif


    if (!p.waitForStarted())
    {
      // execution failed
      OPENMS_LOG_ERROR << p.errorString().toStdString() << std::endl;
  #if defined(Q_WS_MAC)
      OPENMS_LOG_ERROR << "Please check if TOPPAS and TOPPView are located in the same directory" << std::endl;
  #endif
    }
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

  void GUIHelpers::drawText(QPainter& painter, const QStringList& text, const QPoint& where, const QColor col_fg, const QColor col_bg, const QFont& f)
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
      /*
       * QFontMetrics::width() is deprecated after Qt 5.11. Use QFontMetrics::horizontalAdvance()
       */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      width = std::max(width, 4 + metrics.width(text[i]));
#pragma GCC diagnostic pop
    }
    return QRectF(0, 0, width, height);
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

