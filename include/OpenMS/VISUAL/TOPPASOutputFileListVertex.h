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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASOUTPUTFILELISTVERTEX_H
#define OPENMS_VISUAL_TOPPASOUTPUTFILELISTVERTEX_H

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
      @brief A vertex representing an output file list

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASOutputFileListVertex :
    public TOPPASVertex
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASOutputFileListVertex();
    /// Copy constructor
    TOPPASOutputFileListVertex(const TOPPASOutputFileListVertex & rhs);
    /// Destructor
    virtual ~TOPPASOutputFileListVertex();
    /// Assignment operator
    TOPPASOutputFileListVertex & operator=(const TOPPASOutputFileListVertex & rhs);
    /// returns "OutputVertex"
    virtual String getName() const;
    // documented in base class
    virtual void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget);
    // documented in base class
    virtual QRectF boundingRect() const;
    // documented in base class
    virtual QPainterPath shape() const;
    // documented in base class
    virtual void reset(bool reset_all_files = false);
    /// Called when the parent node has finished execution
    virtual void run();
    /// Returns the full directory (including preceding output path as selected by user)
    String getFullOutputDirectory() const;
    /// Returns the directory where the output files are stored
    String getOutputDir() const;
    /// Creates the output directory for this node
    String createOutputDir();
    /// Sets the topological sort number and removes invalidated tmp files
    virtual void setTopoNr(UInt nr);
    /// Opens the folders of the output files
    void openContainingFolder();

public slots:

    //documented in base class
    virtual void inEdgeHasChanged();

signals:

    /// Emitted when an output file was written
    void outputFileWritten(const String & file);
    /// Emitted when the pipeline ending in this vertex is finished
    void iAmDone();

protected:
    static bool copy_(const QString & from, const QString & to); //< STATIC(!) function which calls QFile::copy(); needs to be static, since we need to pass a function pointer (which does not work on member functions)
    // convenience members, not required for operation, but for progress during copying
    int files_written_;       //< files that were already written
    int files_total_;     //< total number of files from upstream
  };
}

#endif
