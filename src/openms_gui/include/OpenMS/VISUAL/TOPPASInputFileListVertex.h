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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASINPUTFILELISTVERTEX_H
#define OPENMS_VISUAL_TOPPASINPUTFILELISTVERTEX_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
      @brief A vertex representing an input file list

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASInputFileListVertex :
    public TOPPASVertex
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASInputFileListVertex();
    /// Constructor
    TOPPASInputFileListVertex(const QStringList & files);
    /// Copy constructor
    TOPPASInputFileListVertex(const TOPPASInputFileListVertex & rhs);
    /// Destructor
    ~TOPPASInputFileListVertex() override;
    /// Assignment operator
    TOPPASInputFileListVertex & operator=(const TOPPASInputFileListVertex & rhs);
    /// returns "InputVertex"
    String getName() const override;
    /// Sets the list of files
    void setFilenames(const QStringList & files);
    /// Starts all tools below this node
    void run() override;
    // documented in base class
    void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget) override;
    // documented in base class
    QRectF boundingRect() const override;
    // documented in base class
    QPainterPath shape() const override;
    /// Checks if the given list of file names is valid
    bool fileNamesValid();
    /// Shows the dialog for editing the files
    void showFilesDialog();
    /// Opens the folders of the input files
    void openContainingFolder();
    /// Returns the key (for applying resources from a resource file)
    const QString & getKey();
    /// Sets the key (for applying resources from a resource file)
    void setKey(const QString & key);

public slots:
    /// Called by an outgoing edge when it has changed
    void outEdgeHasChanged() override;

protected:

    /// The key of this input node (for applying resources from a resource file)
    QString key_;

    /// current working dir, i.e. the last position a file was added from
    QString cwd_;

    ///@name reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent * e) override;
    //@}

  };
}

#endif
