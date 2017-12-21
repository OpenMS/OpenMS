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

#ifndef OPENMS_VISUAL_TOPPASMERGERVERTEX_H
#define OPENMS_VISUAL_TOPPASMERGERVERTEX_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
      @brief A special vertex that allows to merge several inputs.

      A special vertex that allows to merge several inputs. Mergers have two modes: The normal,
      round-based merging mode and a "wait & merge all" mode. In round-based mode, a merger
      first takes the first files of each incoming file list and merges them into a list (which
      has as many elements as the merger has incoming edges).

      In "wait & merge all" mode, the merger first waits for all upstream mergers to finish all
      their merging rounds and then merges all collected files from all merging rounds for all
      incoming edges into one single list and calls the next tool with this list of files as input.

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASMergerVertex :
    public TOPPASVertex
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASMergerVertex();
    /// Constructor
    TOPPASMergerVertex(bool round_based);
    /// Copy constructor
    TOPPASMergerVertex(const TOPPASMergerVertex& rhs);
    /// Destructor
    ~TOPPASMergerVertex() override;
    /// Assignment operator
    TOPPASMergerVertex& operator=(const TOPPASMergerVertex& rhs);
    /// returns "MergerVertex"
    String getName() const override;
    /// check if upstream nodes are finished and call downstream nodes
    void run() override;
    /// Determines whether this merger is merging round based or merging all inputs into one list
    bool roundBasedMode();
    // documented in base class
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) override;
    // documented in base class
    QRectF boundingRect() const override;
    // documented in base class
    QPainterPath shape() const override;
    // documented in base class
    void markUnreachable() override;

public slots:

signals:
    /// Emitted when merging upstream data failed
    void mergeFailed(const QString message);

protected:

    /// Stores whether this merger is merging round based or merging all inputs into one list
    bool round_based_mode_;

    ///@name reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e) override;
    //@}


  };
}

#endif
