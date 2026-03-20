// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASOutputVertex.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <QtCore/QDir>

namespace OpenMS
{
  TOPPASOutputVertex::TOPPASOutputVertex(const TOPPASOutputVertex& rhs):
      TOPPASVertex(rhs),
      output_folder_name_() // leave empty on copy, otherwise we will have conficting output folder names
  {
  }

  TOPPASOutputVertex& TOPPASOutputVertex::operator=(const TOPPASOutputVertex& rhs)
  {
    TOPPASVertex::operator=(rhs);
    output_folder_name_ = ""; // leave empty on copy, otherwise we will have conficting output folder names

    return *this;
  }

  void TOPPASOutputVertex::setOutputFolderName(const QString& name)
  {
    if (output_folder_name_ != name)
    {
      output_folder_name_ = name;
      emit outputFolderNameChanged(); // enable storing of modified pipeline
    }
  }

  const QString& TOPPASOutputVertex::getOutputFolderName() const
  {
    return output_folder_name_;
  }

  void TOPPASOutputVertex::inEdgeHasChanged()
  {
    reset(true);
    qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
    TOPPASVertex::inEdgeHasChanged();
  }

  void TOPPASOutputVertex::openContainingFolder() const
  {
    GUIHelpers::openFolder(getFullOutputDirectory().toQString());
  }

  String TOPPASOutputVertex::getFullOutputDirectory() const
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    auto dir = String(ts->getOutDir()).substitute("\\", "/").ensureLastChar('/') + getOutputDir();
    String clean_dir = QDir::cleanPath(dir.toQString());
    return clean_dir.substitute("\\", "/").ensureLastChar('/');
  }

  String TOPPASOutputVertex::getOutputDir() const
  {
    String dir = "TOPPAS_out/";
    if (output_folder_name_.isEmpty())
    {
      TOPPASEdge* e = *inEdgesBegin();
      if (e == nullptr)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "To open the output folder, an input edge is required to knit a folder name.");
      }
      const TOPPASVertex* tv = e->getSourceVertex();
      // create meaningful output name using vertex + TOPP name + output parameter, e.g. "010-FileConverter-out"
      dir += get3CharsNumber_(topo_nr_) + "-" + tv->getName() + "-" + e->getSourceOutParamName().remove(':');
    }
    else { 
      dir += output_folder_name_;
    }
    dir.ensureLastChar('/');
    return dir;
  }

  String TOPPASOutputVertex::createOutputDir() const
  {
    String full_dir = getFullOutputDirectory();
    if (! File::exists(full_dir))
    {
      if (! File::makeDir(full_dir)) { std::cerr << "Could not create path " << full_dir << std::endl; }
    }
    return full_dir;
  }

  void TOPPASOutputVertex::setTopoNr(UInt nr)
  {
    if (topo_nr_ != nr)
    {
      // topological number changes --> output dir changes --> reset
      reset(true);
      topo_nr_ = nr;
    }
  }

  void TOPPASOutputVertex::reset(bool reset_all_files)
  {
    __DEBUG_BEGIN_METHOD__

    files_total_ = 0;
    files_written_ = 0;

    TOPPASVertex::reset();

    if (reset_all_files)
    {
      // do not actually delete the output files here
      // for TOPPASOutputFolderVertex we do not even know which files were written by the tool
      // and which were already there (we do not want to delete user files)
    }

    __DEBUG_END_METHOD__
  }

  void TOPPASOutputVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
  {
    openContainingFolder();
  }
} // namespace OpenMS
