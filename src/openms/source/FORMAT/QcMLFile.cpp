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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/copy.hpp>

using namespace std;

// TODO fix all the shadowed "const_iterator qpsit"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

namespace OpenMS
{

  QcMLFile::QualityParameter::QualityParameter() :
    name(),
    id(),
    value(),
    cvRef(),
    cvAcc(),
    unitRef(),
    unitAcc(),
    flag()
  {
  }

  QcMLFile::QualityParameter::QualityParameter(const QualityParameter& rhs) :
    name(rhs.name),
    id(rhs.id),
    value(rhs.value),
    cvRef(rhs.cvRef),
    cvAcc(rhs.cvAcc),
    unitRef(rhs.unitRef),
    unitAcc(rhs.unitAcc),
    flag(rhs.flag)
  {
  }

  QcMLFile::QualityParameter& QcMLFile::QualityParameter::operator=(const QualityParameter& rhs)
  {
    if (this != &rhs)
    {
      name = rhs.name;
      id = rhs.id;
      value = rhs.value;
      cvRef = rhs.cvRef;
      cvAcc = rhs.cvAcc;
      unitRef = rhs.unitRef;
      unitAcc = rhs.unitAcc;
      flag = rhs.flag;
    }
    return *this;
  }

  bool QcMLFile::QualityParameter::operator<(const QualityParameter& rhs) const
  {
    return name.toQString() < rhs.name.toQString();
  }

  bool QcMLFile::QualityParameter::operator>(const QualityParameter& rhs) const
  {
    return name.toQString() > rhs.name.toQString();
  }

  bool QcMLFile::QualityParameter::operator==(const QualityParameter& rhs) const
  {
    return name.toQString() == rhs.name.toQString();
  }

  String QcMLFile::QualityParameter::toXMLString(UInt indentation_level) const
  {
    String indent = String(indentation_level, '\t');
    String s = indent;
    s += "<qualityParameter";
    s += " name=\"" + name + "\"" + " ID=\"" + id + "\"" + " cvRef=\"" + cvRef + "\"" + " accession=\"" + cvAcc + "\"";
    if (value != "")
    {
      s += " value=\"" + value + "\"";
    }
    if (unitRef != "")
    {
      s += " unitRef=\"" + unitRef + "\"";
    }
    if (unitAcc != "")
    {
      s += " unitAcc=\"" + unitAcc + "\"";
    }
    if (!flag.empty())
    {
      s += " flag=\"true\"";
    }

    s += "/>\n";

    return s;
  }

  QcMLFile::Attachment::Attachment() :
    name(),
    id(),
    value(),
    cvRef(),
    cvAcc(),
    unitRef(),
    unitAcc(),
    qualityRef(),
    colTypes(),
    tableRows()
  {
  }

  QcMLFile::Attachment::Attachment(const Attachment& rhs) :
    name(rhs.name),
    id(rhs.id),
    value(rhs.value),
    cvRef(rhs.cvRef),
    cvAcc(rhs.cvAcc),
    unitRef(rhs.unitRef),
    unitAcc(rhs.unitAcc),
    binary(rhs.binary),
    qualityRef(rhs.qualityRef),
    colTypes(rhs.colTypes),
    tableRows(rhs.tableRows)
  {
  }

  QcMLFile::Attachment& QcMLFile::Attachment::operator=(const Attachment& rhs)
  {
    if (this != &rhs)
    {
      name = rhs.name;
      id = rhs.id;
      value = rhs.value;
      cvRef = rhs.cvRef;
      cvAcc = rhs.cvAcc;
      unitRef = rhs.unitRef;
      unitAcc = rhs.unitAcc;
      binary = rhs.binary;
      qualityRef = rhs.qualityRef;
      colTypes = rhs.colTypes;
      tableRows = rhs.tableRows;
    }
    return *this;
  }

  bool QcMLFile::Attachment::operator<(const Attachment& rhs) const
  {
    return name.toQString() < rhs.name.toQString();
  }

  bool QcMLFile::Attachment::operator>(const Attachment& rhs) const
  {
    return name.toQString() > rhs.name.toQString();
  }

  bool QcMLFile::Attachment::operator==(const Attachment& rhs) const
  {
    return name.toQString() == rhs.name.toQString();
  }

  String QcMLFile::Attachment::toCSVString(const String separator) const
  {
    String s = "";
    if ((!colTypes.empty()) && (!tableRows.empty()))
    {
      String replacement = "_";
      if (separator == replacement)
      {
        replacement = "$";
      }
      std::vector<String> copy = colTypes;
      for (std::vector<String>::iterator it = copy.begin(); it != copy.end(); ++it)
      {
        it->substitute(separator, replacement);
      }
      s += ListUtils::concatenate(copy, separator).trim();
      s += "\n";
      for (std::vector<std::vector<String> >::const_iterator it = tableRows.begin(); it != tableRows.end(); ++it)
      {
        std::vector<String> copy_row = *it;
        for (std::vector<String>::iterator sit = copy_row.begin(); sit != copy_row.end(); ++sit)
        {
          sit->substitute(separator, replacement);
        }
        s += ListUtils::concatenate(copy_row, separator).trim();
        s += "\n";
      }
    }
    return s;
  }

  //TODO CHANGE TO ACCEPTABLE WAY OF GENERATING A PRETTY XML FILE
  String QcMLFile::Attachment::toXMLString(UInt indentation_level) const
  {
    //TODO manage IDREF to qp internally
    String indent = String(indentation_level, '\t');
    String s = indent;
    s += "<attachment ";
    s += " name=\"" + name + "\"" + " ID=\"" + id + "\"" + " cvRef=\"" + cvRef + "\"" + " accession=\"" + cvAcc + "\"";
    if (value != "")
    {
      s += " value=\"" + value + "\"";
    }
    if (unitRef != "")
    {
      s += " unitRef=\"" + unitRef + "\"";
    }
    if (unitAcc != "")
    {
      s += " unitAcc=\"" + unitAcc + "\"";
    }
    if (qualityRef != "")
    {
      s += " qualityParameterRef=\"" + qualityRef + "\"";
    }

    if (binary != "")
    {
      s += ">\n";
      s += indent + "\t" + "<binary>" + binary + "</binary>\n";
      s += indent + "</attachment>\n";
    }
    else if ((!colTypes.empty()) && (!tableRows.empty()))
    {
      s += ">\n";
      s += "<table>";
      s += indent + "\t" + "<tableColumnTypes>";

      std::vector<String> copy = colTypes;
      for (std::vector<String>::iterator it = copy.begin(); it != copy.end(); ++it)
      {
        it->substitute(String(" "), String("_"));
      }

      s += ListUtils::concatenate(copy, " ").trim();
      s += "</tableColumnTypes>\n";
      for (std::vector<std::vector<String> >::const_iterator it = tableRows.begin(); it != tableRows.end(); ++it)
      {
        s += indent + "\t" + "<tableRowValues>";

        std::vector<String> copy_row = *it;
        for (std::vector<String>::iterator sit = copy_row.begin(); sit != copy_row.end(); ++sit)
        {
          sit->substitute(String(" "), String("_"));
        }

        s += ListUtils::concatenate(*it, " ").trim();
        s += "</tableRowValues>\n";
      }
      s += "</table>";
      s += indent + "</attachment>\n";
    }
    else
    {
      //TODO warning invalid attachment!
      return "";
    }

    return s;
  }

  QcMLFile::QcMLFile() :
    XMLHandler("", "0.7"), XMLFile("/SCHEMAS/qcml.xsd", "0.7"), ProgressLogger() //TODO keep version uptodate
  {
  }

  QcMLFile::~QcMLFile()
  {

  }

  void QcMLFile::addRunQualityParameter(String run_id, QualityParameter qp)
  {
    // TODO warn that run has to be registered!
    std::map<String, std::vector<QcMLFile::QualityParameter> >::const_iterator qpsit = runQualityQPs_.find(run_id); //if 'filename is a ID:'
    if (qpsit != runQualityQPs_.end())
    {
      runQualityQPs_[run_id].push_back(qp);
    }
    else
    {
      std::map<String, String>::const_iterator qpsit = run_Name_ID_map_.find(run_id); //if 'filename' is a name
      if (qpsit != run_Name_ID_map_.end())
      {
        runQualityQPs_[qpsit->second].push_back(qp);
      }
    }
    //TODO redundancy check
  }

  void QcMLFile::addSetQualityParameter(String set_id, QualityParameter qp)
  {
    // TODO warn that set has to be registered!
    std::map<String, std::vector<QcMLFile::QualityParameter> >::const_iterator qpsit = setQualityQPs_.find(set_id); //if 'filename is a ID:'
    if (qpsit != setQualityQPs_.end())
    {
      setQualityQPs_[set_id].push_back(qp);
    }
    else
    {
      std::map<String, String>::const_iterator qpsit = set_Name_ID_map_.find(set_id); //if 'filename' is a name
      if (qpsit != set_Name_ID_map_.end())
      {
        setQualityQPs_[qpsit->second].push_back(qp);
      }
    }
    //TODO redundancy check
  }

  void QcMLFile::addRunAttachment(String run_id, Attachment at)
  {
    runQualityAts_[run_id].push_back(at); //TODO permit AT without a QP (or enable orphan writeout in store),redundancy check
  }

  void QcMLFile::addSetAttachment(String run_id, Attachment at)
  {
    setQualityAts_[run_id].push_back(at); //TODO add file QP to set member
  }

  void QcMLFile::getRunNames(std::vector<String>& ids) const
  {
    ids.clear();
    boost::copy(run_Name_ID_map_ | boost::adaptors::map_keys, std::back_inserter(ids));
  }

  void QcMLFile::getRunIDs(std::vector<String>& ids) const
  {
    ids.clear();
    boost::copy(runQualityQPs_ | boost::adaptors::map_keys, std::back_inserter(ids));
  }

  bool QcMLFile::existsRun(const String filename, bool checkname) const
  {
    std::map<String, std::vector<QcMLFile::QualityParameter> >::const_iterator qpsit = runQualityQPs_.find(filename); //if 'filename is a ID:'
    if (qpsit != runQualityQPs_.end()) //NO, do not!: permit AT without a QP
    {
      return true;
    }
    else if (checkname)
    {
      std::map<String, String>::const_iterator qpsit = run_Name_ID_map_.find(filename); //if 'filename' is a name
      if (qpsit != run_Name_ID_map_.end()) //NO, do not!: permit AT without a QP
      {
        return true;
      }

    }
    return false;
  }

  bool QcMLFile::existsSet(const String filename, bool checkname) const
  {
    std::map<String, std::vector<QcMLFile::QualityParameter> >::const_iterator qpsit = setQualityQPs_.find(filename); //if 'filename is a ID:'
    if (qpsit != setQualityQPs_.end()) //NO, do not!: permit AT without a QP
    {
      return true;
    }
    else if (checkname)
    {
      std::map<String, String>::const_iterator qpsit = set_Name_ID_map_.find(filename); //if 'filename' is a name
      if (qpsit != set_Name_ID_map_.end()) //NO, do not!: permit AT without a QP
      {
        return true;
      }
    }
    return false;
  }

  void QcMLFile::existsRunQualityParameter(const String filename, const String qpname, std::vector<String>& ids) const
  {
    ids.clear();
    std::map<String, std::vector<QcMLFile::QualityParameter> >::const_iterator qpsit = runQualityQPs_.find(filename);
    if (qpsit == runQualityQPs_.end()) //try name mapping if 'filename' is no ID but name
    {
      std::map<String, String>::const_iterator mapsit = run_Name_ID_map_.find(filename);
      if (mapsit != run_Name_ID_map_.end())
      {
        qpsit = runQualityQPs_.find(mapsit->second);
      }
    }
    if (qpsit != runQualityQPs_.end())
    {
      for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        if (qpname == qit->cvAcc)
        {
          ids.push_back(qit->id);
        }
      }
    }
  }

  void QcMLFile::existsSetQualityParameter(const String filename, const String qpname, std::vector<String>& ids) const
  {
    ids.clear();
    std::map<String, std::vector<QcMLFile::QualityParameter> >::const_iterator qpsit = setQualityQPs_.find(filename);
    if (qpsit == setQualityQPs_.end()) //try name mapping if 'filename' is no ID but name
    {
      std::map<String, String>::const_iterator mapsit = set_Name_ID_map_.find(filename);
      if (mapsit != set_Name_ID_map_.end())
      {
        qpsit = setQualityQPs_.find(mapsit->second);
      }
    }
    if (qpsit != setQualityQPs_.end())
    {
      for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        //~ std::cout << qit->name << "setexists" << std::endl;
        //~ std::cout << qpname << "qpname" << std::endl;
        if (qpname == qit->cvAcc)
        {
          ids.push_back(qit->id);
        }
      }
    }
  }

  void QcMLFile::removeQualityParameter(String r, std::vector<String>& ids)
  {
    removeAttachment(r, ids);
    for (Size i = 0; i < ids.size(); ++i)
    {
      std::vector<QcMLFile::QualityParameter>::iterator qit = runQualityQPs_[r].begin();
      while (qit != runQualityQPs_[r].end())
      {
        if (qit->id == ids[i])
        {
          qit = runQualityQPs_[r].erase(qit);
        }
        else
        {
          ++qit;
        }
      }
      qit = setQualityQPs_[r].begin();
      while (qit != setQualityQPs_[r].end())
      {
        if (qit->id == ids[i])
        {
          qit = setQualityQPs_[r].erase(qit);
        }
        else
        {
          ++qit;
        }
      }
    }
  }

  void QcMLFile::removeAttachment(String r, std::vector<String>& ids, String at)
  {
    bool not_all = at.size();
    for (Size i = 0; i < ids.size(); ++i)
    {
      std::vector<QcMLFile::Attachment>::iterator qit = runQualityAts_[r].begin();
      while (qit != runQualityAts_[r].end())
      {
        if (qit->qualityRef == ids[i] && ((qit->name == at) || (!not_all)))
        {
          qit = runQualityAts_[r].erase(qit);
        }
        else
        {
          ++qit;
        }
      }
      qit = setQualityAts_[r].begin();
      while (qit != setQualityAts_[r].end())
      {
        if (qit->qualityRef == ids[i] && ((qit->name == at) || (!not_all)))
        {
          qit = setQualityAts_[r].erase(qit);
        }
        else
        {
          ++qit;
        }
      }
    }
  }

  void QcMLFile::removeAttachment(String r, String at)
  {
    if (existsRun(r))
    {
      std::vector<QcMLFile::Attachment>::iterator qit = runQualityAts_[r].begin();
      //~ cout << "remove from " << r << endl;
      while (qit != runQualityAts_[r].end())
      {
        if (qit->cvAcc == at)
        {
          qit = runQualityAts_[r].erase(qit);
          //~ cout << "remove  " << at << endl;
        }
        else
        {
          ++qit;
        }
      }
    }
    if (existsSet(r))
    {
      std::vector<QcMLFile::Attachment>::iterator qit = setQualityAts_[r].begin();
      while (qit != setQualityAts_[r].end())
      {
        if (qit->cvAcc == at)
        {
          qit = setQualityAts_[r].erase(qit);
        }
        else
        {
          ++qit;
        }
      }
    }
  }

  void QcMLFile::removeAllAttachments(String at)
  {
    for (std::map<String, std::vector<Attachment> >::iterator it = runQualityAts_.begin(); it != runQualityAts_.end(); ++it)
    {
      removeAttachment(it->first, at);
    }
  }

  void QcMLFile::registerRun(const String id, const String name)
  {
    runQualityQPs_[id] = std::vector<QualityParameter>();
    runQualityAts_[id] = std::vector<Attachment>();
    run_Name_ID_map_[name] = id;
  }

  void QcMLFile::registerSet(const String id, const String name, const std::set<String>& names)
  {
    setQualityQPs_[id] = std::vector<QualityParameter>();
    setQualityAts_[id] = std::vector<Attachment>();
    set_Name_ID_map_[name] = id;
    setQualityQPs_members_[id] = names;
  }

  void QcMLFile::merge(const QcMLFile& addendum, String setname)
  {
    //~ runs (and create set if setname not empty)
    for (std::map<String, std::vector<QualityParameter> >::const_iterator it = addendum.runQualityQPs_.begin(); it != addendum.runQualityQPs_.end(); ++it)
    {
      runQualityQPs_[it->first].insert(runQualityQPs_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(runQualityQPs_[it->first].begin(), runQualityQPs_[it->first].end());
      runQualityQPs_[it->first].erase(std::unique(runQualityQPs_[it->first].begin(), runQualityQPs_[it->first].end()), runQualityQPs_[it->first].end());
      if (setname != "")
      {
        setQualityQPs_members_[setname].insert(it->first);
      }
    }
    for (std::map<String, std::vector<Attachment> >::const_iterator it = addendum.runQualityAts_.begin(); it != addendum.runQualityAts_.end(); ++it)
    {
      runQualityAts_[it->first].insert(runQualityAts_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(runQualityAts_[it->first].begin(), runQualityAts_[it->first].end());
      runQualityAts_[it->first].erase(std::unique(runQualityAts_[it->first].begin(), runQualityAts_[it->first].end()), runQualityAts_[it->first].end());
      if (setname != "")
      {
        setQualityQPs_members_[setname].insert(it->first);
      }
    }

    // sets
    //~ TODO sets are not supposed to overlap - throw error if so
    setQualityQPs_members_.insert(addendum.setQualityQPs_members_.begin(), addendum.setQualityQPs_members_.end());
    for (std::map<String, std::vector<QualityParameter> >::const_iterator it = addendum.setQualityQPs_.begin(); it != addendum.setQualityQPs_.end(); ++it)
    {
      setQualityQPs_[it->first].insert(setQualityQPs_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(setQualityQPs_[it->first].begin(), setQualityQPs_[it->first].end());
      setQualityQPs_[it->first].erase(std::unique(setQualityQPs_[it->first].begin(), setQualityQPs_[it->first].end()), setQualityQPs_[it->first].end());
    }
    for (std::map<String, std::vector<Attachment> >::const_iterator it = addendum.setQualityAts_.begin(); it != addendum.setQualityAts_.end(); ++it)
    {
      setQualityAts_[it->first].insert(setQualityAts_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(setQualityAts_[it->first].begin(), setQualityAts_[it->first].end());
      setQualityAts_[it->first].erase(std::unique(setQualityAts_[it->first].begin(), setQualityAts_[it->first].end()), setQualityAts_[it->first].end());
    }
  }

  String QcMLFile::exportAttachment(const String filename, const String qpname) const
  {
    std::map<String, std::vector<QcMLFile::Attachment> >::const_iterator qpsit = runQualityAts_.find(filename);
    if (qpsit == runQualityAts_.end()) //try name mapping if 'filename' is no ID but name
    {
      std::map<String, String>::const_iterator mapsit = run_Name_ID_map_.find(filename);
      if (mapsit != run_Name_ID_map_.end())
      {
        qpsit = runQualityAts_.find(mapsit->second);
      }
    }
    if (qpsit != runQualityAts_.end())
    {
      for (std::vector<QcMLFile::Attachment>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        if ((qpname == qit->name) || (qpname == qit->cvAcc))
        {
          return qit->toCSVString("\t");
          //~ return qit->toXMLString(1);
        }
      }
    }

    // if the return statement wasn't hit from runs maybe it is from sets?
    qpsit = setQualityAts_.find(filename);
    if (qpsit == setQualityAts_.end()) //try name mapping if 'filename' is no ID but name
    {
      std::map<String, String>::const_iterator mapsit = set_Name_ID_map_.find(filename);
      if (mapsit != set_Name_ID_map_.end())
      {
        qpsit = setQualityAts_.find(mapsit->second);
      }
    }
    if (qpsit != setQualityAts_.end())
    {
      for (std::vector<QcMLFile::Attachment>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        if ((qpname == qit->name) || (qpname == qit->cvAcc))
        {
          return qit->toCSVString("\t");
          //~ return qit->toXMLString(1);
        }
      }
    }

    return "";
  }

  String QcMLFile::exportQP(const String filename, const String qpname) const
  {
    std::map<String, std::vector<QcMLFile::QualityParameter> >::const_iterator qpsit = runQualityQPs_.find(filename);
    if (qpsit == runQualityQPs_.end()) //try name mapping if 'filename' is no ID but name
    {
      std::map<String, String>::const_iterator mapsit = run_Name_ID_map_.find(filename);
      if (mapsit != run_Name_ID_map_.end())
      {
        qpsit = runQualityQPs_.find(mapsit->second);
      }
    }
    if (qpsit != runQualityQPs_.end())
    {
      for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        if (qpname == qit->cvAcc)
        {
          return /* "\""+ */ qit->value /* +"\"" */;
        }
      }
    }

    // if the return statement wasn't hit from runs maybe it is from sets?
    qpsit = setQualityQPs_.find(filename);
    if (qpsit == setQualityQPs_.end()) //try name mapping if 'filename' is no ID but name
    {
      std::map<String, String>::const_iterator mapsit = set_Name_ID_map_.find(filename);
      if (mapsit != set_Name_ID_map_.end())
      {
        qpsit = setQualityQPs_.find(mapsit->second);
      }
    }
    if (qpsit != setQualityQPs_.end())
    {
      for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        if (qpname == qit->name)
        {
          return /* "\""+ */ qit->value /* +"\"" */;
        }
      }
    }

    return "N/A";
  }

  String QcMLFile::exportQPs(const String filename, const StringList qpnames) const
  {
    String ret = "";
    for (StringList::const_iterator qit = qpnames.begin(); qit != qpnames.end(); ++qit)
    {
      ret += exportQP(filename, *qit);
      ret += ",";
    }
    return ret;
  }

  String QcMLFile::map2csv(const std::map<String, std::map<String, String> >& cvs_table, const String& separator) const
  {
    String ret = "";
    std::vector<String> cols;
    if (!cvs_table.empty())
    {
      for (std::map<String, String>::const_iterator it = cvs_table.begin()->second.begin(); it != cvs_table.begin()->second.end(); ++it)
      {
        cols.push_back(it->first);
      }
      ret += "qp";
      ret += separator;
      for (std::vector<String>::const_iterator jt = cols.begin(); jt != cols.end(); ++jt)
      {
        ret += *jt;
        ret += separator;
      }
      ret += "\n";
      for (std::map<String, std::map<String, String> >::const_iterator it = cvs_table.begin(); it != cvs_table.end(); ++it)
      {
        ret += it->first;
        ret += separator;
        for (std::vector<String>::const_iterator jt = cols.begin(); jt != cols.end(); ++jt)
        {
          std::map<String, String>::const_iterator found = it->second.find(*jt);
          if (found != it->second.end())
          {
            ret += found->second;
            ret += separator;
          } //TODO else throw error
        }
        ret += "\n";
      }
    }
    return ret;
  }

  String QcMLFile::exportIDstats(const String& filename) const
  {
    std::map<String, std::vector<QualityParameter> >::const_iterator found = setQualityQPs_.find(filename);
    if (found == setQualityQPs_.end()) //try name mapping if 'filename' is no ID but name
    {
      std::map<String, String>::const_iterator mapsit = set_Name_ID_map_.find(filename);
      if (mapsit != set_Name_ID_map_.end())
      {
        found = setQualityQPs_.find(mapsit->second);
      }
    }
    if (found != setQualityQPs_.end())
    {
      std::map<String, std::map<String, String> > cvs_table;
      for (std::vector<QualityParameter>::const_iterator it = found->second.begin(); it != found->second.end(); ++it)
      {
        if (it->cvAcc == "QC:0000043" || it->cvAcc == "QC:0000044" || it->cvAcc == "QC:0000045" || it->cvAcc == "QC:0000046" || it->cvAcc == "QC:0000047")
        {
          cvs_table["id"][it->name.prefix(' ')] = it->value;
        }
        else if (it->cvAcc == "QC:0000053" || it->cvAcc == "QC:0000054" || it->cvAcc == "QC:0000055" || it->cvAcc == "QC:0000056" || it->cvAcc == "QC:0000057")
        {
          cvs_table["ms2"][it->name.prefix(' ')] = it->value;
        }
      }
      if (!cvs_table.empty())
      {
        return map2csv(cvs_table, "\t");
      }
    }

    return "";
  }

  void QcMLFile::collectSetParameter(const String setname, const String qp, std::vector<String>& ret)
  {
    for (std::set<String>::const_iterator it = setQualityQPs_members_[setname].begin(); it != setQualityQPs_members_[setname].end(); ++it)
    {
      for (std::vector<QualityParameter>::const_iterator jt = runQualityQPs_[*it].begin(); jt != runQualityQPs_[*it].end(); ++jt)
      {
        if (jt->cvAcc == qp)
        {
          ret.push_back(jt->value);
        }
      }
    }
  }

  void QcMLFile::load(const String& filename)
  {
    //Filename for error messages in XMLHandler
    file_ = filename;

    runQualityQPs_.clear(); // clear
    runQualityAts_.clear(); // clear
    setQualityQPs_.clear(); // clear
    setQualityAts_.clear(); // clear
    setQualityQPs_members_.clear(); // clear

    parse_(filename, this);
  }

  void QcMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    tag_ = sm_.convert(qname);
    String parent_tag;
    if (!open_tags_.empty())
    {
      parent_tag = open_tags_.back();
    }
    open_tags_.push_back(tag_);

    static set<String> to_ignore;
    if (to_ignore.empty())
    {
      to_ignore.insert("tableColumnTypes"); // will be handled entirely in characters.
      to_ignore.insert("tableRowValues"); // ...
      to_ignore.insert("binary"); // ...
    }

    if (to_ignore.find(tag_) != to_ignore.end())
    {
      return;
    }

    String tmp_str;
    if (tag_ == "qcML")
    {
      startProgress(0, 0, "loading qcML file");
      progress_ = 0;
      setProgress(++progress_);
    }
    else if (tag_ == "runQuality")
    {
      run_id_ = attributeAsString_(attributes, "ID"); //TODO!
      setProgress(++progress_);
      qps_.clear();
      ats_.clear();
      qp_ = QualityParameter();
      at_ = Attachment();
      name_ = "";
      //for the run name wait for the qp with the right cv, otherwise use a uid
    }
    else if (tag_ == "qualityParameter")
    {
      optionalAttributeAsString_(qp_.value, attributes, "value");
      optionalAttributeAsString_(qp_.unitAcc, attributes, "unitAccession");
      optionalAttributeAsString_(qp_.unitRef, attributes, "unitCvRef");
      optionalAttributeAsString_(qp_.flag, attributes, "flag");
      qp_.cvRef = attributeAsString_(attributes, "cvRef");
      qp_.cvAcc = attributeAsString_(attributes, "accession");
      qp_.id = attributeAsString_(attributes, "ID");
      qp_.name = attributeAsString_(attributes, "name");
      if (parent_tag == "runQuality")
      {
        if (qp_.cvAcc == "MS:1000577") //no own qc cv
        {
          name_ = qp_.value;
        }
        //TODO add cvhandling for validation etc
      }
      else //setQuality
      {
        if (qp_.cvAcc == "MS:1000577") //TODO make sure these exist in runs later!
        {
          names_.insert(qp_.value);
        }
        if (qp_.cvAcc == "QC:0000058") //id: MS:1000577 name: raw data file  - with value of the file name of the run
        {
          name_ = qp_.value;
        }
      }
    }
    else if (tag_ == "attachment")
    {
      optionalAttributeAsString_(at_.value, attributes, "value");
      optionalAttributeAsString_(at_.unitAcc, attributes, "unitAccession");
      optionalAttributeAsString_(at_.unitRef, attributes, "unitCvRef");
      at_.cvRef = attributeAsString_(attributes, "cvRef");
      at_.cvAcc = attributeAsString_(attributes, "accession");
      at_.name = attributeAsString_(attributes, "name");
      at_.id = attributeAsString_(attributes, "ID");
      at_.qualityRef = attributeAsString_(attributes, "qualityParameterRef");
    }
    else if (tag_ == "setQuality")
    {
      setProgress(++progress_);
      run_id_ = attributeAsString_(attributes, "ID"); //TODO!
      qps_.clear();
      ats_.clear();
      qp_ = QualityParameter();
      at_ = Attachment();
      name_ = "";
    }
  }

  void QcMLFile::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
    if (tag_ == "tableRowValues")
    {
      String s = sm_.convert(chars);
      s.trim();
      if (!s.empty()) // always two notifications for a row, only the first one contains chars - dunno why
      {
        s.split(" ", row_);
      }
    }
    else if (tag_ == "tableColumnTypes")
    {
      String s = sm_.convert(chars);
      if (!s.empty()) // always two notifications for a row, only the first one contains chars - dunno why
      {
        s.split(" ", header_);
      }
    }
    else if (tag_ == "binary")
    {
      //chars may be split to several chunks => concatenate them
      at_.binary += sm_.convert(chars);
      //~ at_.binary = "bla";
    }
  }

  void QcMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    static set<String> to_ignore;
    if (to_ignore.empty())
    {
      //~ to_ignore.insert("binary");
    }

    tag_ = sm_.convert(qname);

    //determine parent tag
    String parent_tag;
    if (open_tags_.size() > 1)
    {
      parent_tag = *(open_tags_.end() - 2);
    }
    String parent_parent_tag;
    if (open_tags_.size() > 2)
    {
      parent_parent_tag = *(open_tags_.end() - 3);
    }

    //close current tag
    open_tags_.pop_back();

    if (to_ignore.find(tag_) != to_ignore.end())
    {
      return;
    }
    if (tag_ == "tableColumnTypes")
    {
      at_.colTypes.swap(header_);
      header_.clear();
    }
    else if (tag_ == "tableRowValues")
    {
      if (!row_.empty())
      {
        at_.tableRows.push_back(row_);
      }
      row_.clear();
    }
    else if (tag_ == "qualityParameter")
    {
      if (!(qp_.cvAcc == "MS:1000577" && parent_tag == "setQuality")) //set members get treated differently!
      {
        qps_.push_back(qp_);
        qp_ = QualityParameter();
      }
    }
    else if (tag_ == "attachment")
    {
      ats_.push_back(at_);
      at_ = Attachment();
    }
    else if (tag_ == "runQuality")
    {
      if (name_ == "")
      {
        name_ = run_id_;
        //~ name_ = String(UniqueIdGenerator::getUniqueId());
        //TODO give warning that a run should have a name cv!!!
      }
      registerRun(run_id_, name_);
      for (std::vector<QualityParameter>::const_iterator it = qps_.begin(); it != qps_.end(); ++it)
      {
        addRunQualityParameter(run_id_, *it);
      }
      for (std::vector<Attachment>::const_iterator it = ats_.begin(); it != ats_.end(); ++it)
      {
        addRunAttachment(run_id_, *it);
      }
      ats_.clear();
      qps_.clear();
    }
    else if (tag_ == "setQuality")
    {
      if (name_ == "")
      {
        name_ = run_id_;
        //~ name_ = String(UniqueIdGenerator::getUniqueId());
        //TODO give warning that a run should have a name cv!!!
      }
      registerSet(run_id_, name_, names_);
      for (std::vector<QualityParameter>::const_iterator it = qps_.begin(); it != qps_.end(); ++it)
      {
        addSetQualityParameter(run_id_, *it);
      }
      for (std::vector<Attachment>::const_iterator it = ats_.begin(); it != ats_.end(); ++it)
      {
        addSetAttachment(run_id_, *it);
      }
      ats_.clear();
      qps_.clear();
    }
  }

  void QcMLFile::store(const String& filename) const
  {
    //~ startProgress(0, 0, "storing qcML file");
    //~ progress_ = 0;
    //~ setProgress(++progress_);

    //~ file should either contain the complete stylesheet injection (including the stylesheet file preamble, the DOCTYPE definition and the stylesheet itself) or be empty
    std::string xslt = "";
    std::string xslt_ref = "";
    try
    {
      String xslt_file = File::find("XSL/QcML_report_sheet.xsl"); //TODO make this user defined pt.1
      std::ifstream in(xslt_file.c_str());
      xslt = std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      xslt = xslt.erase(0, xslt.find("\n") + 1);
      xslt_ref = "openms-qc-stylesheet"; //TODO make this user defined pt.2
    }
    catch (Exception::FileNotFound &)
    {
      warning(STORE, String("No qcml stylesheet found, result will not be viewable in a browser!"));
    }


    //open stream
    ofstream os(filename.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    os.precision(writtenDigits<double>(0.0));

    //~ setProgress(++progress_);
    //header & xslt

    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    if (!xslt_ref.empty())
    {
        os << "<?xml-stylesheet type=\"text/xml\" href=\"#" << xslt_ref << "\"?>\n";
        os << "<!DOCTYPE catelog [\n"
           << "  <!ATTLIST xsl:stylesheet\n"
           << "  id  ID  #REQUIRED>\n"
           << "  ]>\n";
    }
    os << "<qcML xmlns=\"https://github.com/qcML/qcml\" >\n"; //TODO creation date into schema!!

    //content runs
    std::set<String> keys;
    for (std::map<String, std::vector<QualityParameter> >::const_iterator it = runQualityQPs_.begin(); it != runQualityQPs_.end(); ++it)
    {
      keys.insert(it->first);
    }
    for (std::map<String, std::vector<Attachment> >::const_iterator it = runQualityAts_.begin(); it != runQualityAts_.end(); ++it)
    {
      keys.insert(it->first);
    }

    if (!keys.empty())
    {
      for (std::set<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
      {
        os << "\t<runQuality ID=\"" << String(*it) << "\">\n";
        std::map<String, std::vector<QualityParameter> >::const_iterator qpsit = runQualityQPs_.find(*it);
        if (qpsit != runQualityQPs_.end())
        {
          for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
          {
            os << qit->toXMLString(4);
          }
        }
        std::map<String, std::vector<Attachment> >::const_iterator attit = runQualityAts_.find(*it);
        if (attit != runQualityAts_.end())
        {
          for (std::vector<QcMLFile::Attachment>::const_iterator ait = attit->second.begin(); ait != attit->second.end(); ++ait)
          {
            os << ait->toXMLString(4); //TODO check integrity of reference to qp!
          }
        }
        os << "\t</runQuality>\n";
      }
    }

    //content sets
    keys.clear();
    for (std::map<String, std::vector<QualityParameter> >::const_iterator it = setQualityQPs_.begin(); it != setQualityQPs_.end(); ++it)
    {
      keys.insert(it->first);
    }
    for (std::map<String, std::vector<Attachment> >::const_iterator it = setQualityAts_.begin(); it != setQualityAts_.end(); ++it)
    {
      keys.insert(it->first);
    }

    if (!keys.empty())
    {
      for (std::set<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
      {
        os << "\t<setQuality ID=\"" << String(*it) << "\">\n";
        //~ TODO warn if key has no entries in members_

        //document set members
        std::map<String, std::set<String> >::const_iterator jt = setQualityQPs_members_.find(*it);

        if (jt != setQualityQPs_members_.end())
        {
          for (std::set<String>::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt)
          {
            std::map<String, std::vector<QualityParameter> >::const_iterator rq = runQualityQPs_.find(*kt);
            if (rq != runQualityQPs_.end())
            {
                QcMLFile::QualityParameter qp;
                qp.id = *kt; ///< Identifier
                qp.name = "set name"; ///< Name
                qp.cvRef = "QC"; ///< cv reference
                qp.cvAcc = "QC:0000005";
                for (std::vector<QualityParameter>::const_iterator qit = rq->second.begin(); qit != rq->second.end(); ++qit)
                {
                  ///<qualityParameter name="mzML file" ID="OTT0650-S44-A-Leber_1_run_name" cvRef="MS" accession="MS:1000577" value="OTT0650-S44-A-Leber_1"/>
                  if (qit->cvAcc == "MS:1000577")
                    qp.value = qit->value;
                }
                os << qp.toXMLString(4);
            }
            else
            {
              //TODO warn - no mzML file registered for this runQC
            }
          }
        }

        std::map<String, std::vector<QualityParameter> >::const_iterator qpsit = setQualityQPs_.find(*it);
        if (qpsit != setQualityQPs_.end())
        {
          for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
          {
            os << qit->toXMLString(4);
          }
        }

        std::map<String, std::vector<Attachment> >::const_iterator attit = setQualityAts_.find(*it);
        if (attit != setQualityAts_.end())
        {
          for (std::vector<QcMLFile::Attachment>::const_iterator ait = attit->second.begin(); ait != attit->second.end(); ++ait)
          {
            os << ait->toXMLString(4);
          }
        }
        os << "\t</setQuality>\n";
      }
    }
    os <<  "\t<cvList>\n";
    os <<  "\t<cv uri=\"http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" ID=\"psi_cv_ref\" fullName=\"PSI-MS\" version=\"3.41.0\"/>\n";
    os <<  "\t<cv uri=\"https://github.com/qcML/qcML-development/blob/master/cv/qc-cv.obo\" ID=\"qc_cv_ref\" fullName=\"QC-CV\" version=\"0.1.1\"/>\n";
    os <<  "\t<cv uri=\"http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/phenotype/unit.obo\" ID=\"uo_cv_ref\" fullName=\"unit\" version=\"1.0.0\"/>\n";
    os <<  "\t</cvList>\n";

    if (!xslt_ref.empty())
    {
      os << xslt << "\n";
    }

    os << "</qcML>\n";
  }

}

#pragma clang diagnostic pop

