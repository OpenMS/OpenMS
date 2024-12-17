// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer, Axel Walter $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <QtCore/QFileInfo>

#include <fstream>
#include <set>

using namespace std;

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

  QcMLFile::QualityParameter::QualityParameter(const QualityParameter& rhs) = default;

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
    if (!value.empty())
    {
      s += " value=\"" + value + "\"";
    }
    if (!unitRef.empty())
    {
      s += " unitRef=\"" + unitRef + "\"";
    }
    if (!unitAcc.empty())
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

  QcMLFile::Attachment::Attachment(const Attachment& rhs) = default;

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

  String QcMLFile::Attachment::toCSVString(const String& separator) const
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
    if (!value.empty())
    {
      s += " value=\"" + value + "\"";
    }
    if (!unitRef.empty())
    {
      s += " unitRef=\"" + unitRef + "\"";
    }
    if (!unitAcc.empty())
    {
      s += " unitAcc=\"" + unitAcc + "\"";
    }
    if (!qualityRef.empty())
    {
      s += " qualityParameterRef=\"" + qualityRef + "\"";
    }

    if (!binary.empty())
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
      for (String& it : copy)
      {
        it.substitute(String(" "), String("_"));
      }

      s += ListUtils::concatenate(copy, " ").trim();
      s += "</tableColumnTypes>\n";
      for (std::vector<std::vector<String> >::const_iterator it = tableRows.begin(); it != tableRows.end(); ++it)
      {
        s += indent + "\t" + "<tableRowValues>";

        std::vector<String> copy_row = *it;
        for (String& sit : copy_row)
        {
          sit.substitute(String(" "), String("_"));
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
    XMLHandler("", "0.7"), XMLFile("/SCHEMAS/qcml.xsd", "0.7"), ProgressLogger() //TODO keep version up-to-date
  {
  }

  QcMLFile::~QcMLFile() = default;

  void QcMLFile::addRunQualityParameter(const String& run_id, const QualityParameter& qp)
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

  void QcMLFile::addSetQualityParameter(const String& set_id, const QualityParameter& qp)
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

  void QcMLFile::addRunAttachment(const String& run_id, const Attachment& at)
  {
    runQualityAts_[run_id].push_back(at); //TODO permit AT without a QP (or enable orphan write out in store),redundancy check
  }

  void QcMLFile::addSetAttachment(const String& run_id, const Attachment& at)
  {
    setQualityAts_[run_id].push_back(at); //TODO add file QP to set member
  }

  void QcMLFile::getRunNames(std::vector<String>& ids) const
  {
    ids.clear();
    for (const auto& m : run_Name_ID_map_) ids.push_back(m.first);
  }

  void QcMLFile::getRunIDs(std::vector<String>& ids) const
  {
    ids.clear();
    for (const auto& m : runQualityQPs_) ids.push_back(m.first);
  }

  bool QcMLFile::existsRun(const String& filename, bool checkname) const
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

  bool QcMLFile::existsSet(const String& filename, bool checkname) const
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

  void QcMLFile::existsRunQualityParameter(const String& filename, const String& qpname, std::vector<String>& ids) const
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

  void QcMLFile::existsSetQualityParameter(const String& filename, const String& qpname, std::vector<String>& ids) const
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

  void QcMLFile::removeQualityParameter(const String& r, std::vector<String>& ids)
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

  void QcMLFile::removeAttachment(const String& r, std::vector<String>& ids, const String& at)
  {
    bool not_all = !at.empty();
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

  void QcMLFile::removeAttachment(const String& r, const String& at)
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

  void QcMLFile::removeAllAttachments(const String& at)
  {
    for (std::map<String, std::vector<Attachment> >::iterator it = runQualityAts_.begin(); it != runQualityAts_.end(); ++it)
    {
      removeAttachment(it->first, at);
    }
  }

  void QcMLFile::registerRun(const String& id, const String& name)
  {
    runQualityQPs_[id] = std::vector<QualityParameter>();
    runQualityAts_[id] = std::vector<Attachment>();
    run_Name_ID_map_[name] = id;
  }

  void QcMLFile::registerSet(const String& id, const String& name, const std::set<String>& names)
  {
    setQualityQPs_[id] = std::vector<QualityParameter>();
    setQualityAts_[id] = std::vector<Attachment>();
    set_Name_ID_map_[name] = id;
    setQualityQPs_members_[id] = names;
  }

  void QcMLFile::merge(const QcMLFile& addendum, const String& setname)
  {
    //~ runs (and create set if setname not empty)
    for (std::map<String, std::vector<QualityParameter> >::const_iterator it = addendum.runQualityQPs_.begin(); it != addendum.runQualityQPs_.end(); ++it)
    {
      runQualityQPs_[it->first].insert(runQualityQPs_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(runQualityQPs_[it->first].begin(), runQualityQPs_[it->first].end());
      runQualityQPs_[it->first].erase(std::unique(runQualityQPs_[it->first].begin(), runQualityQPs_[it->first].end()), runQualityQPs_[it->first].end());
      if (!setname.empty())
      {
        setQualityQPs_members_[setname].insert(it->first);
      }
    }
    for (std::map<String, std::vector<Attachment> >::const_iterator it = addendum.runQualityAts_.begin(); it != addendum.runQualityAts_.end(); ++it)
    {
      runQualityAts_[it->first].insert(runQualityAts_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(runQualityAts_[it->first].begin(), runQualityAts_[it->first].end());
      runQualityAts_[it->first].erase(std::unique(runQualityAts_[it->first].begin(), runQualityAts_[it->first].end()), runQualityAts_[it->first].end());
      if (!setname.empty())
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

  String QcMLFile::exportAttachment(const String& filename, const String& qpname) const
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

  String QcMLFile::exportQP(const String& filename, const String& qpname) const
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

  String QcMLFile::exportQPs(const String& filename, const StringList& qpnames) const
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
      for (const QualityParameter& it : found->second)
      {
        if (it.cvAcc == "QC:0000043" || it.cvAcc == "QC:0000044" || it.cvAcc == "QC:0000045" || it.cvAcc == "QC:0000046" || it.cvAcc == "QC:0000047")
        {
          cvs_table["id"][it.name.prefix(' ')] = it.value;
        }
        else if (it.cvAcc == "QC:0000053" || it.cvAcc == "QC:0000054" || it.cvAcc == "QC:0000055" || it.cvAcc == "QC:0000056" || it.cvAcc == "QC:0000057")
        {
          cvs_table["ms2"][it.name.prefix(' ')] = it.value;
        }
      }
      if (!cvs_table.empty())
      {
        return map2csv(cvs_table, "\t");
      }
    }

    return "";
  }

  void QcMLFile::collectSetParameter(const String& setname, const String& qp, std::vector<String>& ret)
  {
    for (std::set<String>::const_iterator it = setQualityQPs_members_[setname].begin(); it != setQualityQPs_members_[setname].end(); ++it)
    {
      for (const QualityParameter& jt : runQualityQPs_[*it])
      {
        if (jt.cvAcc == qp)
        {
          ret.push_back(jt.value);
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
      if (name_.empty())
      {
        name_ = run_id_;
        //~ name_ = String(UniqueIdGenerator::getUniqueId());
        //TODO give warning that a run should have a name cv!!!
      }
      registerRun(run_id_, name_);
      for (const QualityParameter& it : qps_)
      {
        addRunQualityParameter(run_id_, it);
      }
      for (const Attachment& it : ats_)
      {
        addRunAttachment(run_id_, it);
      }
      ats_.clear();
      qps_.clear();
    }
    else if (tag_ == "setQuality")
    {
      if (name_.empty())
      {
        name_ = run_id_;
        //~ name_ = String(UniqueIdGenerator::getUniqueId());
        //TODO give warning that a run should have a name cv!!!
      }
      registerSet(run_id_, name_, names_);
      for (const QualityParameter& it : qps_)
      {
        addSetQualityParameter(run_id_, it);
      }
      for (const Attachment& it : ats_)
      {
        addSetAttachment(run_id_, it);
      }
      ats_.clear();
      qps_.clear();
    }
  }

  float calculateSNmedian(const MSSpectrum& spec, bool norm = true)
  {
    if (spec.empty())
    {
      return 0;
    }
    vector<UInt> intensities;
    for (auto& pt : spec)
    {
      intensities.push_back(pt.getIntensity());
    }
    float median = Math::median(intensities.begin(), intensities.end());
    
    float maxi = spec.back().getIntensity();
    if (!norm)
    {
      float sn_by_max2median = maxi / median;
      return sn_by_max2median;
    }

    float sign_int= 0;
    float nois_int = 0;
    size_t sign_cnt= 0;
    size_t nois_cnt = 0;
    for (const Peak1D& pt : spec)
    {
      if (pt.getIntensity() <= median)
      {
        ++nois_cnt;
        nois_int += pt.getIntensity();
      }
      else
      {
        ++sign_cnt;
        sign_int += pt.getIntensity();
      }
    }
    if (sign_cnt == 0 || nois_cnt == 0 || nois_int <= 0)
    {
      return 0;
    }
    return (sign_int / sign_cnt) / (nois_int / nois_cnt);
  }

  void QcMLFile::collectQCData(vector<ProteinIdentification>& prot_ids,
                               vector<PeptideIdentification>& pep_ids,
                               const FeatureMap& feature_map,
                               const ConsensusMap& consensus_map,
                               const String& inputfile_raw, 
                               const bool remove_duplicate_features,
                               const MSExperiment& exp)
  {
      // fetch vocabularies
      ControlledVocabulary cv;
      cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));
      cv.loadFromOBO("QC", File::find("/CV/qc-cv-legacy.obo"));
      //-------------------------------------------------------------
      // MS acquisition
      //------------------------------------------------------------
      String base_name = QFileInfo(QString::fromStdString(inputfile_raw)).baseName();

      UInt min_mz = std::numeric_limits<UInt>::max();
      UInt max_mz = 0;
      std::map<Size, UInt> mslevelcounts;
      
      registerRun(base_name,base_name); //TODO use UIDs
      //---base MS aquisition qp
      String msaq_ref = base_name + "_msaq";
      QcMLFile::QualityParameter qp;
      qp.id = msaq_ref; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000004";
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "mzML file"; ///< Name
      }
      addRunQualityParameter(base_name, qp);
      
      //---file origin qp
      qp = QcMLFile::QualityParameter();
      qp.name = "mzML file"; ///< Name
      qp.id = base_name + "_run_name"; ///< Identifier
      qp.cvRef = "MS"; ///< cv reference
      qp.cvAcc = "MS:1000577";
      qp.value = base_name;
      addRunQualityParameter(base_name, qp);
      
      qp = QcMLFile::QualityParameter();
      qp.name = "instrument model"; ///< Name
      qp.id = base_name + "_instrument_name"; ///< Identifier
      qp.cvRef = "MS"; ///< cv reference
      qp.cvAcc = "MS:1000031";
      qp.value = exp.getInstrument().getName();
      addRunQualityParameter(base_name, qp);    

      qp = QcMLFile::QualityParameter();
      qp.name = "completion time"; ///< Name
      qp.id = base_name + "_date"; ///< Identifier
      qp.cvRef = "MS"; ///< cv reference
      qp.cvAcc = "MS:1000747";
      qp.value = exp.getDateTime().getDate();
      addRunQualityParameter(base_name, qp);

      //---precursors and SN
      QcMLFile::Attachment at;
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000044";
      at.qualityRef = msaq_ref;
      at.id = base_name + "_precursors"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "precursors"; ///< Name
      }

      at.colTypes.emplace_back("MS:1000894_[sec]");  // RT
      at.colTypes.emplace_back("MS:1000040");  // MZ
      at.colTypes.emplace_back("MS:1000041");  // charge
      at.colTypes.emplace_back("S/N");  // S/N
      at.colTypes.emplace_back("peak count");  // peak count
      
      for (Size i = 0; i < exp.size(); ++i)
      {
        mslevelcounts[exp[i].getMSLevel()]++;
        if (exp[i].getMSLevel() == 2)
        {
          if (exp[i].getPrecursors().front().getMZ() < min_mz)
          {
            min_mz = exp[i].getPrecursors().front().getMZ();
          }
          if (exp[i].getPrecursors().front().getMZ() > max_mz)
          {
            max_mz = exp[i].getPrecursors().front().getMZ();
          }
          std::vector<String> row;
          row.emplace_back(exp[i].getRT());
          row.emplace_back(exp[i].getPrecursors().front().getMZ());
          row.emplace_back(exp[i].getPrecursors().front().getCharge());
          row.emplace_back(calculateSNmedian(exp[i]));
          row.emplace_back(exp[i].size());
          at.tableRows.push_back(row);
        }
      }
      addRunAttachment(base_name, at);

      //---aquisition results qp
      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000006"; ///< cv accession for "aquisition results"
      qp.id = base_name + "_ms1aquisition"; ///< Identifier
      qp.value = String(mslevelcounts[1]);
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "number of ms1 spectra"; ///< Name
      }
      addRunQualityParameter(base_name, qp);
      

      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000007"; ///< cv accession for "aquisition results"
      qp.id = base_name + "_ms2aquisition"; ///< Identifier
      qp.value = String(mslevelcounts[2]);
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "number of ms2 spectra"; ///< Name
      }
      addRunQualityParameter(base_name, qp);

      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000008"; ///< cv accession for "aquisition results"
      qp.id = base_name + "_Chromaquisition"; ///< Identifier
      qp.value = String(exp.getChromatograms().size());
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "number of chromatograms"; ///< Name
      }
      addRunQualityParameter(base_name, qp);
      
      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000009";
      at.qualityRef = msaq_ref;
      at.id = base_name + "_mzrange"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "MS MZ aquisition ranges"; ///< Name
      }

      at.colTypes.emplace_back("QC:0000010"); //MZ
      at.colTypes.emplace_back("QC:0000011"); //MZ
      std::vector<String> rowmz;
      rowmz.emplace_back(min_mz);
      rowmz.emplace_back(max_mz);
      at.tableRows.push_back(rowmz);
      addRunAttachment(base_name, at);

      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000012";
      at.qualityRef = msaq_ref;
      at.id = base_name + "_rtrange"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "MS RT aquisition ranges"; ///< Name
      }

      at.colTypes.emplace_back("QC:0000013"); //MZ
      at.colTypes.emplace_back("QC:0000014"); //MZ
      std::vector<String> rowrt;
      rowrt.emplace_back(exp.begin()->getRT());
      rowrt.emplace_back(exp.getSpectra().back().getRT());
      at.tableRows.push_back(rowrt);
      addRunAttachment(base_name, at);
      

      //---ion current stability ( & tic ) qp
      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000022";
      at.qualityRef = msaq_ref;
      at.id = base_name + "_tics"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "MS TICs"; ///< Name
      }

      at.colTypes.emplace_back("MS:1000894_[sec]");
      at.colTypes.emplace_back("MS:1000285");
      Size below_10k = 0;
      std::vector<OpenMS::Chromatogram> chroms = exp.getChromatograms();
      if (!chroms.empty()) //real TIC from the mzML
      {
        for (Size t = 0; t < chroms.size(); ++t)
        {
          if (chroms[t].getChromatogramType() == ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM)
          {
            for (Size i = 0; i < chroms[t].size(); ++i)
            {
              double sum = chroms[t][i].getIntensity();
              if (sum < 10000)
              {
                ++below_10k;
              }
              std::vector<String> row;
              row.emplace_back(chroms[t][i].getRT() * 60);
              row.emplace_back(sum);
              at.tableRows.push_back(row);
            }
            break;  // what if there are more than one? should generally not be though ...
          }
        }
        addRunAttachment(base_name, at);

        qp = QcMLFile::QualityParameter();
        qp.id = base_name + "_ticslump"; ///< Identifier
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000023";
        qp.value = String((100 / exp.size()) * below_10k);
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "percentage of tic slumps"; ///< Name
        }
        addRunQualityParameter(base_name, qp);
      }

      // -- reconstructed TIC or RIC from the MS1 intensities
      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000056";
      at.qualityRef = msaq_ref;
      at.id = base_name + "_rics"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "MS RICs"; ///< Name
      }

      at.colTypes.emplace_back("MS:1000894_[sec]");
      at.colTypes.emplace_back("MS:1000285");
      at.colTypes.emplace_back("S/N");
      at.colTypes.emplace_back("peak count");
      Size prev = 0;
      below_10k = 0;
      Size jumps = 0;
      Size drops = 0;
      Size fact = 10;
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() == 1)
        {
          UInt sum = 0;
          for (Size j = 0; j < exp[i].size(); ++j)
          {
            sum += exp[i][j].getIntensity();
          }
          if (prev > 0 && sum > fact * prev)  // no jumps after complete drops (or [re]starts)
          {
            ++jumps;
          }
          else if (sum < fact*prev)
          {
            ++drops;
          }
          if (sum < 10000)
          {
            ++below_10k;
          }
          prev = sum;
          std::vector<String> row;
          row.emplace_back(exp[i].getRT());
          row.emplace_back(sum);
          row.emplace_back(calculateSNmedian(exp[i]));
          row.emplace_back(exp[i].size());
          at.tableRows.push_back(row);
        }
      }
      addRunAttachment(base_name, at);

      qp = QcMLFile::QualityParameter();
      qp.id = base_name + "_ricslump"; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000057";
      qp.value = String((100 / exp.size()) * below_10k);
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "percentage of ric slumps"; ///< Name
      }
      addRunQualityParameter(base_name, qp);

      qp = QcMLFile::QualityParameter();
      qp.id = base_name + "_ricjump"; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000059";
      qp.value = String(jumps);
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "IS-1A"; ///< Name
      }
      addRunQualityParameter(base_name, qp);

      qp = QcMLFile::QualityParameter();
      qp.id = base_name + "_ricdump"; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000060";
      qp.value = String(drops);
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "IS-1B"; ///< Name
      }
      addRunQualityParameter(base_name, qp);

      //---injection times MSn
      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000018";
      at.qualityRef = msaq_ref;
      at.id = base_name + "_ms2inj"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "MS2 injection time"; ///< Name
      }

      at.colTypes.emplace_back("MS:1000894_[sec]");
      at.colTypes.emplace_back("MS:1000927");
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() > 1 && !exp[i].getAcquisitionInfo().empty())
        {
          for (Size j = 0; j < exp[i].getAcquisitionInfo().size(); ++j)
          {
            if (exp[i].getAcquisitionInfo()[j].metaValueExists("MS:1000927"))
            {
              std::vector<String> row;
              row.emplace_back(exp[i].getRT());
              row.emplace_back(exp[i].getAcquisitionInfo()[j].getMetaValue("MS:1000927"));
              at.tableRows.push_back(row);
            }
          }
        }
      }
      if (!at.tableRows.empty())
      {
        addRunAttachment(base_name, at);
      }

      //-------------------------------------------------------------
      // MS  id
      //------------------------------------------------------------
      if (!prot_ids.empty() && !pep_ids.empty())
      {
        ProteinIdentification::SearchParameters params = prot_ids[0].getSearchParameters();
        vector<String> var_mods = params.variable_modifications;
        //~ boost::regex re("(?<=[KR])(?=[^P])");
      
        String msid_ref = base_name + "_msid";
        QcMLFile::QualityParameter qp;
        qp.id = msid_ref; ///< Identifier
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000025";
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "MS identification result details"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        at = QcMLFile::Attachment();
        at.cvRef = "QC"; ///< cv reference
        at.cvAcc = "QC:0000026";
        at.qualityRef = msid_ref;
        at.id = base_name + "_idsetting"; ///< Identifier
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
          at.name = term.name; ///< Name
        }
        catch (...)
        {
          at.name = "MS id settings"; ///< Name
        }
        
        at.colTypes.emplace_back("MS:1001013"); //MS:1001013 db name  MS:1001016 version  MS:1001020 taxonomy
        at.colTypes.emplace_back("MS:1001016");
        at.colTypes.emplace_back("MS:1001020");
        std::vector<String> row;
        row.emplace_back(prot_ids.front().getSearchParameters().db);
        row.emplace_back(prot_ids.front().getSearchParameters().db_version);
        row.emplace_back(prot_ids.front().getSearchParameters().taxonomy);
        at.tableRows.push_back(row);
        addRunAttachment(base_name, at);


        UInt spectrum_count = 0;
        Size peptide_hit_count = 0;
        Size protein_hit_count = 0;
        set<String> peptides;
        set<String> proteins;
        Size missedcleavages = 0;
        for (Size i = 0; i < pep_ids.size(); ++i)
        {
          if (!pep_ids[i].empty())
          {
            ++spectrum_count;
            peptide_hit_count += pep_ids[i].getHits().size();
            const vector<PeptideHit>& temp_hits = pep_ids[i].getHits();
            for (Size j = 0; j < temp_hits.size(); ++j)
            {
              peptides.insert(temp_hits[j].getSequence().toString());
            }
          }
        }
        for (set<String>::iterator it = peptides.begin(); it != peptides.end(); ++it)
        {
          for (String::const_iterator st = it->begin(); st != it->end() - 1; ++st)
          {
            if (*st == 'K' || *st == 'R')
            {
              ++missedcleavages;
            }
          }
        }

        for (Size i = 0; i < prot_ids.size(); ++i)
        {
          protein_hit_count += prot_ids[i].getHits().size();
          const vector<ProteinHit>& temp_hits = prot_ids[i].getHits();
          for (Size j = 0; j < temp_hits.size(); ++j)
          {
            proteins.insert(temp_hits[j].getAccession());
          }
        }
        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000037"; ///< cv accession
        qp.id = base_name + "_misscleave"; ///< Identifier
        qp.value = missedcleavages;
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "total number of missed cleavages"; ///< Name
        }
        addRunQualityParameter(base_name, qp);

        
        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000032"; ///< cv accession
        qp.id = base_name + "_totprot"; ///< Identifier
        qp.value = protein_hit_count;
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "total number of identified proteins"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000033"; ///< cv accession
        qp.id = base_name + "_totuniqprot"; ///< Identifier
        qp.value = String(proteins.size());
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "total number of uniquely identified proteins"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000029"; ///< cv accession
        qp.id = base_name + "_psms"; ///< Identifier
        qp.value = String(spectrum_count);
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "total number of PSM"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000030"; ///< cv accession
        qp.id = base_name + "_totpeps"; ///< Identifier
        qp.value = String(peptide_hit_count);
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "total number of identified peptides"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000031"; ///< cv accession
        qp.id = base_name + "_totuniqpeps"; ///< Identifier
        qp.value = String(peptides.size());
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "total number of uniquely identified peptides"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        at = QcMLFile::Attachment();
        at.cvRef = "QC"; ///< cv reference
        at.cvAcc = "QC:0000038";
        at.qualityRef = msid_ref;
        at.id = base_name + "_massacc"; ///< Identifier
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
          at.name = term.name; ///< Name
        }
        catch (...)
        {
          at.name = "delta ppm tables";
        }
        
        //~ delta ppm QC:0000039 RT MZ uniqueness ProteinID MS:1000885 target/decoy Score PeptideSequence MS:1000889 Annots string Similarity Charge UO:0000219 TheoreticalWeight UO:0000221 Oxidation_(M)
        at.colTypes.emplace_back("RT");
        at.colTypes.emplace_back("MZ");
        at.colTypes.emplace_back("Score");
        at.colTypes.emplace_back("PeptideSequence");
        at.colTypes.emplace_back("Charge");
        at.colTypes.emplace_back("TheoreticalWeight");
        at.colTypes.emplace_back("delta_ppm");
  //      at.colTypes.push_back("S/N");
        for (UInt w = 0; w < var_mods.size(); ++w)
        {
          at.colTypes.push_back(String(var_mods[w]).substitute(' ', '_'));
        }

        std::vector<double> deltas;
        //~ prot_ids[0].getSearchParameters();
        for (PeptideIdentification& pep_id : pep_ids)
        {
          if (!pep_id.getHits().empty())
          {
            std::vector<String> row;
            row.emplace_back(pep_id.getRT());
            row.emplace_back(pep_id.getMZ());
            PeptideHit tmp = pep_id.getHits().front();  //N.B.: depends on score & sort
            vector<UInt> pep_mods;
            for (UInt w = 0; w < var_mods.size(); ++w)
            {
              pep_mods.push_back(0);
            }
            for (const Residue& z : tmp.getSequence())
            {
              Residue res = z;
              String temp;
              if (res.isModified() && res.getModificationName() != "Carbamidomethyl")
              {
                temp = res.getModificationName() + " (" + res.getOneLetterCode()  + ")";
                //cout<<res.getModification()<<endl;
                for (UInt w = 0; w < var_mods.size(); ++w)
                {
                  if (temp == var_mods[w])
                  {
                    //cout<<temp;
                    pep_mods[w] += 1;
                  }
                }
              }
            }

            row.emplace_back(tmp.getScore());
            row.push_back(tmp.getSequence().toString().removeWhitespaces());
            row.emplace_back(tmp.getCharge());
            double mz = tmp.getSequence().getMZ(tmp.getCharge());
            row.emplace_back(mz);
            double dppm = (pep_id.getMZ()-mz)/(mz*(double)1e-6);
            row.emplace_back(dppm);
            deltas.push_back(dppm);
            for (UInt w = 0; w < var_mods.size(); ++w)
            {
              row.emplace_back(pep_mods[w]);
            }
            at.tableRows.push_back(row);
          }
        }
        addRunAttachment(base_name, at);
        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000040"; ///< cv accession
        qp.id = base_name + "_mean_delta"; ///< Identifier
        qp.value = String(OpenMS::Math::mean(deltas.begin(), deltas.end()));
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "mean delta ppm"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000041"; ///< cv accession
        qp.id = base_name + "_median_delta"; ///< Identifier
        qp.value = String(OpenMS::Math::median(deltas.begin(), deltas.end(), false));
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "median delta ppm"; ///< Name
        }
        addRunQualityParameter(base_name, qp);


        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000035"; ///< cv accession
        qp.id = base_name + "_ratio_id"; ///< Identifier
        qp.value = String(double(pep_ids.size()) / double(mslevelcounts[2]));
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "id ratio"; ///< Name
        }
        addRunQualityParameter(base_name, qp);
      }

      //-------------------------------------------------------------
      // MS quantitation
      //------------------------------------------------------------
      String msqu_ref = base_name + "_msqu";
      if (!feature_map.empty())
      {
        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000045"; ///< cv accession
        qp.id = msqu_ref; ///< Identifier
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "MS quantification result details"; ///< Name
        }
        addRunQualityParameter(base_name, qp);
        
        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000046"; ///< cv accession
        qp.id = base_name + "_feature_count"; ///< Identifier
        qp.value = String(feature_map.size());
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "number of features"; ///< Name
        }
        addRunQualityParameter(base_name, qp);      
      }

      if (!feature_map.empty() && !remove_duplicate_features)
      {
        QcMLFile::Attachment at;
        at = QcMLFile::Attachment();
        at.cvRef = "QC"; ///< cv reference
        at.cvAcc = "QC:0000047";
        at.qualityRef = msqu_ref;
        at.id = base_name + "_features"; ///< Identifier
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
          at.name = term.name; ///< Name
        }
        catch (...)
        {
          at.name = "features"; ///< Name
        }
        at.colTypes.emplace_back("MZ");
        at.colTypes.emplace_back("RT");
        at.colTypes.emplace_back("Intensity");
        at.colTypes.emplace_back("Charge");
        at.colTypes.emplace_back("Quality");
        at.colTypes.emplace_back("FWHM");
        at.colTypes.emplace_back("IDs");
        UInt fiter = 0;
        UInt ided = 0;
        //ofstream out(outputfile_name.c_str());
        while (fiter < feature_map.size())
        {
          std::vector<String> row;
          row.emplace_back(feature_map[fiter].getMZ());
          row.emplace_back(feature_map[fiter].getRT());
          row.emplace_back(feature_map[fiter].getIntensity());
          row.emplace_back(feature_map[fiter].getCharge());
          row.emplace_back(feature_map[fiter].getOverallQuality());
          row.emplace_back(feature_map[fiter].getWidth());
          row.emplace_back(feature_map[fiter].getPeptideIdentifications().size());
          if (!feature_map[fiter].getPeptideIdentifications().empty())
          {
            ++ided;
          }
          fiter++;
          at.tableRows.push_back(row);
        }     
        addRunAttachment(base_name, at);

        qp = QcMLFile::QualityParameter();
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:0000058"; ///< cv accession
        qp.id = base_name + "_idfeature_count"; ///< Identifier
        qp.value = ided;
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
          qp.name = term.name; ///< Name
        }
        catch (...)
        {
          qp.name = "number of identified features"; ///< Name
        }
        addRunQualityParameter(base_name, qp);
      }
      else if (!feature_map.empty() && remove_duplicate_features)
      {
        QcMLFile::Attachment at;
        at = QcMLFile::Attachment();
        at.cvRef = "QC"; ///< cv reference
        at.cvAcc = "QC:0000047";
        at.qualityRef = msqu_ref;
        at.id = base_name + "_features"; ///< Identifier
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
          at.name = term.name; ///< Name
        }
        catch (...)
        {
          at.name = "features"; ///< Name
        }
        
        at.colTypes.emplace_back("MZ");
        at.colTypes.emplace_back("RT");
        at.colTypes.emplace_back("Intensity");
        at.colTypes.emplace_back("Charge");
        FeatureMap map_out;
        UInt fiter = 0;
        while (fiter < feature_map.size())
        {
          FeatureMap map_tmp;
          for (UInt k = fiter; k <= feature_map.size(); ++k)
          {
            if (abs(feature_map[fiter].getRT() - feature_map[k].getRT()) < 0.1)
            {
              //~ cout << fiter << endl;
              map_tmp.push_back(feature_map[k]);
            }
            else
            {
              fiter = k;
              break;
            }
          }
          map_tmp.sortByMZ();
          UInt retif = 1;
          map_out.push_back(map_tmp[0]);
          while (retif < map_tmp.size())
          {
            if (abs(map_tmp[retif].getMZ() - map_tmp[retif - 1].getMZ()) > 0.01)
            {
              cout << "equal RT, but mass different" << endl;
              map_out.push_back(map_tmp[retif]);
            }
            retif++;
          }
        }
        addRunAttachment(base_name, at);
      }
      if (!consensus_map.empty())
      {
        at = QcMLFile::Attachment();
        qp.name = "consensuspoints"; ///< Name
        //~ qp.id = base_name + "_consensuses"; ///< Identifier
        qp.cvRef = "QC"; ///< cv reference
        qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession "feature mapper results"

        at.colTypes.emplace_back("Native_spectrum_ID");
        at.colTypes.emplace_back("DECON_RT_(sec)");
        at.colTypes.emplace_back("DECON_MZ_(Th)");
        at.colTypes.emplace_back("DECON_Intensity");
        at.colTypes.emplace_back("Feature_RT_(sec)");
        at.colTypes.emplace_back("Feature_MZ_(Th)");
        at.colTypes.emplace_back("Feature_Intensity");
        at.colTypes.emplace_back("Feature_Charge");
        for (ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit != consensus_map.end(); ++cmit)
        {
          const ConsensusFeature& CF = *cmit;
          for (ConsensusFeature::const_iterator cfit = CF.begin(); cfit != CF.end(); ++cfit)
          {
            std::vector<String> row;
            FeatureHandle FH = *cfit;
            row.emplace_back(CF.getMetaValue("spectrum_native_id"));
            row.emplace_back(CF.getRT()); row.emplace_back(CF.getMZ());
            row.emplace_back(CF.getIntensity());
            row.emplace_back(FH.getRT());
            row.emplace_back(FH.getMZ());
            row.emplace_back(FH.getCharge());
            at.tableRows.push_back(row);
          }
        }
        addRunAttachment(base_name, at);
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
      xslt = xslt.erase(0, xslt.find('\n') + 1);
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
        os << R"(<?xml-stylesheet type="text/xml" href="#)" << xslt_ref << "\"?>\n";
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
                for (const QualityParameter& qit : rq->second)
                {
                  if (qit.cvAcc == "MS:1000577")
                  {
                    qp.value = qit.value;
                  }
                }
                os << qp.toXMLString(4);
            }
            else
            {
              //TODO warn - no mzML file registered for this run
            }
          }
        }

        std::map<String, std::vector<QualityParameter> >::const_iterator qpsit = setQualityQPs_.find(*it);
        if (qpsit != setQualityQPs_.end())
        {
          for (const QcMLFile::QualityParameter& qit : qpsit->second)
          {
            os << qit.toXMLString(4);
          }
        }

        std::map<String, std::vector<Attachment> >::const_iterator attit = setQualityAts_.find(*it);
        if (attit != setQualityAts_.end())
        {
          for (const QcMLFile::Attachment& ait : attit->second)
          {
            os << ait.toXMLString(4);
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
