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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QcMLFile.h>
#include <fstream>
#include <iostream>
#include <algorithm>

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
    unitAcc()
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
    colTypes(rhs.colTypes),
    tableRows(rhs.tableRows)
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
      colTypes = rhs.colTypes;
      tableRows = rhs.tableRows;
    }
    return *this;
  }

  bool QcMLFile::QualityParameter::operator< (const QualityParameter & rhs) const
  {
    return name.toQString() < rhs.name.toQString();
  }

  bool QcMLFile::QualityParameter::operator> (const QualityParameter & rhs) const
  {
    return name.toQString() > rhs.name.toQString();
  }

  bool QcMLFile::QualityParameter::operator== (const QualityParameter & rhs) const
  {
    return name.toQString() == rhs.name.toQString();
  }

  String QcMLFile::QualityParameter::toXMLString(UInt indentation_level) const
  {
    String indent = String(indentation_level, '\t');
    String s = indent;
    s += "<QualityParameter";
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

    if ((!colTypes.empty()) && (!tableRows.empty()))
    {
      s += ">\n";
      s += indent + "\t" + "<TableColumnTypes>";
      s += StringList(colTypes).concatenate(" ").trim();
      s += "</TableColumnTypes>\n";
      for (std::vector< std::vector<String> >::const_iterator it = tableRows.begin(); it != tableRows.end(); ++it)
      {
        s += indent + "\t" + "<TableRowValues>";
        s += StringList(*it).concatenate(" ").trim();
        s += "</TableRowValues>\n";
      }
      s += indent + "</QualityParameter>\n";
    }
    else
    {
      s += "/>\n";
    }
    return s;
  }


  String QcMLFile::QualityParameter::toCSVString(const String separator) const
  {
    String s = "";
    if ((!colTypes.empty()) && (!tableRows.empty()))
    {
      s += StringList(colTypes).concatenate(separator).trim();
      s += "\n";
      for (std::vector< std::vector<String> >::const_iterator it = tableRows.begin(); it != tableRows.end(); ++it)
      {
        s += StringList(*it).concatenate(separator).trim();
        s += "\n";
      }
    }
    return s;
  }

  QcMLFile::Attachment::Attachment() :
    name(),
    value(),
    cvRef(),
    cvAcc(),
    unitRef(),
    unitAcc()
  {
  }

  QcMLFile::Attachment::Attachment(const Attachment& rhs) :
    name(rhs.name),
    value(rhs.value),
    cvRef(rhs.cvRef),
    cvAcc(rhs.cvAcc),
    unitRef(rhs.unitRef),
    unitAcc(rhs.unitAcc),
    binary(rhs.binary),
    qualityRef(rhs.qualityRef)
  {
  }

  QcMLFile::Attachment& QcMLFile::Attachment::operator=(const Attachment& rhs)
  {
    if (this != &rhs)
    {
      name = rhs.name;
      value = rhs.value;
      cvRef = rhs.cvRef;
      cvAcc = rhs.cvAcc;
      unitRef = rhs.unitRef;
      unitAcc = rhs.unitAcc;
      binary = rhs.binary;
      qualityRef = rhs.qualityRef;
    }
    return *this;
  }

  bool QcMLFile::Attachment::operator< (const Attachment& rhs) const
  {
    return name.toQString() < rhs.name.toQString();
  }

  bool QcMLFile::Attachment::operator> (const Attachment& rhs) const
  {
    return name.toQString() > rhs.name.toQString();
  }

  bool QcMLFile::Attachment::operator== (const Attachment& rhs) const
  {
    return name.toQString() == rhs.name.toQString();
  }

  String QcMLFile::Attachment::toXMLString(UInt indentation_level) const
  {
    String indent = String(indentation_level, '\t');
    String s = indent;
    s += "<Attachment ";
    s += " name=\"" + name + "\"" + " cvRef=\"" + cvRef + "\"" + " accession=\"" + cvAcc + "\"";
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
      s += indent + "</Attachment>\n";
    }
    else
    {
      s += "/>\n";
    }
    return s;
  }

  QcMLFile::QcMLFile():
    XMLHandler("", "0.3"), XMLFile("/SCHEMAS/qcml.xsd", "0.3"), ProgressLogger()
  {
  }

  QcMLFile::~QcMLFile()
  {

  }

  void QcMLFile::addQualityParameter(String r, QualityParameter qp)
  {
    runQualityQPs_[r].push_back(qp); //TODO redundancy check
  }

  void QcMLFile::removeQualityParameter(String r, String qp)
  {
    std::vector<QcMLFile::QualityParameter>::iterator qit = runQualityQPs_[r].begin();
    while (qit != runQualityQPs_[r].end())
    {
      if (qit->name == qp)
      {
        qit = runQualityQPs_[r].erase(qit);
      }
      else
      {
        ++qit;
      }
    }
  }

  void QcMLFile::merge(const QcMLFile & addendum)
  {
    for (std::map<String, std::vector< QualityParameter > >::const_iterator it = addendum.runQualityQPs_.begin(); it != addendum.runQualityQPs_.end(); ++it)
    {
      runQualityQPs_[it->first].insert(runQualityQPs_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(runQualityQPs_[it->first].begin(), runQualityQPs_[it->first].end() );
      runQualityQPs_[it->first].erase(std::unique(runQualityQPs_[it->first].begin(), runQualityQPs_[it->first].end()), runQualityQPs_[it->first].end() );
    }
    for (std::map<String, std::vector< Attachment > >::const_iterator it = addendum.runQualityAts_.begin(); it != addendum.runQualityAts_.end(); ++it)
    {
      runQualityAts_[it->first].insert(runQualityAts_[it->first].end(), it->second.begin(), it->second.end());
      std::sort(runQualityAts_[it->first].begin(), runQualityAts_[it->first].end() );
      runQualityAts_[it->first].erase(std::unique(runQualityAts_[it->first].begin(), runQualityAts_[it->first].end()), runQualityAts_[it->first].end() );

    }
  }

  void QcMLFile::addAttachment(String r, Attachment at)
  {
    runQualityAts_[r].push_back(at); //TODO redundancy check
  }

  String QcMLFile::exportQualityParameter(const String filename, const String qpname) const
  {
    std::map<String, std::vector< QcMLFile::QualityParameter > >::const_iterator qpsit = runQualityQPs_.find(filename);

    if (qpsit != runQualityQPs_.end())
    {
      for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        if (qpname == qit->name)
        {
          return qit->toCSVString("\t");
          //~ return qit->toXMLString(1);
        }
      }
    }
    return "";
  }

  String QcMLFile::existsQualityParameter(const String filename, const String qpname) const
  {
    String found = "";
    std::map<String, std::vector< QcMLFile::QualityParameter > >::const_iterator qpsit = runQualityQPs_.find(filename);

    if (qpsit != runQualityQPs_.end())
    {
      for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
      {
        if (qpname == qit->name)
        {
          return qit->id;
        }
      }
    }
    return found;
  }

  void QcMLFile::load(const String & filename)
  {
    //~ siz = 0;
    //Filename for error messages in XMLHandler
    file_ = filename;

    runQualityQPs_.clear(); // clear
    runQualityAts_.clear(); // clear

    parse_(filename, this);

    //reset more members
    //?
  }

  void QcMLFile::startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes)
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
      to_ignore.insert("TableColumnTypes"); // will be handled entirely in characters.
      to_ignore.insert("TableRowValues"); // ...
      to_ignore.insert("binary"); // ...
    }

    if (to_ignore.find(tag_) != to_ignore.end())
    {
      return;
    }

    String tmp_str;
    if (tag_ == "MzQualityML")
    {
      startProgress(0, 0, "loading qcML file");
      progress_ = 0;
      setProgress(++progress_);
    }
    else if (tag_ == "RunQuality")
    {
      setProgress(++progress_);
      qps_.clear();
      ats_.clear();
      qp_ = QualityParameter();
      at_ = Attachment();
      name_ = "";
      //for the run name wait for the qp with the right cv, otherwise use a uid
    }
    else if (tag_ == "QualityParameter")
    {
      if (parent_tag == "RunQuality")
      {
        optionalAttributeAsString_(qp_.value, attributes, "value");
        optionalAttributeAsString_(qp_.unitAcc, attributes, "unitAccession");
        optionalAttributeAsString_(qp_.unitRef, attributes, "unitCvRef");
        qp_.cvRef = attributeAsString_(attributes, "cvRef");
        qp_.cvAcc = attributeAsString_(attributes, "accession");
        qp_.id = attributeAsString_(attributes, "ID");
        qp_.name = attributeAsString_(attributes, "name");
        if (qp_.cvAcc == "MS:1000577") //id: MS:1000577 name: raw data file  - with value of the file name of the run
        {
          name_ = qp_.value;
        }
        //TODO add cvhandling for validation etc
      }
      else //SetQuality
      {
        //TODO
      }
    }
    else if (tag_ == "Attachment")
    {
        optionalAttributeAsString_(at_.value, attributes, "value");
        optionalAttributeAsString_(at_.unitAcc, attributes, "unitAccession");
        optionalAttributeAsString_(at_.unitRef, attributes, "unitCvRef");
        at_.cvRef = attributeAsString_(attributes, "cvRef");
        at_.cvAcc = attributeAsString_(attributes, "accession");
        at_.name = attributeAsString_(attributes, "name");
    }
    else if (tag_ == "SetQuality")
    {
      setProgress(++progress_);
      //TODO
    }
  }

  void QcMLFile::characters(const XMLCh * const chars, const XMLSize_t /*length*/)
  {
    if (tag_ == "TableRowValues")
    {
      String s = sm_.convert(chars);
      s.trim();
      if (!s.empty()) // always two notifications for a row, only the first one contains chars - dunno why
      {
        s.split(" ", row_);
      }
    }
    else if (tag_ == "TableColumnTypes")
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
    if (tag_ == "TableColumnTypes")
    {
      qp_.colTypes.swap(header_);
      header_.clear();
    }
    else if (tag_ == "TableRowValues")
    {
      if (!row_.empty())
      {
        qp_.tableRows.push_back(row_);
      }
      row_.clear();
    }
    else if (tag_ == "QualityParameter")
    {
      qps_.push_back(qp_);
      qp_ = QualityParameter();
    }
    else if (tag_ == "Attachment")
    {
      ats_.push_back(at_);
      at_ = Attachment();
    }
    else if (tag_ == "RunQuality")
    {
      if (name_ == "")
      {
        name_ = String(UniqueIdGenerator::getUniqueId());
        //TODO give warning that a run should have a name cv!!!
      }
      for (std::vector<QualityParameter>::const_iterator it = qps_.begin(); it != qps_.end(); ++it)
      {
        addQualityParameter(name_, *it);
      }
      for (std::vector<Attachment>::const_iterator it = ats_.begin(); it != ats_.end(); ++it)
      {
        addAttachment(name_, *it);
      }
      ats_.clear();
      qps_.clear();
    }
  }

  void QcMLFile::store(const String & filename) const
  {
    //~ startProgress(0, 0, "storing qcML file");
    //~ progress_ = 0;
    //~ setProgress(++progress_);

    //open stream
    ofstream os(filename.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    os.precision(writtenDigits<DoubleReal>());

    //~ setProgress(++progress_);
    //header & xslt
    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    os << "<?xml-stylesheet type=\"text/xml\" href=\"#stylesheet\"?>\n";
    os << "<!DOCTYPE catelog [\n"
       << "  <!ATTLIST xsl:stylesheet\n"
       << "  id  ID  #REQUIRED>\n"
       << "  ]>\n";
    os << "<MzQualityMLType>\n"; //TODO creation date into schema!!
    os << "<xsl:stylesheet id=\"stylesheet\" version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">\n"
       << "<xsl:template match=\"/\">\n"
       << "  <html>\n"
       << "  <body>\n"
       << "   <h2>The Quality Parameters</h2>\n"
       << "   <xsl:for-each select=\"MzQualityMLType/RunQuality/QualityParameter\">\n"
       << "    <xsl:value-of select=\"@name\" />: <xsl:value-of select=\"@value\" />\n"
       << "    <table border=\"1\">\n"
       << "    <tr bgcolor=\"#9acd32\">\n"
       << "     <xsl:call-template name=\"output-header\">\n"
       << "      <xsl:with-param name=\"list\"><xsl:value-of select=\"TableColumnTypes\" /></xsl:with-param>\n"
       << "     </xsl:call-template>\n"
       << "    </tr>\n"
       << "    <xsl:for-each select=\"TableRowValues\">\n"
       << "     <tr>\n"
       << "     <xsl:call-template name=\"output-row\">\n"
       << "      <xsl:with-param name=\"list\"><xsl:value-of select=\".\" /></xsl:with-param>\n"
       << "     </xsl:call-template></tr></xsl:for-each>\n"
       << "     </table><br/>\n"
       << "    </xsl:for-each>\n"
       << "   <h2>The Quality Plots</h2>\n"
       << "    <xsl:for-each select=\"MzQualityMLType/RunQuality/Attachment\">\n"
       << "     <img>\n"
       << "      <xsl:attribute name=\"src\">\n"
       << "       data:image/png;base64,<xsl:value-of select=\"binary\" />\n"
       << "      </xsl:attribute>\n"
       << "     </img> <br/>\n"
       << "    </xsl:for-each>\n"
       << "  </body>\n"
       << "  </html>\n"
       << "</xsl:template>\n"
       << "<xsl:template name=\"output-header\">\n"
       << "    <xsl:param name=\"list\" />\n"
       << "    <xsl:variable name=\"newlist\" select=\"concat(normalize-space($list), ' ')\" />\n"
       << "    <xsl:variable name=\"first\" select=\"substring-before($newlist, ' ')\" />\n"
       << "    <xsl:variable name=\"remaining\" select=\"substring-after($newlist, ' ')\" />\n"
       << "    <th><xsl:value-of select=\"$first\" /></th>\n"
       << "    <xsl:if test=\"$remaining\">\n"
       << "        <xsl:call-template name=\"output-header\">\n"
       << "            <xsl:with-param name=\"list\" select=\"$remaining\" />\n"
       << "        </xsl:call-template>\n"
       << "    </xsl:if>\n"
       << "</xsl:template>\n"
       << "<xsl:template name=\"output-row\">\n"
       << "    <xsl:param name=\"list\" />\n"
       << "    <xsl:variable name=\"newlist\" select=\"concat(normalize-space($list), ' ')\" />\n"
       << "    <xsl:variable name=\"first\" select=\"substring-before($newlist, ' ')\" />\n"
       << "    <xsl:variable name=\"remaining\" select=\"substring-after($newlist, ' ')\" />\n"
       << "    <td><xsl:value-of select=\"$first\" /></td>\n"
       << "    <xsl:if test=\"$remaining\">\n"
       << "        <xsl:call-template name=\"output-row\">\n"
       << "            <xsl:with-param name=\"list\" select=\"$remaining\" />\n"
       << "        </xsl:call-template>\n"
       << "    </xsl:if>\n"
       << "</xsl:template>\n"
       << "</xsl:stylesheet>\n";

    //content
    std::set<String> keys;
    for (std::map<String, std::vector< QualityParameter > >::const_iterator it = runQualityQPs_.begin(); it != runQualityQPs_.end(); ++it)
    {
      keys.insert(it->first);
    }
    for (std::map<String, std::vector< Attachment > >::const_iterator it = runQualityAts_.begin(); it != runQualityAts_.end(); ++it)
    {
      keys.insert(it->first);
    }

    if (!keys.empty())
    {
      for (std::set<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
      {
        os << "\t<RunQuality>\n";
        std::map<String, std::vector< QualityParameter > >::const_iterator qpsit = runQualityQPs_.find(*it);
        if (qpsit != runQualityQPs_.end())
        {
          for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qpsit->second.begin(); qit != qpsit->second.end(); ++qit)
          {
            os << qit->toXMLString(4);
          }
        }
        std::map<String, std::vector< Attachment > >::const_iterator attit = runQualityAts_.find(*it);
        if (attit != runQualityAts_.end())
        {
          for (std::vector<QcMLFile::Attachment>::const_iterator ait = attit->second.begin(); ait != attit->second.end(); ++ait)
          {
            os << ait->toXMLString(4);
          }
        }
        os << "\t</RunQuality>\n";
      }
    }
    //TODO SetQuality
    os << "</MzQualityMLType>\n";
  }

}
/*
<xsl:stylesheet id="stylesheet" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/">
  <html>
  <body>
   <h2>The Quality Parameters</h2>
		<xsl:for-each select="MzQualityMLType/RunQuality/QualityParameter">
			<xsl:value-of select="@name" />: <xsl:value-of select="@value" />
			<table border="1">
				<tr bgcolor="#9acd32">
   <xsl:call-template name="output-header">
     <xsl:with-param name="list"><xsl:value-of select="TableColumnTypes" /></xsl:with-param>
   </xsl:call-template>
       </tr>
       <tr>
   <xsl:call-template name="output-row">
     <xsl:with-param name="list"><xsl:value-of select="TableRowValues" /></xsl:with-param>
   </xsl:call-template>
       </tr>
     </table><br/>
		</xsl:for-each>
   <h2>The Quality Plots</h2>
      <xsl:for-each select="MzQualityMLType/RunQuality/Attachment">
        <img>
          <xsl:attribute name="src">
          data:image/png;base64,<xsl:value-of select="binary" />
          </xsl:attribute>
        </img> <br/>
     </xsl:for-each>
 </body>
 </html>
</xsl:template>

<xsl:template name="output-header">
    <xsl:param name="list" />
    <xsl:variable name="newlist" select="concat(normalize-space($list), ' ')" />
    <xsl:variable name="first" select="substring-before($newlist, ' ')" />
    <xsl:variable name="remaining" select="substring-after($newlist, ' ')" />
    <th><xsl:value-of select="$first" /></th>
    <xsl:if test="$remaining">
        <xsl:call-template name="output-header">
            <xsl:with-param name="list" select="$remaining" />
        </xsl:call-template>
    </xsl:if>
</xsl:template>
<xsl:template name="output-row">
    <xsl:param name="list" />
    <xsl:variable name="newlist" select="concat(normalize-space($list), ' ')" />
    <xsl:variable name="first" select="substring-before($newlist, ' ')" />
    <xsl:variable name="remaining" select="substring-after($newlist, ' ')" />
    <td><xsl:value-of select="$first" /></td>
    <xsl:if test="$remaining">
        <xsl:call-template name="output-row">
            <xsl:with-param name="list" select="$remaining" />
        </xsl:call-template>
    </xsl:if>
</xsl:template>
</xsl:stylesheet>*/
