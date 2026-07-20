// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <filesystem>
#include <cerrno>
#ifdef OPENMS_WINDOWSPLATFORM
  #include <io.h>
  #include <share.h>
  #include <sys/stat.h>
  #include <fcntl.h>
#else
  #include <fcntl.h>
  #include <unistd.h>
#endif

namespace
{
  /**
    @brief Create a new, empty file next to @p target and return its path (empty on failure).

    The file is created exclusively, i.e. the call fails if the path already exists. That matters
    because the temporary lives in whatever directory the user's file is in: with a plain
    std::ofstream another process could pre-create the path as a symlink -- the name is derived from
    a timestamp and PID and is therefore guessable -- and have us truncate whatever it points at.
  */
  std::filesystem::path createExclusiveTempFile(const std::filesystem::path& target)
  {
    for (int attempt = 0; attempt < 32; ++attempt)
    {
      std::filesystem::path tmp = target;
      tmp += "." + OpenMS::File::getUniqueName(false) + ".tmp";
#ifdef OPENMS_WINDOWSPLATFORM
      int fd = -1;
      errno_t err = _wsopen_s(&fd, tmp.wstring().c_str(), _O_CREAT | _O_EXCL | _O_WRONLY | _O_BINARY, _SH_DENYNO, _S_IREAD | _S_IWRITE);
      if (err == 0 && fd != -1)
      {
        _close(fd);
        return tmp;
      }
      if (err != EEXIST)
      {
        break;
      }
#else
      int fd = ::open(tmp.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0666);
      if (fd != -1)
      {
        ::close(fd);
        return tmp;
      }
      if (errno != EEXIST)
      {
        break;
      }
#endif
    }
    return {};
  }
} // namespace

namespace OpenMS
{

  std::string writeXMLEscape(const std::string& to_escape)
  {
    return Internal::XMLHandler::writeXMLEscape(to_escape);
  }

  ParamXMLFile::ParamXMLFile() :
    XMLFile("/SCHEMAS/Param_1_8_0.xsd", "1.8.0")
  {
  }

  void ParamXMLFile::store(const std::string& filename, const Param& param) const
  {
    if (filename == "-")
    {
      writeXMLToStream(&std::cout, param);
      return;
    }

    namespace fs = std::filesystem;

    // Serialise into a temporary file and rename it over the target, so that a failure part way
    // through (full disk, I/O error, crash) leaves the previous, valid file untouched. Opening the
    // target directly would truncate it before we know the new content can be written at all -- for
    // an editor that file is often the user's only copy.
    std::error_code ec;
    // resolve symlinks so we replace the file that is linked to, not the link itself.
    // canonical() follows a whole chain of links; on failure (e.g. a dangling link) we keep the
    // path as given and simply create it.
    fs::path target(filename);
    if (fs::is_symlink(target, ec))
    {
      fs::path resolved = fs::canonical(target, ec);
      if (!ec)
      {
        target = resolved;
      }
      ec.clear();
    }
    // Refuse targets that must not be replaced by a rename. Opening the target directly used to
    // reject these implicitly; without the check we would happily remove a directory or overwrite a
    // write-protected file, because on POSIX rename() is governed by the directory's permissions.
    if (fs::exists(target, ec))
    {
      if (fs::is_directory(target, ec) || !File::writable(target.string()))
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }
    ec.clear();

    // keep the temporary next to the target so the rename stays within one filesystem
    fs::path tmp = createExclusiveTempFile(target);
    if (tmp.empty())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    {
      std::ofstream os(tmp, std::ofstream::out);
      if (!os)
      {
        fs::remove(tmp, ec);
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      try
      {
        writeXMLToStream(&os, param);
      }
      catch (...)
      {
        os.close();
        fs::remove(tmp, ec);
        throw;
      }
      os.close();
      if (!os) // flush/close failed, e.g. the disk filled up
      {
        fs::remove(tmp, ec);
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    // Carry over the permissions of the file we are replacing. Note that replacing by rename gives
    // the target a new inode, so ownership, ACLs, extended attributes and hard links to it are not
    // preserved -- the accepted trade-off for not truncating the only copy on a failed write (Qt's
    // QSaveFile behaves the same way).
    if (fs::exists(target, ec))
    {
      auto perms = fs::status(target, ec).permissions();
      if (!ec)
      {
        fs::permissions(tmp, perms, ec);
      }
    }

    fs::rename(tmp, target, ec);
    if (ec)
    { // POSIX replaces the target atomically. Windows refuses when it exists, so retry after
      // removing it -- but only for a plain file we already know we may overwrite, and never as a
      // blanket reaction to an arbitrary rename failure, which would delete the user's data and
      // then still fail.
      std::error_code probe;
      if (fs::is_regular_file(target, probe) && fs::remove(target, probe))
      {
        fs::rename(tmp, target, ec);
      }
    }
    if (ec)
    {
      std::error_code ignored;
      fs::remove(tmp, ignored);
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void ParamXMLFile::writeXMLToStream(std::ostream* os_ptr, const Param& param) const
  {
    // Note: For a long time the handling of 'getTrace()' was vulnerable to an unpruned tree (a path of nodes, but no entries in them), i.e.
    //       too many closing tags are written to the INI file, but no opening ones.
    //       This never mattered here, as removeAll() was fixed to prune the tree.
    // TODO: Nowadays this should be fixed and removeAll() might not be necessary.

    std::ostream& os = *os_ptr;

    os.precision(writtenDigits<double>(0.0));

    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    os << "<PARAMETERS version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/Param_1_8_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    std::string indentation = "  ";
    Param::ParamIterator it = param.begin();
    while (it != param.end())
    {
      //write opened/closed nodes
      const std::vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
      for (const Param::ParamIterator::TraceInfo& it2 : trace)
      {
        if (it2.opened) //opened node
        {
          std::string d = it2.description;
          //StringUtils::substitute(d, '"','\'');
          StringUtils::substitute(d, "\n", "#br#");
          //StringUtils::substitute(d, "<","&lt;");
          //StringUtils::substitute(d, ">","&gt;");
          os << indentation  << "<NODE name=\"" << writeXMLEscape(it2.name) << "\" description=\"" << writeXMLEscape(d) << "\">" << "\n";
          indentation += "  ";
        }
        else //closed node
        {
          indentation.resize(indentation.size() - 2);
          os << indentation << "</NODE>" << "\n";
        }
      }

      //write item
      if (it->value.valueType() != ParamValue::EMPTY_VALUE)
      {
        // we create a temporary copy of the tag list, since we remove certain tags while writing,
        // that will be represented differently in the xml
        std::set<std::string> tag_list = it->tags;
        ParamValue::ValueType value_type = it->value.valueType();
        bool stringParamIsFlag = false;

        //write opening tag
        switch (value_type)
        {
        case ParamValue::INT_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << R"(" type="int")";
          break;

        case ParamValue::DOUBLE_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << R"(" type="double")";
          break;

        case ParamValue::STRING_VALUE:
          if (tag_list.contains("input file"))
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << R"(" type="input-file")";
            tag_list.erase("input file");
          }
          else if (tag_list.contains("output file"))
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << R"(" type="output-file")";
            tag_list.erase("output file");
          }
          else if (tag_list.contains("output prefix"))
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << R"(" type="output-prefix")";
            tag_list.erase("output prefix");
          }

          else if (it->valid_strings.size() == 2 &&
          it->valid_strings[0] == "true" && it->valid_strings[1] == "false" &&
          it->value == "false")
          {
            stringParamIsFlag = true;
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << Internal::encodeTab(writeXMLEscape(it->value.toString())) << R"(" type="bool")";
          }
          else
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << Internal::encodeTab(writeXMLEscape(it->value.toString())) << R"(" type="string")";
          }
          break;

        case ParamValue::STRING_LIST:
          if (tag_list.contains("input file"))
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="input-file")";
            tag_list.erase("input file");
          }
          else if (tag_list.contains("output file"))
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="output-file")";
            tag_list.erase("output file");
          }
          else
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="string")";
          }
          break;

        case ParamValue::INT_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="int")";
          break;

        case ParamValue::DOUBLE_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="double")";
          break;

        default:
          break;
        }

        //replace all critical characters in description
        std::string d = it->description;
        //StringUtils::substitute(d, "\"","'");
        StringUtils::substitute(d, "\n", "#br#");
        //StringUtils::substitute(d, "<","&lt;");
        //StringUtils::substitute(d, ">","&gt;");
        os << " description=\"" << writeXMLEscape(d) << "\"";

        // required
        if (tag_list.contains("required"))
        {
          os << " required=\"true\"";
          tag_list.erase("required");
        }
        else
        {
          os << " required=\"false\"";
        }

        // advanced
        if (tag_list.contains("advanced"))
        {
          os << " advanced=\"true\"";
          tag_list.erase("advanced");
        }
        else
        {
          os << " advanced=\"false\"";
        }

        // tags
        if (!tag_list.empty())
        {
          std::string list;
          for (std::set<std::string>::const_iterator tag_it = tag_list.begin(); tag_it != tag_list.end(); ++tag_it)
          {
            if (!list.empty())
              list += ",";
            list += *tag_it;
          }
          os << " tags=\"" << writeXMLEscape(list) << "\"";
        }

        //restrictions
        // for boolean Flags they are implicitly given
        if (!stringParamIsFlag)
        {
          std::string restrictions;
          switch (value_type)
          {
            case ParamValue::INT_VALUE:
            case ParamValue::INT_LIST:
            {
              bool min_set = (it->min_int != -std::numeric_limits<Int>::max());
              bool max_set = (it->max_int != std::numeric_limits<Int>::max());
              if (max_set || min_set)
              {
                if (min_set)
                {
                  restrictions +=StringUtils::toStr(it->min_int);
                }
                restrictions += ':';
                if (max_set)
                {
                  restrictions +=StringUtils::toStr(it->max_int);
                }
              }
            }
              break;

            case ParamValue::DOUBLE_VALUE:
            case ParamValue::DOUBLE_LIST:
            {
              bool min_set = (it->min_float != -std::numeric_limits<double>::max());
              bool max_set = (it->max_float != std::numeric_limits<double>::max());
              if (max_set || min_set)
              {
                if (min_set)
                {
                  restrictions +=StringUtils::toStr(it->min_float);
                }
                restrictions += ':';
                if (max_set)
                {
                  restrictions +=StringUtils::toStr(it->max_float);
                }
              }
            }
              break;

            case ParamValue::STRING_VALUE:
            case ParamValue::STRING_LIST:
              if (!it->valid_strings.empty())
              {
                restrictions = StringUtils::concatenate(it->valid_strings, ",");
              }
              break;

            default:
              break;
          }
          // for files we store the restrictions as supported_formats
          if (!restrictions.empty())
          {
            if (it->tags.contains("input file") 
              || it->tags.contains("output file")
              || it->tags.contains("output prefix"))
            {
              os << " supported_formats=\"" << writeXMLEscape(restrictions) << "\"";
            }
            else
            {
              os << " restrictions=\"" << writeXMLEscape(restrictions) << "\"";
            }
          }
        }

        //finish opening tag
        switch (value_type)
        {
        case ParamValue::INT_VALUE:
        case ParamValue::DOUBLE_VALUE:
        case ParamValue::STRING_VALUE:
          os << " />" <<  "\n";
          break;

        case ParamValue::STRING_LIST:
        {
          os << ">" <<  "\n";
          const std::vector<std::string>& list = it->value;
          for (Size i = 0; i < list.size(); ++i)
          {
            os << indentation << "  <LISTITEM value=\"" << Internal::encodeTab(writeXMLEscape(list[i])) << "\"/>" << "\n";
          }
          os << indentation << "</ITEMLIST>" << "\n";
        }
        break;

        case ParamValue::INT_LIST:
        {
          os << ">" <<  "\n";
          const IntList& list = it->value;
          for (Size i = 0; i < list.size(); ++i)
          {
            os << indentation << "  <LISTITEM value=\"" << list[i] << "\"/>" << "\n";
          }
          os << indentation << "</ITEMLIST>" << "\n";
        }
        break;

        case ParamValue::DOUBLE_LIST:
        {
          os << ">" <<  "\n";
          const DoubleList& list = it->value;
          for (Size i = 0; i < list.size(); ++i)
          {
            os << indentation << "  <LISTITEM value=\"" << list[i] << "\"/>" << "\n";
          }
          os << indentation << "</ITEMLIST>" << "\n";
        }
        break;

        default:
          break;
        }
      }
      ++it;
    }

    // if we had tags ...
    if (param.begin() != param.end())
    {
      //close remaining tags
      const std::vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
      for (std::vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
      {
        Size ss = indentation.size();
        indentation.resize(ss - 2);
        os << indentation << "</NODE>" << "\n";
      }
    }

    os << "</PARAMETERS>" << std::endl; // forces a flush
  }

  void ParamXMLFile::load(const std::string& filename, Param& param)
  {
    Internal::ParamXMLHandler handler(param, filename, schema_version_);
    parse_(filename, &handler);

    // Reject documents that are not parameter files. parse_() validates nothing and the handler
    // ignores unknown elements, so without this any well-formed XML "loads" as an empty Param and
    // reports success -- in an editor that means opening an unrelated file appears to work and
    // saving then replaces it. See OpenMS/OpenMS#9768.
    // Checking for a PARAMETERS element rather than a PARAMETERS *root* is deliberate: CTD files
    // nest it inside a 'tool' root and are loaded through here as well.
    if (!handler.hasSeenParametersElement())
    {
      // stray ITEM elements may already have been imported, so do not leave the caller's Param
      // half-populated if it swallows the exception
      param.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                  "Not a parameter file: no 'PARAMETERS' element found.");
    }
  }

}
