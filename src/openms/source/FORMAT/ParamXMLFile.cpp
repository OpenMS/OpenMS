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
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <cerrno>
#include <fcntl.h> // common to both platforms
#ifdef OPENMS_WINDOWSPLATFORM
  #include <io.h>
  #include <share.h>
  #include <sys/stat.h>
#else
  #include <unistd.h>
#endif

namespace
{
  /// An exclusively-created temporary file: its path together with the still-open write descriptor.
  struct ExclusiveTempFile
  {
    std::filesystem::path path;
    int fd = -1; ///< -1 means creation failed
  };

  /**
    @brief Create a new, empty file next to @p target, opened exclusively, and return its path
           together with the still-open write descriptor (@c fd == -1 on failure).

    The file is created with @c O_EXCL, i.e. the call fails if the path already exists. That matters
    because the temporary lives in whatever directory the user's file is in, and its name is derived
    from a timestamp and PID and is therefore guessable: without @c O_EXCL another process could
    pre-create the path as a symlink and have us truncate whatever it points at.

    The descriptor is deliberately kept open so the caller can write through it directly. Closing it
    and reopening the path with a fresh std::ofstream would reintroduce exactly that race -- another
    process could swap the name for a symlink between the two opens and redirect our write. Content
    is written in binary mode (LF line endings on all platforms), matching the repository's
    paramXML/INI reference files and the binary output of the base XMLFile writer.

    When @p restrictive is true (we are replacing an existing file) the file is created owner-only so
    the (possibly sensitive) parameters are not briefly world/group-readable while they are being
    written; the caller restores the target's permissions afterwards. When false (a new file) the
    normal umask-based default permissions are used, so a freshly created file is not unexpectedly
    restricted to owner-only.
  */
  ExclusiveTempFile createExclusiveTempFile(const std::filesystem::path& target, [[maybe_unused]] bool restrictive)
  {
    for (int attempt = 0; attempt < 32; ++attempt)
    {
      std::filesystem::path tmp = target;
      tmp += "." + OpenMS::File::getUniqueName(false) + ".tmp";
#ifdef OPENMS_WINDOWSPLATFORM
      int fd = -1;
      errno_t err = _wsopen_s(&fd, tmp.wstring().c_str(), _O_CREAT | _O_EXCL | _O_WRONLY | _O_BINARY, _SH_DENYRW, _S_IREAD | _S_IWRITE);
      if (err == 0 && fd != -1)
      {
        return {tmp, fd};
      }
      if (err != EEXIST)
      {
        break;
      }
#else
      // 0600 while replacing (owner-only during write); 0666 for a new file so the process umask
      // determines the final permissions as it would for an ordinary new file.
      const int mode = restrictive ? 0600 : 0666;
      int fd = ::open(tmp.c_str(), O_CREAT | O_EXCL | O_WRONLY, mode);
      if (fd != -1)
      {
        return {tmp, fd};
      }
      if (errno != EEXIST)
      {
        break;
      }
#endif
    }
    return {};
  }

  /**
    @brief Close a descriptor obtained from createExclusiveTempFile(). Returns 0 on success, -1 on error.

    The result must be checked: on NFS or when a quota is hit, the write error is not reported by the
    individual write() calls (the data is only buffered) but surfaces here at close time. Ignoring it
    would let an incomplete temporary be renamed over the previously-valid target.
  */
  int closeDescriptor(int fd)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    return _close(fd);
#else
    return ::close(fd);
#endif
  }

  /**
    @brief Non-mutating check that an @em existing path is writable.

    Deliberately not OpenMS::File::writable(): that probes a non-existent path by creating and then
    removing it, and if its internal existence check hits a transient stat error it treats the path
    as non-existent and opens it with a truncating std::ofstream -- which would destroy the very
    target we are trying to preserve. access()/_waccess_s only query permissions and never mutate;
    any error is reported as "not writable" (fail closed).
  */
  bool isExistingPathWritable(const std::filesystem::path& p)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    return _waccess_s(p.wstring().c_str(), 2) == 0; // 2 = write permission
#else
    return ::access(p.c_str(), W_OK) == 0;
#endif
  }

  /// Write the whole buffer to @p fd, looping over partial writes. Returns false on any error.
  bool writeAllToDescriptor(int fd, const char* data, std::size_t size)
  {
    std::size_t written = 0;
    while (written < size)
    {
      const std::size_t remaining = size - written;
#ifdef OPENMS_WINDOWSPLATFORM
      const unsigned int chunk = remaining > 0x40000000u ? 0x40000000u : static_cast<unsigned int>(remaining);
      const int n = _write(fd, data + written, chunk);
#else
      const ssize_t n = ::write(fd, data + written, remaining);
#endif
      if (n < 0)
      {
        if (errno == EINTR)
        {
          continue; // interrupted before any byte was written -- retry
        }
        return false;
      }
      if (n == 0)
      {
        return false; // no progress on a regular file -- treat as failure rather than spin forever
      }
      written += static_cast<std::size_t>(n);
    }
    return true;
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
    // through (full disk, I/O error) leaves the previous, valid file untouched. Opening the target
    // directly would truncate it before we know the new content can be written at all -- for an
    // editor that file is often the user's only copy. (Note: this is not a durability guarantee
    // across an OS/power crash -- there is no fsync of the file or its parent directory.)
    std::error_code ec;

    // Resolve symlinks so we replace the file the link points at, not any link in the chain. We
    // cannot rely on canonical() alone: it fails on a dangling link. Following the chain by hand
    // lets us create a not-yet-existing target while leaving every intermediate link intact -- each
    // relative hop is resolved against the directory of the link it came from. A cycle or an
    // over-long chain is rejected. Inspection errors fail closed (throw): silently continuing could
    // let a transient I/O error cause us to replace a symlink or special file we could not classify.
    fs::path target(filename);
    fs::file_status target_status;
    {
      int hops = 0;
      while (true)
      {
        target_status = fs::symlink_status(target, ec); // does not follow the final link
        if (ec)
        {
          if (ec == std::errc::no_such_file_or_directory)
          { // the path (or a parent component) does not exist: a new file to be created. Some
            // std::filesystem implementations report this via ec rather than file_type::not_found.
            target_status = fs::file_status(fs::file_type::not_found);
            ec.clear();
            break;
          }
          // any other error means we cannot classify the target -- fail closed rather than risk
          // replacing something we could not inspect
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                              "could not determine the type of the target");
        }
        if (target_status.type() != fs::file_type::symlink)
        {
          break; // reached a non-symlink terminal (regular, special, directory, or not-found)
        }
        if (++hops > 40) // more hops than any real chain -- almost certainly a cycle
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                              "the target is a cyclic or excessively long symlink chain");
        }
        fs::path link_target = fs::read_symlink(target, ec);
        if (ec)
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                              "could not read the target symlink");
        }
        target = link_target.is_absolute() ? link_target : target.parent_path() / link_target;
      }
    }
    // Only a regular file (or a not-yet-existing path) may be created or replaced by the rename.
    // Reject a directory, a special file (FIFO, socket, device) or a write-protected regular file:
    // the old stream-based code rejected these implicitly by failing to open them for writing, and
    // replacing them by rename would silently destroy them. The type comes from the cached
    // symlink_status above (no extra lookup); the writability check uses a non-mutating access()
    // that fails closed, so neither can fail open and destroy the target.
    const bool target_exists = fs::exists(target_status);
    if (target_exists && (!fs::is_regular_file(target_status) || !isExistingPathWritable(target)))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "the target exists but is not a writable regular file");
    }
    ec.clear();

    // Serialise the whole document up front. Doing this before touching the target means a
    // serialisation failure cannot leave a half-written temporary behind, and lets us write the
    // result through the exclusive descriptor in one go.
    std::ostringstream buffer;
    writeXMLToStream(&buffer, param);
    if (!buffer) // a stream-buffer/allocation failure can set badbit without throwing
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "failed to serialise the parameters");
    }
    const std::string data = buffer.str();

    // keep the temporary next to the target so the rename stays within one filesystem. When we are
    // replacing an existing file, create the temporary owner-only so its (possibly sensitive)
    // content is not briefly world/group-readable; when creating a new file, use the normal default
    // permissions (umask-based) so a new file is not unexpectedly restricted to 0600.
    ExclusiveTempFile tmp = createExclusiveTempFile(target, /* restrictive = */ target_exists);
    if (tmp.fd == -1)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "could not create a temporary file next to the target");
    }

    // Write through the descriptor we still hold from the exclusive create -- never by reopening
    // the path, which would let another process redirect the write via a symlink. The close result
    // is part of the check: on NFS / at a quota the write error surfaces only at close time.
    const bool write_ok = writeAllToDescriptor(tmp.fd, data.data(), data.size());
    const bool close_ok = (closeDescriptor(tmp.fd) == 0);
    if (!write_ok || !close_ok) // e.g. the disk filled up
    {
      fs::remove(tmp.path, ec);
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "could not write the temporary file");
    }

    // Carry over the permissions of the file we are replacing (the temporary was created owner-only
    // in that case). Note that replacing by rename gives the target a new inode, so ownership, ACLs,
    // extended attributes and hard links to it are not preserved -- the accepted trade-off for not
    // truncating the only copy on a failed write (Qt's QSaveFile behaves the same way). A failed
    // permission copy is left best-effort: the file keeps its owner-only temporary permissions,
    // which is safe (does not grant any group/other mode bits).
    if (target_exists)
    {
      auto perms = fs::status(target, ec).permissions();
      if (!ec)
      {
        fs::permissions(tmp.path, perms, ec);
      }
    }

    fs::rename(tmp.path, target, ec);
#ifdef OPENMS_WINDOWSPLATFORM
    bool original_removed = false;
    if (ec == std::errc::file_exists)
    { // POSIX (and Windows MoveFileEx with replacement) replace the target atomically; only a
      // Windows configuration that refuses to rename onto an existing file reaches here. Retry
      // after removing the target -- but only for that specific error, and only for a plain file we
      // already established we may overwrite, never as a blanket reaction to an arbitrary failure.
      std::error_code probe;
      if (fs::is_regular_file(target, probe) && fs::remove(target, probe))
      {
        original_removed = true;
        fs::rename(tmp.path, target, ec);
      }
    }
    if (ec && original_removed)
    { // The original was already deleted and the replacement could not be moved into place.
      // Removing tmp now would destroy both copies -- the very data loss this mechanism exists to
      // prevent. Leave the freshly written content under the temporary name and point at it.
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
          "the previous file was removed but the new content could not be put in place; it has been preserved as '" + tmp.path.string() + "'");
    }
#endif
    if (ec)
    {
      // The original (if any) is still intact, so discarding the freshly written temporary is safe.
      std::error_code ignored;
      fs::remove(tmp.path, ignored);
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "could not replace the target with the temporary file");
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
