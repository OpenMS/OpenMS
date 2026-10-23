// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/UnimodXMLHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <future>
#include <thread>

using namespace std;

namespace OpenMS
{

  UnimodXMLFile::UnimodXMLFile() :
    Internal::XMLFile()
  {

  }

  UnimodXMLFile::~UnimodXMLFile() = default;

  void UnimodXMLFile::load(const std::string& filename, vector<ResidueModification*> & modifications)
  {
    std::string file = File::find(filename);

    if (loadInSlices_(file, modifications))
    {
      return;
    }
    Internal::UnimodXMLHandler handler(modifications, file);
    parse_(file, &handler);
  }

  namespace
  {
    /// Are the end tags of \<umod:mod\> in [@p from, @p to) of @p text the only text like them (no comments, CDATA
    /// sections, processing instructions), and does every specificity have a position the handler knows?
    bool sliceable(const std::string& text, size_t from, size_t to)
    {
      for (size_t p = text.find('<', from); p < to; p = text.find('<', p + 1))
      {
        if (text[p + 1] == '!' || text[p + 1] == '?') return false;
      }
      // a position the handler does not know is reported per specificity: keep such messages in file order
      const std::string position_attr = "position=\"";
      for (size_t p = text.find(position_attr, from); p < to; p = text.find(position_attr, p + 1))
      {
        if (!std::isspace(static_cast<unsigned char>(text[p - 1]))) continue; // e.g. 'composition="'
        const size_t value = p + position_attr.size();
        const size_t value_end = text.find('"', value);
        if (value_end == std::string::npos) return false;
        const std::string pos = text.substr(value, value_end - value);
        if (pos != "Anywhere" && pos != "Protein N-term" && pos != "Protein C-term" && pos != "Any C-term" && pos != "Any N-term") return false;
      }
      return true;
    }
  }

  bool UnimodXMLFile::loadInSlices_(const std::string & file, vector<ResidueModification*> & modifications)
  {
    // unimod.xml (2.5 MB) takes most of the start-up of a tool that needs the ModificationsDB
    const Size max_slices = 8;
    const Size min_slice_bytes = 256 * 1024;
    const Size n_threads = std::max(1u, std::thread::hardware_concurrency());
    std::string text;
    {
      std::ifstream in(file, std::ios::binary | std::ios::ate);
      if (!in) return false;
      const std::streamoff size = in.tellg();
      if (size <= 0) return false;
      text.resize(size_t(size));
      in.seekg(0);
      if (!in.read(text.data(), size)) return false;
    }
    const std::string open_tag = "<umod:modifications>", close_tag = "</umod:modifications>", mod_end = "</umod:mod>";
    const size_t open = text.find(open_tag);
    if (open == std::string::npos) return false;
    const size_t begin = open + open_tag.size();
    const size_t end = text.find(close_tag, begin);
    if (end == std::string::npos) return false;
    const Size n_slices = std::min({max_slices, n_threads, Size((end - begin) / min_slice_bytes)});
    if (n_slices < 2) return false;
    std::vector<size_t> cuts{begin};
    for (Size k = 1; k < n_slices; ++k)
    {
      const size_t p = text.find(mod_end, std::max(cuts.back(), begin + (end - begin) * k / n_slices));
      if (p == std::string::npos || p + mod_end.size() >= end) break;
      cuts.push_back(p + mod_end.size());
    }
    cuts.push_back(end);
    if (cuts.size() < 3) return false;

    // each slice: everything before the section + its modifications + everything after the section; the end tags it
    // was cut at are those of elements if no slice holds a comment, CDATA section or processing instruction
    std::vector<std::vector<ResidueModification*>> parts(cuts.size() - 1);
    auto parse_slice = [&](Size i)
    {
      if (!sliceable(text, cuts[i], cuts[i + 1])) return false;
      // an error is reported by parsing the whole file again, not with the line numbers of the slice
      Logger::LogSinkGuard quiet(getThreadLocalLogFatal(), std::cerr);
      const std::string buffer = text.substr(0, begin) + text.substr(cuts[i], cuts[i + 1] - cuts[i]) + text.substr(end);
      Internal::UnimodXMLHandler handler(parts[i], file);
      UnimodXMLFile().parseBuffer_(buffer, &handler);
      return true;
    };
    std::vector<std::future<bool>> later;
    for (Size i = 1; i < parts.size(); ++i)
    {
      later.push_back(std::async(std::launch::async | std::launch::deferred, parse_slice, i));
    }
    bool ok = true;
    try
    {
      ok = parse_slice(0);
    }
    catch (...)
    {
      ok = false;
    }
    for (auto& f : later)
    {
      try
      {
        ok = f.get() && ok;
      }
      catch (...)
      {
        ok = false;
      }
    }
    for (auto& part : parts)
    {
      if (ok)
      {
        modifications.insert(modifications.end(), part.begin(), part.end());
      }
      else
      {
        for (ResidueModification* m : part) delete m;
      }
    }
    return ok;
  }

} // namespace OpenMS
