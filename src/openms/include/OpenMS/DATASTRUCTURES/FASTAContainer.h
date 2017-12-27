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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_FASTACONTAINER_H
#define OPENMS_DATASTRUCTURES_FASTACONTAINER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <functional>
#include <fstream>
#include <memory>
#include <utility>
#include <vector>

namespace OpenMS
{

  struct TFI_File; ///< template parameter for file-based FASTA access
  struct TFI_Vector; ///< template parameter for vector-based FASTA access

  /**
  @brief This class allows for a chunk-wise single linear read over a (large) FASTA file, 
         with spurious (since potentially slow) access to earlier entries which are currently
         not in the active chunk.
  
  Internally uses FASTAFile class to read single sequences.

  FASTAContainer supports two template specializations FASTAContainer<TFI_File> and FASTAContainer<TFI_Vector>.
  
  FASTAContainer<TFI_File> will make FASTA entries available chunk-wise from start to end by loading it from a FASTA file.
  This avoids having to load the full file into memory. While loading, the container will
  memorize the file offsets of each entry, allowing to read an arbitrary i'th entry again from disk.
  If possible, only entries from the currently cached chunk should be queried, otherwise access will be slow.
  
  FASTAContainer<TFI_Vector> simply takes an existing vector of FASTAEntries and provides the same interface
  (with a potentially huge speed benefit over FASTAContainer<TFI_File> since it does not need disk access, but at the cost of memory).

  If an algorithm searches through a FASTA file linearly, you can use FASTAContainer<TFI_File> to pre-load a small chunk
  and start working, while loading the next chunk in a background thread and swap it in when the active chunk 
  was processed.

  */
template<typename TBackend>
class FASTAContainer; // prototype

/**
  @brief FASTAContainer<TFI_File> will make FASTA entries available chunk-wise from start to end by loading it from a FASTA file.
  This avoids having to load the full file into memory. While loading, the container will
  memorize the file offsets of each entry, allowing to read an arbitrary i'th entry again from disk.
  If possible, only entries from the currently cached chunk should be queried, otherwise access will be slow.

  Internally uses FASTAFile class to read single sequences.
*/
template<>
class FASTAContainer<TFI_File>
{
public:
  FASTAContainer() = delete;

  /// C'tor with FASTA filename
  FASTAContainer(const String& FASTA_file)
    : f_(),
    offsets_(),
    data_fg_(),
    data_bg_(),
    chunk_offset_(0)
  {
    f_.readStart(FASTA_file);
  }

  /// how many entries were read and got swapped out already
  size_t getChunkOffset() const
  {
    return chunk_offset_;
  }

  /** @brief Swaps in the background cache of entries, read previously via @p cacheChunk()
      
      If you call this function without a prior call to @p cacheChunk(), the cache will be empty.
      @return true if cache contains data; false if empty
      @note Should be invoked by a single thread, followed by a barrier to sync access of subsequent calls to chunkAt()
  */
  bool activateCache()
  {
    chunk_offset_ += data_fg_.size();
    data_fg_.swap(data_bg_);
    data_bg_.clear(); // just in case someone calls activateCache() multiple times...
    return !data_fg_.empty();
  }

  /** @brief Prefetch a new cache in the background, with up to @p suggestedSize entries (or fewer upon reaching EOF)

     Call @p activateCache() afterwards to make the data available via @p chunkAt() or @p readAt().
     @param suggested_size Number of FASTA entries to read from disk
     @return true if new data is available; false if background data is empty
  */
  bool cacheChunk(int suggested_size)
  {
    data_bg_.clear();
    data_bg_.reserve(suggested_size);
    FASTAFile::FASTAEntry p;
    for (int i = 0; i < suggested_size; ++i)
    {
      std::streampos spos = f_.position();
      if (!f_.readNext(p)) break;
      data_bg_.push_back(std::move(p));
      offsets_.push_back(spos);
    }
    return !data_bg_.empty();
  }

  /// number of entries in active cache
  size_t chunkSize() const
  {
    return data_fg_.size();
  }

  /** @brief Retrieve a FASTA entry at cache position @p pos (fast)
      
      Requires prior call to activateCache().
      Index @p pos must be smaller than chunkSize().

      @note: can be used by multiple threads at a time (until activateCache() is called)
  */
  const FASTAFile::FASTAEntry& chunkAt(size_t pos) const
  {
    return data_fg_[pos];
  }

  /** @brief Retrieve a FASTA entry at global position @pos (must not be behind the currently active chunk, but can be smaller)

    This query is fast, if @pos hits the currently active chunk, and slow (read from disk) for
    earlier entries. Can be used before reaching the end of the file,
    since it will reset the file position after its done reading (if reading from disk is required), but
    must not be used for entries beyond the active chunk (unseen data).
    
    @param protein Return value
    @param pos Absolute entry number in FASTA file
    @return true if reading was successful; false otherwise (e.g. EOF)
    @throw Exception::IndexOverflow if @p pos is beyond active chunk
    @note: not multi-threading safe (use chunkAt())!
  */
  bool readAt(FASTAFile::FASTAEntry& protein, size_t pos)
  {
    // check if position is currently cached...
    if (chunk_offset_ <= pos && pos < chunk_offset_ + chunkSize())
    {
      protein = data_fg_[pos - chunk_offset_];
      return true;
    }
    // read anew from disk...
    if (pos >= offsets_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pos, offsets_.size());
    }
    std::streampos spos = f_.position(); // save old position
    if (!f_.setPosition(offsets_[pos])) return false;
    bool r = f_.readNext(protein);
    f_.setPosition(spos); // restore old position
    return r;
  }

  /// is the FASTA file empty?
  bool empty() const
  { // trusting the FASTA file can be read...
    return f_.atEnd() && offsets_.empty();
  }

  /** @brief NOT the number of entries in the FASTA file, but merely the number of already read entries (since we do not know how many are still to come)

      @note Data in the background cache is included here, i.e. access to size()-1 using readAt() might be slow 
            if activateCache() was not called yet.
  */
  size_t size() const
  {
    return offsets_.size();
  }

private:
  FASTAFile f_; ///< FASTA file connection
  std::vector<std::streampos> offsets_; ///< internal byte offsets into FASTA file for random access reading of previous entries.
  std::vector<FASTAFile::FASTAEntry> data_fg_; ///< active (foreground) data
  std::vector<FASTAFile::FASTAEntry> data_bg_; ///< prefetched (background) data; will become the next active data
  size_t chunk_offset_; ///< number of entries before the current chunk
};

/**
@brief 
FASTAContainer<TFI_Vector> simply takes an existing vector of FASTAEntries and provides the same interface
with a potentially huge speed benefit over FASTAContainer<TFI_File> since it does not need disk access, but at the cost of memory.

*/
template<>
class FASTAContainer<TFI_Vector>
{
public:
  FASTAContainer() = delete;

  /** @brief C'tor for already existing data (by reference). 
  
   An internal reference will be kept. Make sure the data is not deleted during the lifetime of FASTAContainer

  */
  FASTAContainer(const std::vector<FASTAFile::FASTAEntry>& data)
    : data_(data)
  {
  }

  /// always 0, since this specialization requires/supports no chunking
  size_t getChunkOffset() const
  {
    return 0;
  }

  /** @brief no-op (since data is already fully available as vector)

     @return true only on the first call; false on subsequent calls
  */
  bool activateCache()
  {
    static int count = 0;
    ++count;
    return (count == 1); // only true on first call.
  }

  /** @brief no-op (since data is already fully available as vector)
     @return true only on the first call; false on subsequent calls
  */
  bool cacheChunk(int /*suggested_size*/)
  {
    // NOOP, since we already have all the data...
    static int count = 0;
    ++count;
    return (count == 1); // only true on first call.
  }

  /** @brief active data spans the full range, i.e. size of container
      
      @return the size of the underlying vector
  */
  size_t chunkSize() const
  {
    return data_.size();
  }

  /// fast access to chunked (i.e. all) entries
  const FASTAFile::FASTAEntry& chunkAt(size_t pos) const
  {
    return data_[pos];
  }

  /// fast access to an entry
  bool readAt(FASTAFile::FASTAEntry& protein, size_t pos) const
  {
    protein = data_[pos];
    return true;
  }

  /// calls empty() on underlying vector
  bool empty() const
  {
    return data_.empty();
  }

  /// calls size() on underlying vector
  size_t size() const
  {
    return data_.size();
  }

private:
  const std::vector<FASTAFile::FASTAEntry>& data_; ///< reference to existing data
};

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_FASTACONTAINER_H
