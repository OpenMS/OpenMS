// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <memory>
#include <string>

namespace OpenMS::Internal
{
  /// Compile-time UTF-16 string literal, usable wherever a @c XMLCh* is expected
  /// (XMLCh is guaranteed to be char16_t-compatible; see StringManager.cpp).
  #define CONST_XMLCH(s) (u ## s)

  //Adapted from https://www.codeproject.com/articles/99551/redux-raii-adapter-for-xerces
  //Copyright 2010 Orjan Westin
  //Under BSD license
  //========================================================================================================
  template<typename T>
  class OPENMS_DLLAPI shared_xerces_ptr
  {
    // Function to release Xerces data type with a release member function
    template<typename U>
    static void doRelease_(U* item)
    {
      // Only release this if it has no owner
      if (nullptr == item->getOwnerDocument())
        item->release();
    }

    static void doRelease_(char* item);
    static void doRelease_(char16_t* item);

    // The actual data we're holding
    std::shared_ptr<T> item_;
  public:
    // Default constructor
    shared_xerces_ptr() = default;
    // Assignment constructor
    shared_xerces_ptr(T* item)
        : item_(item, doRelease_ )
    {}
    // Assignment of data to guard
    shared_xerces_ptr& operator=(T* item)
    {
      assign(item);
      return *this;
    }
    // Give up hold on data
    void reset()
    {
      item_.reset();
    }
    // Release currently held data, if any, to hold another
    void assign(T* item)
    {
      item_.reset(item, doRelease_ );
    }
    // Get pointer to the currently held data, if any
    T* get()
    {
      return item_.get();
    }
    const T* get() const
    {
      return item_.get();
    }
    // Return true if no data is held
    bool is_released() const
    {
      return (nullptr == item_.get());
    }
  };

  template <typename T>
  class OPENMS_DLLAPI unique_xerces_ptr
  {
  private:

    template<typename U>
    static void doRelease_(U*& item)
    {
      // Only release this if it has no parent (otherwise
      // parent will release it)
      if (nullptr == item->getOwnerDocument())
        item->release();
    }

    static void doRelease_(char*& item);
    static void doRelease_(char16_t*& item);

    T* item_;

  public:

    // Hide copy constructor and assignment operator
    unique_xerces_ptr(const unique_xerces_ptr<T>&) = delete;
    unique_xerces_ptr& operator=(const unique_xerces_ptr<T>&) = delete;

    unique_xerces_ptr()
        : item_(nullptr)
    {}

    explicit unique_xerces_ptr(T* i)
        : item_(i)
    {}

    ~unique_xerces_ptr()
    {
      xerces_release();
    }

    unique_xerces_ptr(unique_xerces_ptr<T>&& other) noexcept
        : item_(nullptr)
    {
      this->swap(other);
    }

    void swap(unique_xerces_ptr<T>& other) noexcept
    {
      std::swap(item_, other.item_);
    }

    // Assignment of data to guard (not chainable)
    void operator=(T* i)
    {
      reassign(i);
    }

    // Release held data (i.e. delete/free it)
    void xerces_release()
    {
      if (!is_released())
      {
        // Use type-specific release mechanism
        doRelease_(item_);
        item_ = nullptr;
      }
    }

    // Give up held data (i.e. return data without releasing)
    T* yield()
    {
      T* tempItem = item_;
      item_ = nullptr;
      return tempItem;
    }

    // Release currently held data, if any, to hold another
    void assign(T* i)
    {
      xerces_release();
      item_ = i;
    }

    // Get pointer to the currently held data, if any
    T* get() const
    {
      return item_;
    }

    // Return true if no data is held
    bool is_released() const
    {
      return (nullptr == item_);
    }
  };

  //========================================================================================================

  /**
   * @brief Helper class for XML parsing that handles the conversions of Xerces strings
   *
   * It provides the convert() function which internally calls
   * XMLString::transcode and ensures that the memory is released properly
   * through XMLString::release internally. It returns a std::string or
   * std::basic_string<char16_t> to the caller who takes ownership of the data.
   *
   * The public interface uses @c char16_t (UTF-16 code units), which is
   * guaranteed to match Xerces' @c XMLCh, so this header stays Xerces-free.
  */
  class OPENMS_DLLAPI StringManager
  {
    typedef std::basic_string<char16_t> XercesString;

    /// Converts from a narrow-character string to a wide-character string.
    static unique_xerces_ptr<char16_t> fromNative_(const char* str);

    /// Converts from a narrow-character string to a wide-character string.
    static unique_xerces_ptr<char16_t> fromNative_(const std::string& str);

    /// Converts from a wide-character string to a narrow-character string.
    static std::string toNative_(const char16_t* str);

    /// Converts from a wide-character string to a narrow-character string.
    static std::string toNative_(const unique_xerces_ptr<char16_t>& str);

protected:
    /// Compresses eight 8x16bit Chars in char16_t* to 8x8bit Chars by cutting upper byte
    static void compress64_ (const char16_t * input_it, char* output_it);

public:
    /// Constructor
    StringManager();

    /// Destructor
    ~StringManager();

    /// Calculates the length of a char16_t* string using SIMDe
    // https://github.com/OpenMS/OpenMS/issues/8122
    static Size strLength(const char16_t* input_ptr);

    /// Transcode the supplied C string to a UTF-16 string
    static XercesString convert(const char * str);

    /// Transcode the supplied C++ string to a UTF-16 string
    static XercesString convert(const std::string & str);

    /// Transcode the supplied C string to a UTF-16 string pointer
    static unique_xerces_ptr<char16_t> convertPtr(const char * str);

    /// Transcode the supplied C++ string to a UTF-16 string pointer
    static unique_xerces_ptr<char16_t> convertPtr(const std::string & str);

    /// Transcode the supplied char16_t* to a std::string
    static std::string convert(const char16_t * str);

    /// Parse an integer from a UTF-16 numeric string (as Xerces' XMLString::parseInt would).
    static Int parseInt(const char16_t * str);

    /// Checks if supplied chars in char16_t* can be encoded with ASCII (i.e. the upper byte of each char is 0)
    static bool isASCII(const char16_t * chars, const Size length);

    /**
     * @brief Transcodes the supplied char16_t* and appends it to the OpenMS String
     *
     * @note Assumes that the char16_t* only contains ASCII characters
     *
    */
    static void appendASCII(const char16_t * str, const Size length, std::string & result);
  };

} // namespace OpenMS::Internal
