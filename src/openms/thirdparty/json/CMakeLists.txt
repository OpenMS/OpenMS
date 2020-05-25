# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS-json")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)


set(_JSON_SOURCES
adl_serializer.hpp
detail
detail/conversions
detail/conversions/from_json.hpp
detail/conversions/to_chars.hpp
detail/conversions/to_json.hpp
detail/exceptions.hpp
detail/input
detail/input/binary_reader.hpp
detail/input/input_adapters.hpp
detail/input/json_sax.hpp
detail/input/lexer.hpp
detail/input/parser.hpp
detail/iterators
detail/iterators/internal_iterator.hpp
detail/iterators/iteration_proxy.hpp
detail/iterators/iter_impl.hpp
detail/iterators/json_reverse_iterator.hpp
detail/iterators/primitive_iterator.hpp
detail/json_pointer.hpp
detail/json_ref.hpp
detail/macro_scope.hpp
detail/macro_unscope.hpp
detail/meta
detail/meta/cpp_future.hpp
detail/meta/detected.hpp
detail/meta/is_sax.hpp
detail/meta/type_traits.hpp
detail/meta/void_t.hpp
detail/output
detail/output/binary_writer.hpp
detail/output/output_adapters.hpp
detail/output/serializer.hpp
detail/value_t.hpp
json.hpp
json_fwd.hpp
)

# export include directory
set(JSON_INCLUDE_DIRECTORY
  ${PROJECT_SOURCE_DIR}
  CACHE
  INTERNAL "json include directories" FORCE
)

set(JSON_SOURCES )
foreach(_file ${_JSON_SOURCES})
  # append source files
  set(
    JSON_SOURCES
    ${JSON_SOURCES}
    ${PROJECT_SOURCE_DIR}/nlohmann/${_file}
  )
endforeach()

