// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Roest $
// $Authors: Philippe M. Groarke, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once
#include <type_traits>

#define OPENMS_TEST_DEFAULT_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_default_constructible_v<t>, \
				#t " : must be default constructible"); \
		return std::is_default_constructible_v<t>; \
	}

#define OPENMS_TEST_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_default_constructible_v<t>, \
				#t " : must be trivially default constructible"); \
		return std::std::is_trivially_default_constructible_v<t>; \
	}

#define OPENMS_TEST_COPY_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_copy_constructible_v<t>, \
				#t " : must be copy constructible"); \
		return std::is_copy_constructible_v<t>; \
	}

#define OPENMS_TEST_NOT_COPY_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(!std::is_copy_constructible_v<t>, \
				#t " : must not be copy constructible"); \
		return !std::is_copy_constructible_v<t>; \
	}

#define OPENMS_TEST_TRIVIALLY_COPY_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_copy_constructible_v<t>, \
				#t " : must be trivially copy constructible"); \
		return std::is_trivially_copy_constructible_v<t>; \
	}

#define OPENMS_TEST_COPY_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_copy_assignable_v<t>, \
				#t " : must be copy assignable"); \
		return std::is_copy_assignable_v<t>; \
	}

#define OPENMS_TEST_NOT_COPY_ASSIGNABLE(t) \
	[]() { \
		static_assert(!std::is_copy_assignable_v<t>, \
				#t " : must not be copy assignable"); \
		return !std::is_copy_assignable_v<t>; \
	}

#define OPENMS_TEST_TRIVIALLY_COPY_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_trivially_copy_assignable_v<t>, \
				#t " : must be trivially copy assignable"); \
		return std::is_trivially_copy_assignable_v<t>; \
	}

#define OPENMS_TEST_MOVE_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_move_constructible_v<t>, \
				#t " : must be move constructible"); \
		return std::is_move_constructible_v<t>; \
	}

#define OPENMS_TEST_TRIVIALLY_MOVE_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_move_constructible_v<t>, \
				#t " : must be trivially move constructible"); \
		return std::is_trivially_move_constructible_v<t>; \
	}

#define OPENMS_TEST_NOTHROW_MOVE_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_nothrow_move_constructible_v<t>, \
				#t " : must be nothrow move constructible"); \
		return std::is_nothrow_move_constructible_v<t>; \
	}

#define OPENMS_TEST_MOVE_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_move_assignable_v<t>, \
				#t " : must be move assignable"); \
		return std::is_move_assignable_v<t>; \
	}

#define OPENMS_TEST_TRIVIALLY_MOVE_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_trivially_move_assignable_v<t>, \
				#t " : must be trivially move assignable"); \
		return std::is_trivially_move_assignable_v<t>; \
	}

#define OPENMS_TEST_DESTRUCTIBLE(t) \
	[]() { \
		static_assert( \
				std::is_destructible_v<t>, #t " : must be destructible"); \
		return std::is_destructible_v<t>; \
	}

#define OPENMS_TEST_TRIVIALLY_DESTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_destructible_v<t>, \
				#t " : must be trivially destructible"); \
		return std::is_trivially_destructible_v<t>; \
	}

#define OPENMS_TEST_TRIVIALLY_COPYABLE(t) \
	[]() { \
		static_assert(std::is_trivially_copyable_v<t>, \
				#t " : must be trivially copyable"); \
		return std::is_trivially_copyable_v<t>; \
	}

#define OPENMS_TEST_FULFILLS_RULE_OF_5(t) \
	static_assert(OPENMS_TEST_DESTRUCTIBLE(t)() && OPENMS_TEST_COPY_CONSTRUCTIBLE(t)() \
					&& OPENMS_TEST_MOVE_CONSTRUCTIBLE(t)() && OPENMS_TEST_COPY_ASSIGNABLE(t)() \
					&& OPENMS_TEST_MOVE_ASSIGNABLE(t)(), \
			#t " : doesn't fulfill rule of 5")

#define OPENMS_TEST_FULFILLS_RULE_OF_6(t) \
	static_assert(OPENMS_TEST_DESTRUCTIBLE(t)() && OPENMS_TEST_DEFAULT_CONSTRUCTIBLE(t)() \
					&& OPENMS_TEST_COPY_CONSTRUCTIBLE(t)() \
					&& OPENMS_TEST_MOVE_CONSTRUCTIBLE(t)() && OPENMS_TEST_COPY_ASSIGNABLE(t)() \
					&& OPENMS_TEST_MOVE_ASSIGNABLE(t)(), \
			#t " : doesn't fulfill rule of 5")

// is_trivially_copyable broken everywhere
#define OPENMS_TEST_FULFILLS_FAST_VECTOR(t) \
	static_assert((std::is_trivially_copy_constructible_v< \
						   t> && std::is_trivially_destructible_v<t>) \
					|| std::is_nothrow_move_constructible_v<t>, \
			#t " : doesn't fulfill fast vector (trivially copy constructible " \
			   "and trivially destructible, or nothrow move constructible)")

#define OPENMS_TEST_FULFILLS_MOVE_ONLY(t) \
	static_assert(OPENMS_TEST_NOT_COPY_CONSTRUCTIBLE(t)() \
					&& OPENMS_TEST_MOVE_CONSTRUCTIBLE(t)() \
					&& OPENMS_TEST_NOT_COPY_ASSIGNABLE(t)() \
					&& OPENMS_TEST_MOVE_ASSIGNABLE(t)(), \
			#t " : doesn't fulfill move only")

