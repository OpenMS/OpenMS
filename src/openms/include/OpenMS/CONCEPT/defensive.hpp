/*
			DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
					Version 2, December 2004

 Copyright (C) 2018 Philippe Groarke <philippe.groarke@gmail.com>

 Everyone is permitted to copy and distribute verbatim or modified
 copies of this license document, and changing it is allowed as long
 as the name is changed.

			DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. You just DO WHAT THE FUCK YOU WANT TO.
*/
#pragma once
#include <type_traits>

#define FEA_DEFAULT_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_default_constructible_v<t>, \
				#t " : must be default constructible"); \
		return std::is_default_constructible_v<t>; \
	}

#define FEA_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_default_constructible_v<t>, \
				#t " : must be trivially default constructible"); \
		return std::std::is_trivially_default_constructible_v<t>; \
	}

#define FEA_COPY_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_copy_constructible_v<t>, \
				#t " : must be copy constructible"); \
		return std::is_copy_constructible_v<t>; \
	}

#define FEA_NOT_COPY_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(!std::is_copy_constructible_v<t>, \
				#t " : must not be copy constructible"); \
		return !std::is_copy_constructible_v<t>; \
	}

#define FEA_TRIVIALLY_COPY_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_copy_constructible_v<t>, \
				#t " : must be trivially copy constructible"); \
		return std::is_trivially_copy_constructible_v<t>; \
	}

#define FEA_COPY_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_copy_assignable_v<t>, \
				#t " : must be copy assignable"); \
		return std::is_copy_assignable_v<t>; \
	}

#define FEA_NOT_COPY_ASSIGNABLE(t) \
	[]() { \
		static_assert(!std::is_copy_assignable_v<t>, \
				#t " : must not be copy assignable"); \
		return !std::is_copy_assignable_v<t>; \
	}

#define FEA_TRIVIALLY_COPY_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_trivially_copy_assignable_v<t>, \
				#t " : must be trivially copy assignable"); \
		return std::is_trivially_copy_assignable_v<t>; \
	}

#define FEA_MOVE_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_move_constructible_v<t>, \
				#t " : must be move constructible"); \
		return std::is_move_constructible_v<t>; \
	}

#define FEA_TRIVIALLY_MOVE_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_move_constructible_v<t>, \
				#t " : must be trivially move constructible"); \
		return std::is_trivially_move_constructible_v<t>; \
	}

#define FEA_NOTHROW_MOVE_CONSTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_nothrow_move_constructible_v<t>, \
				#t " : must be nothrow move constructible"); \
		return std::is_nothrow_move_constructible_v<t>; \
	}

#define FEA_MOVE_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_move_assignable_v<t>, \
				#t " : must be move assignable"); \
		return std::is_move_assignable_v<t>; \
	}

#define FEA_TRIVIALLY_MOVE_ASSIGNABLE(t) \
	[]() { \
		static_assert(std::is_trivially_move_assignable_v<t>, \
				#t " : must be trivially move assignable"); \
		return std::is_trivially_move_assignable_v<t>; \
	}

#define FEA_DESTRUCTIBLE(t) \
	[]() { \
		static_assert( \
				std::is_destructible_v<t>, #t " : must be destructible"); \
		return std::is_destructible_v<t>; \
	}

#define FEA_TRIVIALLY_DESTRUCTIBLE(t) \
	[]() { \
		static_assert(std::is_trivially_destructible_v<t>, \
				#t " : must be trivially destructible"); \
		return std::is_trivially_destructible_v<t>; \
	}

#define FEA_TRIVIALLY_COPYABLE(t) \
	[]() { \
		static_assert(std::is_trivially_copyable_v<t>, \
				#t " : must be trivially copyable"); \
		return std::is_trivially_copyable_v<t>; \
	}

#define FEA_FULFILLS_RULE_OF_5(t) \
	static_assert(FEA_DESTRUCTIBLE(t)() && FEA_COPY_CONSTRUCTIBLE(t)() \
					&& FEA_MOVE_CONSTRUCTIBLE(t)() && FEA_COPY_ASSIGNABLE(t)() \
					&& FEA_MOVE_ASSIGNABLE(t)(), \
			#t " : doesn't fulfill rule of 5")

#define FEA_FULFILLS_RULE_OF_6(t) \
	static_assert(FEA_DESTRUCTIBLE(t)() && FEA_DEFAULT_CONSTRUCTIBLE(t)() \
					&& FEA_COPY_CONSTRUCTIBLE(t)() \
					&& FEA_MOVE_CONSTRUCTIBLE(t)() && FEA_COPY_ASSIGNABLE(t)() \
					&& FEA_MOVE_ASSIGNABLE(t)(), \
			#t " : doesn't fulfill rule of 5")

// is_trivially_copyable broken everywhere
#define FEA_FULFILLS_FAST_VECTOR(t) \
	static_assert((std::is_trivially_copy_constructible_v< \
						   t> && std::is_trivially_destructible_v<t>) \
					|| std::is_nothrow_move_constructible_v<t>, \
			#t " : doesn't fulfill fast vector (trivially copy constructible " \
			   "and trivially destructible, or nothrow move constructible)")

#define FEA_FULFILLS_MOVE_ONLY(t) \
	static_assert(FEA_NOT_COPY_CONSTRUCTIBLE(t)() \
					&& FEA_MOVE_CONSTRUCTIBLE(t)() \
					&& FEA_NOT_COPY_ASSIGNABLE(t)() \
					&& FEA_MOVE_ASSIGNABLE(t)(), \
			#t " : doesn't fulfill move only")
