// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_BENCHMARK_H
#define OPENMS_CONCEPT_BENCHMARK_H

#include <OpenMS/config.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <iostream>

/** 
	@defgroup Benchmark Benchmarking macros

  @brief Macros used by the benchmark.
	
	@ingroup Concept
	
	@{
*/

/**
	Start a new benchmark section.
	The argument weight determines the weighting factor of the section.
	
	@hideinitializer
*/
#define START_SECTION(name, weight) \
	BENCHMARK::section_time = BENCHMARK::timer.getCPUTime();\
	BENCHMARK::section_name = #name;\
	BENCHMARK::section_weight = weight;


/**
	End of a benchmark section.
	
	@hideinitializer
*/
#define END_SECTION \
	BENCHMARK::timer.stop();\
	BENCHMARK::section_time = BENCHMARK::timer.getCPUTime() - BENCHMARK::section_time;\
	if (BENCHMARK::verbose > 0)\
	{\
		std::cout << BENCHMARK::section_name << ": " \
		  << BENCHMARK::section_time << " s"\
			<< " (weight = " << BENCHMARK::section_weight << ")" << std::endl;\
	}\
	BENCHMARK::total_time += BENCHMARK::section_time * BENCHMARK::section_weight;\


/**	
	@brief Status output.
	
	Print debugging information if called with -v.

	@hideinitializer
*/
#define STATUS(a) \
	if (BENCHMARK::verbose > 0)\
	{\
		std::cout << "  status: " << a << std::endl;\
	}


/**	
	@brief Start the timer.
	
	This macro is used to determine the running time
	of a set of commands.  It may be used in benchmarks and requires
	a prior invocation of the START_BENCHMARK() macro.  All commands
	that are between the START_TIMER() and the STOP_TIMER() command
	contribute to the overall running time of the benchmark.

	@hideinitializer
*/
#define START_TIMER \
	BENCHMARK::timer.start();\


/**	
	@brief Stop the timer.

	This macro is used to determine the running time of a set of
	commands.  It may be used in benchmarks and requires a prior
	invocation of the START_BENCHMARK() and START_TIMER() macros.  All
	commands that are between the START_TIMER() and the STOP_TIMER()
	command contribute to the overall running time of the benchmark.

	@hideinitializer
*/
#define STOP_TIMER \
	BENCHMARK::timer.stop();

/**	
	@brief Program body for the benchmark.
		
	The parameter <tt>weight</tt> determines the overall weight of
	this test in the accumulated benchmark (OpenMSStones).

	@hideinitializer
*/
#define START_BENCHMARK(class_name, overall_weight, version)\
/* define a special namespace for all internal variables */\
/* to avoid potential collisions                         */\
namespace BENCHMARK {\
	int						verbose = 0;\
	bool					all_tests = true;\
	int						exception = 0;\
	std::string		exception_name = "";\
	const char*		version_string = version;\
	std::string		section_name = "";\
	float					section_weight = 1.0;\
	float					weight = overall_weight;\
	float					total_time;\
	float					section_time;\
	OpenMS::StopWatch		timer;\
}\
\
\
int main(int argc, char **argv)\
{\
\
	if (argc == 2) {\
		if (!strcmp(argv[1], "-v"))\
			BENCHMARK::verbose = 1;\
	};\
\
	if ((argc > 2) || ((argc == 2) && (BENCHMARK::verbose == 0))) {\
		std::cerr << "Execute a benchmark for the " #class_name " class." << std::endl;\
		std::cerr << "Overall weight of the test: " << BENCHMARK::weight << std::endl;\
\
		std::cerr << "On successful operation, the total CPU time (in seconds)," << std::endl;\
		std::cerr << "is printed." << std::endl;\
		std::cerr << "If called with an argument of -v, " << argv[0] << " detailed" << std::endl;\
		std::cerr << "information about individual benchmarks is printed." << std::endl;\
		return 1;\
	}\
\
	if (BENCHMARK::verbose > 0)\
		std::cout << "Version: " << BENCHMARK::version_string << std::endl;\
\
	try {\

/**	
	@brief End of the test program

	@hideinitializer	
*/
#define END_BENCHMARK \
	/* global try block */\
	}\
	/* catch FileNotFound exceptions to print out the file name */\
	catch (OpenMS::Exception::FileNotFound& e)\
	{\
		BENCHMARK::all_tests = false;\
  	if (BENCHMARK::verbose > 1)\
		{\
			if (BENCHMARK::exception == 1) /* dummy to avoid compiler warnings */\
				BENCHMARK::exception++;\
    	std::cout << std::endl << "    (caught exception of type ";\
			std::cout << e.getName();\
			if ((e.getLine() > 0) && (!(e.getFile() == "")))\
				std::cout << " outside a benchmark block, which was thrown in line " << e.getLine() << " of file " << e.getFile();\
			std::cout << " while looking for file " << e.getFilename();\
			std::cout << " - unexpected!) " << std::endl;\
		}\
  }\
	/* catch OpenMS exceptions to retrieve additional information */\
	catch (OpenMS::Exception::Base& e)\
	{\
		BENCHMARK::all_tests = false;\
  	if (BENCHMARK::verbose > 1)\
		{\
			if (BENCHMARK::exception == 1) /* dummy to avoid compiler warnings */\
				BENCHMARK::exception++;\
    	std::cout << std::endl << "    (caught exception of type ";\
			std::cout << e.getName();\
			if ((e.getLine() > 0) && (!(e.getFile() == "")))\
				std::cout << " outside a benchmark block, which was thrown in line " << e.getLine() << " of file " << e.getFile();\
			std::cout << " - unexpected!) " << std::endl;\
		}\
  }\
	/* catch all non-OpenMS exceptions */\
	catch (...)\
	{\
		BENCHMARK::all_tests = false;\
  	if (BENCHMARK::verbose > 1)\
		{\
    	std::cout << std::endl << "    (caught unidentified and unexpected exception outside a benchmark block!) " << std::endl;\
		}\
	}\
\
	/* check for exit code */\
	if (!BENCHMARK::all_tests)\
	{\
		std::cout << "(" << BENCHMARK::weight * BENCHMARK::total_time << ")" << std::endl;\
		return 1;\
	} else {\
		std::cout << BENCHMARK::weight * BENCHMARK::total_time << std::endl;\
		return 0;\
	}\
}\

/** @} */ // end of Benchmark

#endif // OPENMS_CONCEPT_BENCHMARK_H
