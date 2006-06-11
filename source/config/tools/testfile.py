#!/usr/bin/env python

# script to search for common problems in BALL-files
# you can specify a specific problem to look for.
# to do so add a number as second argument.
# 0 is the first general problem
# 100 is the first problem in header-files

import sys, string, re

if len(sys.argv) <2:
	print 'zu wenige argumente'
	sys.exit(-1)

class test:
	errors=0
	line=''
	lastline=''
	linenr= 0 
	carriage=0
	testnr= 0
	f = open(sys.argv[1])

	# expressions for use in all files.
	exp = []
	# error messages
	msg = []

	exp +=[re.compile('\015')]
	msg +=['#00 DOS carriage']

	exp +=[re.compile('\?\?\?')]
	msg +=['#01 code problems']

	exp +=[re.compile('cout')]
	msg +=['#02 no cout in BALL!']

	exp +=[re.compile('cerr')]
	msg +=['#03 no cerr in BALL!']

	exp +=[re.compile('[\s(]+int[\s\(\)\&\*]+')]
	msg +=['#04 integer values are bad!']

	exp +=[re.compile('[\s(]+long[\s\(\)\&\*]+')]
	msg +=['#05 long values are']

	exp +=[re.compile('[\s(]+short[\s\(\)\&\*]+')]
	msg +=['#06 short values are bad!']

	exp +=[re.compile(';{2}')]
	msg +=['#07 ;;']

	exp +=[re.compile('}[\s]*else[\s]*{')]
	msg +=['#08 } else {']

	exp +=[re.compile('throw[\s]*Exception::NotImplemented')]
	msg +=['#09 Exception Not Implemented']

	exp +=[re.compile('\([\s]*Exception::NotImplemented')]
	msg +=['#10 Exception Not Implemented in throw specifier']

	exp +=[re.compile('[^:]std::endl')]
	msg +=['#11 missing :: before std::endl']

	exp +=[re.compile('\([\s]*bool[\s]*\)')]
	msg +=['#12 superflous bool cast']

	exp +=[re.compile('const[\s]*float[\s]*&')]
	msg +=['#13 no const float references']

	exp +=[re.compile('const[\s]*double[\s]*&')]
	msg +=['#14 no const double references']

	exp += [re.compile('^M')]
	msg +=['#15 msdos carriage return']

	#	99 tab info line missing

	
	# expressions for use with header-files
	exp_header = []
	msg_header = []

	exp_header +=[re.compile('///[\s]*\Z')]
	msg_header +=['#100 empty comment']

	exp_header +=[re.compile('/\*\*[\s]*\Z')]
	msg_header +=['#101 empty comment']

	exp_header +=[re.compile('@exception[\s]*NotImplemented')]
	msg_header +=['#102 no usefull information']

	exp_header +=[re.compile('@param[\s]*{')]
	msg_header +=['#103 standard problem => tex error']

	exp_header +=[re.compile('@return[\s]*{')]
	msg_header +=['#104 standard problem => tex error']
	
	
	#test if check is empty																		
	exp_check = [	
		re.compile('RESULT'),
		re.compile('\A[\s]*\Z'),
		re.compile('BAUSTELLE')
	]

		
	# expressions for use in test-files
	#exp_test = []																				#not yet an idea

			
	def write(self, error_code):
		print '\n--------------------------------- ' + `self.linenr` + ' --------------- ' + error_code
		print string.strip(self.line[0:-1]),


	def ende(self):
		self.f.close()
		if self.errors > 0:
			print
		sys.exit(self.errors)


	def test_line(self, x):
		#test for DOS-carriage, print it just once
		if self.exp[0].search(x, 0) and self.carriage==0:
			self.carriage=1
			self.write(self.msg[0])
		for i in range(0, len(self.exp)):
			if self.exp[i].search(x, 0):
				self.errors = self.errors + 1	
				self.write(self.msg[i])
	
	
	def getLine(self):
		self.linenr = self.linenr + 1
		self.lastline = self.line
		self.line = self.f.readline()
		if not self.line: return 0
		self.test_line(self.line)
		return 1

	
	#test if file start with tab info lines
	def BALL_TABINFO_TEST(self):
		tabInfoLines=(
			'// -*- Mode: C++; tab-width: 2; -*-',
			'// vi: set ts=2:',
			'//')
		#compare first 3 lines with tab info lines
		for reg in tabInfoLines:
			self.linenr = self.linenr + 1
			self.lastline = self.line
			self.line = self.f.readline()
			if not self.line: break;
			if string.find(self.line, reg) == -1: 
				self.write('#99 tab info lines missing')
				break
		return 1


	def BALL_TEST_TEST(self):
		while self.getLine() == 1:
			if self.line[:5] =='CHECK':
				self.getLine()
				for i in range(len(self.exp_check)):
					if self.exp_check[i].search(self.line, 0):
						self.errors = self.errors + 1	
						self.write('test ' + `i`)		
		self.ende()


	def BALL_HEADER_TEST(self):
		while self.getLine() == 1:
			for i in range(len(self.exp_header)):
				if self.exp_header[i].search(self.line, 0):
					self.errors = self.errors + 1	
					self.write(self.msg_header[i])
		self.ende()


	def BALL_ALL_TEST(self):
		while self.getLine() == 1:
			pass
		self.ende()

	def ONE_TEST(self):
		if self.testnr < 100:
			if self.testnr > len(self.exp):
				print 'no valid number of test given, aborting...'
				sys.exit(-1)
			myexp=self.exp[self.testnr]
		if self.testnr > 99 and self.testnr < 200:
			if self.testnr - 100 > len(self.exp_header):
				print 'no valid number of test given, aborting...'
				sys.exit(-1)
			myexp=self.exp_header[self.testnr - 100]
		
		self.line = self.f.readline()
		while self.line:
			self.linenr = self.linenr + 1
			if myexp.search(self.line, 0):
				self.errors = self.errors + 1
				self.write(sys.argv[2])
			self.lastline = self.line
			self.line = self.f.readline()
		self.ende()

	def debug(self):
		while self.getLine() == 1:
			# Header tests
			for i in range(len(self.exp_header)):
				if self.exp_header[i].search(self.line, 0):
					self.errors = self.errors + 1
					print '  H ', i
					
			# Testfile tests	
			for i in range(len(self.exp_test)):
				if self.exp_test[i].search(self.line, 0):
					self.errors = self.errors + 1
					print '  T ', i


	def __init__(self):
		try:
			self.testnr=int(sys.argv[2])
			if self.testnr>299 or self.testnr<0:
				print 'no valid number of test given, aborting...'
				self.testnr=0
		except:
			self.testnr=0
		# test for a given problem
		if self.testnr != 0:
			self.ONE_TEST()
			sys.exit(self.errors)
		
		# file for debuging this script
		if string.find(sys.argv[1], 'test.file') != -1:
			self.debug()
			sys.exit(self.errors)
			
		# tab info lines test
		self.BALL_TABINFO_TEST()

		# general testing
		self.BALL_ALL_TEST()
		
		# test for tabulator config files
		
		
		# file is a TEST-FILE
		if string.find(sys.argv[1], '_test.C') != -1:
			self.BALL_TEST_TEST()
			sys.exit(self.errors)
			
		# file is a HEADER-FILE
		if string.find(sys.argv[1], '.h') != -1:
			self.BALL_HEADER_TEST()
			sys.exit(self.errors)

		sys.exit(self.errors)



t = test()
sys.exit(t.errors)
