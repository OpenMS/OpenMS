#!/usr/bin/env ruby

# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

require 'getoptlong'
require 'erb'

def getTemplatePath(file)
	script_dir = File.dirname(__FILE__)
	return "#{script_dir}/templates/#{file}"
end

def getTemplate(file)
	return File.read(getTemplatePath(file))
end

def isOpenMS(openms_path)
  isOpenMS = TRUE
  if not File.directory?(openms_path)
    puts "Given path is not a directory. Abort!"
    isOpenMS = FALSE
  elsif not File.exist?("#{openms_path}/CMakeLists.txt")
    puts "Given path is not the OpenMS root. Abort!"
    isOpenMS = FALSE
  elsif not File.exist?("#{openms_path}/src/openms/source/")
    puts "Given path is not the OpenMS root. Abort!"
    isOpenMS = FALSE
  elsif not File.exist?("#{openms_path}/src/openms/include/OpenMS/")
    puts "Given path is not the OpenMS root. Abort!"
    isOpenMS = FALSE
  end

  return isOpenMS
end

def libExists(openms_path, lib_name)
	libExists = TRUE

	if not File.directory?("#{openms_path}/src/#{lib_name}")
		puts "The given lib #{lib_name} doesn't exist. Abort!"
		libExists = FALSE
	elsif not File.exist?("#{openms_path}/src/#{lib_name}/CMakeLists.txt")
		puts "The given lib #{lib_name} is not a valid OpenMS sub-lib. Abort!"
		libExists = FALSE
	end

	return libExists
end

def getLicense(openms_path)
	license = File.read("#{openms_path}/LICENSE")
	return license.gsub(/^/, "// ").chop()
end

def create_header_sources(openms_path, lib_name, path)
	# create directory
	new_dir = "#{openms_path}/src/#{lib_name}/include/OpenMS/#{path}"
	Dir.mkdir(new_dir)

	# template variables
	header_directory = path
	header_escaped_path = path.gsub("/", "\\\\\\\\")
	template = ERB.new(getTemplate("header_sources_cmake"))

	doc = template.result(binding)
	sources_cmake = "#{new_dir}/sources.cmake"
	File.open(sources_cmake, 'w') {|f| f.write(doc) }

	puts "Added new directory to the OpenMS build system. Please register manually in cmake/includes.cmake by adding the following line"
	puts "include(#{sources_cmake.sub(openms_path, "").sub("/","")})"
end

def create_source_sources(openms_path, lib_name, path)
	# create directory
	new_dir = "#{openms_path}/src/#{lib_name}/source/#{path}"
	Dir.mkdir(new_dir)

	# template variables
	source_directory = path
	source_escaped_path= path.gsub("/", "\\\\\\\\")
	template = ERB.new(getTemplate("source_sources_cmake"))

	doc = template.result(binding)
	sources_cmake = "#{new_dir}/sources.cmake"
	File.open(sources_cmake, 'w') {|f| f.write(doc) }

	puts "Added new directory to the OpenMS build system. Please register manually in cmake/includes.cmake by adding the following line"
	puts "include(#{sources_cmake.sub(openms_path, "").sub("/","")})"
end

def create_sources(openms_path, lib_name, path, clazz, maintainer)

	license = getLicense(openms_path)

	# write the header file
	header_guard = path.gsub("/", "_") + "_" + clazz.upcase + "_H"
	template = ERB.new(getTemplate("header"))

	doc = template.result(binding)
	header_file = "#{openms_path}/src/#{lib_name}/include/OpenMS/#{path}/#{clazz}.h"
	File.open(header_file, 'w') {|f| f.write(doc) }

	# write the source file
	header = "#{path}/#{clazz}.h"
	template = ERB.new(getTemplate("source"))

	doc = template.result(binding)
	source_file = "#{openms_path}/src/#{lib_name}/source/#{path}/#{clazz}.cpp"
	File.open(source_file, 'w') {|f| f.write(doc) }
end

def register_file(sources_cmake, filename)
	contentsArray = File.readlines(sources_cmake)

	block = false
	insert_idx = -1
	contentsArray.each_with_index {|val, index|
		if block
			if val.chop().strip == ")"
				insert_idx = index
				break
			end
		else
			if val.chop().strip == "set(sources_list" or val.chop().strip == "set(sources_list_h"
				block = true
			end
		end
	}

	contentsArray.insert(insert_idx, filename)
	File.open(sources_cmake, 'w') do |file|
	  file.puts contentsArray
	end
end

def register(openms_path, lib_name, path, clazz)
	include_cmake = "#{openms_path}/src/#{lib_name}/include/OpenMS/#{path}/sources.cmake"
	source_cmake = "#{openms_path}/src/#{lib_name}/source/#{path}/sources.cmake"

	register_file(include_cmake, "#{clazz}.h")
	register_file(source_cmake, "#{clazz}.cpp")
end

##### INSTANTIATE CMD PARSER
opts = GetoptLong.new(
  [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
  [ '--maintainer', '-v', GetoptLong::OPTIONAL_ARGUMENT ]
)

maintainer = ""

usage = "#{File.basename($0)} --maintainer \"Maintainer Line\" OPENMS_SOURCE LIBNAME CLAZZNAME_W_PATH

-h, --help:
  show help
-v, --verbose:
  increase verbosity
--maintainer:
	set the maintainer of the generated class
"

opts.each do |opt, arg|
  case opt
    when '--help'
      puts usage.to_s
    when '--maintainer'
      maintainer = arg
    when '-v'
      $DEBUG = true
  end
end

if ARGV.length != 3
	puts "Missing openms path, lib name and class name. Abort!"
	exit 1
elsif maintainer == ""
  puts "Missing maintainer/author information. Abort!"
  exit 1
end

openms_path = ARGV.shift
lib_name = ARGV.shift
clazz_w_path = ARGV.shift

if not isOpenMS(openms_path)
  exit 1
end

if not libExists(openms_path, lib_name)
	exit 1
end

# check if the required path already exists
path_elements=clazz_w_path.split("/")

# extract path/clazz name
clazz=path_elements[path_elements.length-1]
path=path_elements[0..path_elements.length-2] * "/"

# create missing directories
if not File.exist?("#{openms_path}/src/#{lib_name}/include/OpenMS/#{path}")
	create_header_sources(openms_path, lib_name, path)
end

if not File.exist?("#{openms_path}/src/#{lib_name}/source/#{path}")
	create_source_sources(openms_path, lib_name, path)
end

# create header / source
create_sources(openms_path, lib_name, path, clazz, maintainer)

# register in sources.cmake
register(openms_path, lib_name, path, clazz)

puts "Successfully added class \"#{clazz}\" to OpenMS. Do not forget to add the corresponding test."
