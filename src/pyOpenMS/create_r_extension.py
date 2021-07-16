# input-encoding: latin-1
from __future__ import print_function

# windows ?
import sys
iswin = sys.platform == "win32"

# make sure we only log errors and not info/debug ...
from logging import basicConfig
from logging import CRITICAL, ERROR, WARNING, INFO, DEBUG
basicConfig(level=21)

import config
from env import (QT_QMAKE_VERSION_INFO, OPEN_MS_BUILD_TYPE, OPEN_MS_SRC,
                 OPEN_MS_CONTRIB_BUILD_DIRS, OPEN_MS_LIB, OPEN_SWATH_ALGO_LIB, SUPERHIRN_LIB,
                 OPEN_MS_BUILD_DIR, MSVS_RTLIBS, OPEN_MS_VERSION,
                 Boost_MAJOR_VERSION, Boost_MINOR_VERSION)

IS_DEBUG = OPEN_MS_BUILD_TYPE.upper() == "DEBUG"

if iswin and IS_DEBUG:
    raise Exception("building pyopenms on windows in debug mode not tested yet.")

# use autowrap to generate Cython and .cpp file for wrapping OpenMS:
import autowrap.Main
import autowrap.CodeGenerator
import autowrap.DeclResolver
import glob
import pickle
import os.path
import os
import shutil

classdocu_base = "http://www.openms.de/current_doxygen/html/"
autowrap.CodeGenerator.special_class_doc = "Documentation is available at " + classdocu_base + "class%(namespace)s_1_1%(cpp_name)s.html\n"
autowrap.DeclResolver.default_namespace = "OpenMS"

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while len(out) < num:
      out.append(seq[int(last):int(last + avg)])
      last += avg

    # append the rest to the last element (if there is any)
    out[-1].extend( seq[int(last):] )
    return out

j = os.path.join

# Threads and Modules
PY_NUM_THREADS = 1
PY_NUM_MODULES = 1

src_pyopenms = j(OPEN_MS_SRC, "src/pyOpenMS")
pxd_files = glob.glob(src_pyopenms + "/pxds/*.pxd")
addons = glob.glob(src_pyopenms + "/R-addons/*.R")
converters = [j(src_pyopenms, "R-converters")]

persisted_data_path = "include_dir.bin"

extra_cimports = []

# We need to parse them all together but keep the association about which class
# we found in which file (as they often need to be analyzed together)
decls, instance_map = autowrap.parse(pxd_files, ".", num_processes=int(PY_NUM_THREADS))

# Perform mapping
pxd_decl_mapping = {}
for de in decls:
    tmp = pxd_decl_mapping.get(de.cpp_decl.pxd_path, [])
    tmp.append(de)
    pxd_decl_mapping[ de.cpp_decl.pxd_path] = tmp


# Split into chunks based on pxd files and store the mapping to decls, addons
# and actual pxd files in a hash. We need to produce the exact number of chunks
# as setup.py relies on it as well.
pxd_files_chunk = chunkIt(list(pxd_decl_mapping.keys()), int(PY_NUM_MODULES))

# Sanity checks: we should find all of our chunks and not have lost files
if len(pxd_files_chunk) != int(PY_NUM_MODULES):
    raise Exception("Internal Error: number of chunks not equal to number of modules")
if sum([len(ch) for ch in pxd_files_chunk]) != len(pxd_decl_mapping):
    raise Exception("Internal Error: chunking lost files")

mnames = ["pyopenms_%s" % (k+1) for k in range(int(PY_NUM_MODULES))]
allDecl_mapping = {}
for pxd_f, m in zip(pxd_files_chunk, mnames):
    tmp_decls = []
    for f in pxd_f:
        tmp_decls.extend( pxd_decl_mapping[f] )

    allDecl_mapping[m] =  {"decls" : tmp_decls, "addons" : [] , "files" : pxd_f}


# Deal with addons, make sure the addons are added to the correct class generator.
is_added = [False for k in addons]
for modname in mnames:

    for k,a in enumerate(addons):
        # Match addon basename to pxd basename
        for pfile in allDecl_mapping[modname]["files"]:
            if os.path.basename(a).split(".")[0] == os.path.basename(pfile).split(".")[0]:
                allDecl_mapping[modname]["addons"].append(a)
                is_added[k] = True

        if is_added[k]:
            continue

        # Also match by class name (sometimes one pxd contains multiple classes
        # and the addon is named after one of them)
        for dclass in allDecl_mapping[modname]["decls"]:
            if os.path.basename(a) == dclass.name + ".R":
                allDecl_mapping[modname]["addons"].append(a)
                is_added[k] = True

# add any addons that did not get added anywhere else
for k, got_added in enumerate(is_added):
    if not got_added:
        # add to all modules
        for m in mnames:
            allDecl_mapping[m]["addons"].append( addons[k] )

def doRCodeGeneration(modname, allDecl_mapping, instance_map, converters):
    # create the ropenms/ and R/ directory.
    src_code_dir = j(OPEN_MS_SRC,"src/ropenms/R")

    try:
        os.makedirs(src_code_dir)
    except OSError:
        pass
    
    m_filename = j(src_code_dir,"%s.R" % modname)
    cimports, manual_code = autowrap.Main.collect_manual_code(allDecl_mapping[modname]["addons"])
    autowrap.Main.register_converters(converters)
    # Need to specify binding "R" as default is "cython"
    autowrap_include_dirs = autowrap.generate_code(allDecl_mapping[modname]["decls"], instance_map,
                                                        target=m_filename, debug=False, manual_code=manual_code,
                                                        extra_cimports=cimports,
                                                        include_boost=False, include_numpy=True, allDecl=allDecl_mapping, target_binding="R")
    allDecl_mapping[modname]["inc_dirs"] = autowrap_include_dirs

    # 
    # copy the tests folder from ropenms-attachments to ropenms
    #
    dest = j(OPEN_MS_SRC,"src/ropenms/tests")

    if os.path.isdir(dest):
        shutil.rmtree(dest, ignore_errors=True)
    
    src = j(OPEN_MS_SRC,"src/pyOpenMS/ropenms-attachments/tests")
    shutil.copytree(src,dest)
    return autowrap_include_dirs


for modname in mnames:
    autowrap_include_dirs = doRCodeGeneration(modname, allDecl_mapping, instance_map, converters)
    pickle.dump(autowrap_include_dirs, open(persisted_data_path, "wb"))

argzip = [ (modname, allDecl_mapping[modname]["inc_dirs"]) for modname in mnames]

print("created ropenms package")