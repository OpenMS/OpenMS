import numpy
import os
import glob
import pyopenms

toTest = set()

# pyopenms in the ignore list as it is represented twice in
# the pyopenms package
ignore = ["numpy", "np", "re", "os", "types", "sysinfo", "pyopenms"]

for clz_name, clz in pyopenms.__dict__.items():
    if clz_name in ignore or clz_name.startswith("__"):
        continue

    if not hasattr(clz, "__dict__"):
        continue
    for method_name, method in clz.__dict__.items():

        if method_name.startswith("_") and not method_name.startswith("__"):
            continue

        if method_name in [ "__doc__", "__new__", "__file__", "__name__",
                "__package__", "__builtins__", "__copy__"]:
            continue

        toTest.add("%s.%s" % (clz_name, method_name))

def parse_doc(item, collection):
    if item.__doc__ is not None:
       it = iter(item.__doc__.split("\n"))
       for line in it:
            if not "@tests" in line:
                continue
            for line in it:
                line = line.strip()
                if "@end" in line:
                    break
                if not line:
                    continue
                clz, method = line.split(".")
                if clz=="":
                    clz = oldclzz
                fullname = "%s.%s" % (clz, method)
                if fullname.endswith("()"):
                    print fullname, "declared with parentesis, fix it"
                    fullname = fullname[:-2]
                collection.add(fullname)
                oldclzz = clz


def collectRecursed(obj, collection):
    if hasattr(obj, "__dict__"):
        for name, item in obj.__dict__.items():
                if name.upper().startswith("TEST") or\
                    name.upper().startswith("_TEST"):
                    parse_doc(item, collection)
                    collectRecursed(item, collection)

declaredAsTested = set()

for p in glob.glob("tests/unittests/test*.py"):
    module_name= p[:-3].replace("/",".").replace(os.sep, ".")

    module = __import__(module_name).unittests
    collectRecursed(module, declaredAsTested)

missing = toTest-declaredAsTested

if missing:
    print
    print len(missing), "tests/test declarations  are missing !"
    for name in sorted(missing):
        print "    ", name

toMuch = declaredAsTested-toTest

if toMuch:
    print
    print len(toMuch), "tests/test declarations do not fit:"
    for name in toMuch:
        print "    ", name
