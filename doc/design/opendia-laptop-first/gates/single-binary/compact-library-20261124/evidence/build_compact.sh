#!/bin/bash
set -ex
cd /home/singjc/.jcode/scratch/compact-20261124
FLAGS=$(ninja -C /home/singjc/Documents/github/OpenMS-r1-full-ninja-build -t commands lib/libOpenDIA_impl.a | grep OpenDIAWorkflow.cpp.o | head -1 | sed 's/^.*c++ //; s/ -o .*//; s/-MD -MT [^ ]* -MF [^ ]*//')
for tool in compact_heap compact_probe; do
c++ $FLAGS -c "$tool.cpp" -o "$tool.o"
c++ -O3 -fopenmp "$tool.o" -o "$tool" -Wl,-rpath,/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib:/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/onnxruntime-linux-x64-1.23.2/lib /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib/libOpenDIA_impl.a /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib/libOpenMS_CLI.so /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib/libOpenMS.so /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib/libOpenSwathAlgo.so /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib/libsqlite3.a -Wl,-rpath-link,/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/onnxruntime-linux-x64-1.23.2/lib:/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
done
