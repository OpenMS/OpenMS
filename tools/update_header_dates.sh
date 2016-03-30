find . -name "*.h" -o -name "*.doxygen" -o -name "*.txt" -o -name "*.cpp" | xargs sed -i 's/and Freie Universitaet Berlin 2002-20../and Freie Universitaet Berlin 2002-2016/g'
