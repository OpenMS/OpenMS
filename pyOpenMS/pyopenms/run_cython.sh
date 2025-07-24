AUTOWRAP=$(python -c "import autowrap; print autowrap.__path__[0]")
cython --cplus pyopenms.pyx -I ../pxds  -I $AUTOWRAP/data_files -I $AUTOWRAP/data_files/boost
