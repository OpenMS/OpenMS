#include "tabulator.h"


#define PUSH_TO_ARRAY(array, input) \
{ \
    if (confs_no == current_size) \
    { \
        current_size *= 2; \
        array = (double *) realloc(array, current_size * sizeof(double)); \
    } \
    array[confs_no] = input; \
    confs_no++; \
}

#define MEMCPY_TO_ARRAY(array, input) \
{ \
    if (confs_no == current_size) \
    { \
        current_size *= 2; \
        array = (double *) realloc(array, current_size * sizeof(double)); \
    } \
    memcpy(array+confs_tbl_idx, generator->get_conf_signature(), allDimSizeOfInt) \
    confs_tbl_idx += generator->getAllDim(); \
}



template class Tabulator<IsoThresholdGenerator>;
template class Tabulator<IsoLayeredGenerator>;
