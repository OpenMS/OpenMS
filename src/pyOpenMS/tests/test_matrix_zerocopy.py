import numpy as np
import pyopenms as poms
import pytest

def test_matrix_double_zerocopy():
    # 1. Create a NumPy array (64-bit float to match MatrixDouble)
    np_arr = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64, order='F')

    # 2. Map it to OpenMS using your new fromNdArray logic
    openms_mat = poms.MatrixDouble.fromNdArray(np_arr)

    # 3. Get the view back using your new get_matrix_as_view logic
    view_arr = openms_mat.get_matrix_as_view()

    # 4. VERIFY ZERO-COPY (Bi-directional sharing)
    openms_mat.setValue(0, 0, 99.0)
    assert view_arr[0, 0] == 99.0

    # Change in NumPy should reflect in C++
    view_arr[1, 1] = 42.0
    assert openms_mat.getValue(1, 1) == 42.0

    # 5. Verify Memory Layout (Strides)
    assert view_arr.flags['F_CONTIGUOUS'] == True

if __name__ == "__main__":
    test_matrix_double_zerocopy()
    print(" Zero-copy mapping verified")