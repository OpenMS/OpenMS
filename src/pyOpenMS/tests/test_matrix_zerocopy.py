import numpy as np
import pyopenms as poms
import pytest

def test_matrix_double_zerocopy():
    # 1. Create a NumPy array (64-bit float to match MatrixDouble)
    np_arr = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64, order='F')

    # 2. Map it to OpenMS using your new fromNdArray logic
    openms_mat = poms.MatrixDouble.fromNdArray(np_arr)

    # 3. Get the view back using your new get_matrix_view logic
    view_arr = openms_mat.get_matrix_view()

    # 4. VERIFY ZERO-COPY & MEMORY LAYOUT (Bi-directional sharing)
    openms_mat.setValue(1, 0, 99.0)
    assert view_arr[1, 0] == 99.0
    view_arr[0, 1] = 42.0
    assert openms_mat.getValue(0, 1) == 42.0

    # 5. Verify Memory Layout (Strides)
    assert view_arr.flags['F_CONTIGUOUS'] == True

def test_matrix_view_empty():
    """get_matrix_view() should return empty array for empty matrix."""
    mat = poms.MatrixDouble()
    arr = mat.get_matrix_view()
    assert isinstance(arr, np.ndarray)
    assert arr.size == 0

if __name__ == "__main__":
    test_matrix_double_zerocopy()
    test_matrix_view_empty()
    print("Zero-copy mapping verified")