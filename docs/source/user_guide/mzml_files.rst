:term:`mzML` Files
==================

.. NOTE::

    This is an advanced section that dives deep into the :term:`mzML` format and we
    will investigate the file format in greater detail.  The intricacies of the
    :term:`mzML` file format are all handled by pyOpenMS internally
    and this section is only intended for the interested reader

Specifically, we will look at :term:`mzML` stores raw spectral data and how this data
is encoded in the XML format. The :term:`mzML` standard is developed by the HUPO-PSI
committee and can be read on the `official mzML website
<http://www.psidev.info/mzML>`_. It describes how to store the meta data and
the raw data for spectra and chromatograms. In short, the standard uses XML to
encode all meta data and stores the raw data using `Base64 encoding
<https://en.wikipedia.org/wiki/Base64>`_. 

Binary Encoding
---------------

:index:`To proceed <Binary encoding (mzML)>`, we will download an example file:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master"
    urlretrieve(gh + "/src/data/tiny.mzML", "test.mzML")

Let's investigate the file ``test.mzML`` and look at line 197:

.. code-block:: python
    :linenos:

    print(open("test.mzML").readlines()[197].strip())
.. code-block:: output

    <binary>AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQAAAAAAAACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkA=</binary>

We see that line 197 in the ``test.mzML`` file contains the ``binary`` XML tag
that contains a long datastring that starts with ``AAAAA`` and ends with
``AAMkA=``. This is the raw spectrum data encoded using
`Base64 <https://en.wikipedia.org/wiki/Base64>`_. We can confirm this 
by looking at some more context (lines 193 to 199 of the file):

.. code-block:: python
    :linenos:

    print("".join(open("test.mzML").readlines()[193:199]))
.. code-block:: output

    <binaryDataArray encodedLength="108" dataProcessingRef="CompassXtract_x0020_processing">
      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
      <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
      <binary>AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQAAAAAAAACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkA=</binary>
    </binaryDataArray>

We can now see that the surrounding XML tags describe how to decode the data,
namely we see that the data is describes the m/z array and is uncompressed 64
bit data. We can now open the file with pyOpenMS and print the corresponding
array which is from the second spectrum in the file:

.. code-block:: python
    :linenos:

    exp = oms.MSExperiment()
    oms.MzMLFile().load("test.mzML", exp)

    print(exp.getSpectrum(1).get_peaks()[0])

.. code-block:: output

    [ 0.  2.  4.  6.  8. 10. 12. 14. 16. 18.]

We now see that the data encoded describes 10 m/z data points that are equally
spaced in intervals of two, starting from 0 m/z and ending at 18 m/z (note:
this is a synthetic dataset).

Base64 Encoding
---------------

From the :term:`mzML` standard, we know that the array is :index:`base64 <Base64
encoding>` encoded and we can now try to decode this data ourselves. We will
first use pure Python functions :

.. code-block:: python
    :linenos:

    encoded_data = (
        b"AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQ"
        + b"AAAAAAAACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkA="
    )

    import base64, struct

    raw_data = base64.decodebytes(encoded_data)
    out = struct.unpack("<%sd" % (len(raw_data) // 8), raw_data)
    # struct.unpack('<%sf' % (len(raw_data) // 4), raw_data) # for 32 bit data
    print(out)
.. code-block:: output

    (0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0)

The code above uses the ``base64`` package on line 5 to decode the encoded data
to raw binary data. On line 6, we use the ``struct`` package to transform the
raw binary data to 64-bit floating point values. Note that ``<%sd`` is used for
64 bit data and ``<%sf`` for 32 bit data.

Alternatively, we could also use pyOpenMS to decode the same data:

.. code-block:: python
    :linenos:

    encoded_data = (
        b"AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQ"
        + b"AAAAAAAACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkA="
    )

    out = []
    oms.Base64().decode64(
        encoded_data, oms.Base64.ByteOrder.BYTEORDER_LITTLEENDIAN, out, False
    )
    print(out)
.. code-block:: output

    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]

This allows us thus to manually decode the data. We can use pyOpenMS to encode and decode 32 and 64 bit values:


.. code-block:: python
    :linenos:

    encoded_data = (
        b"AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQ"
        + b"AAAAAAAACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkA="
    )

    out = []
    oms.Base64().decode64(
        encoded_data, oms.Base64.ByteOrder.BYTEORDER_LITTLEENDIAN, out, False
    )
    print(out)

    data = oms.String()
    oms.Base64().encode64(out, oms.Base64.ByteOrder.BYTEORDER_LITTLEENDIAN, data, False)
    print(data)

    oms.Base64().encode64(out, oms.Base64.ByteOrder.BYTEORDER_LITTLEENDIAN, data, True)
    print(data)

    data = oms.String()
    oms.Base64().encode32(out, oms.Base64.ByteOrder.BYTEORDER_LITTLEENDIAN, data, False)
    print(data)

    oms.Base64().encode32(out, oms.Base64.ByteOrder.BYTEORDER_LITTLEENDIAN, data, True)
    print(data)

.. code-block:: output

    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
    b'AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQAAAAAAAACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkA='
    b'eJxjYEABDhBKAEpLQGkFKK0CpTWgtA6UNoDSRg4AZlQDYw=='
    b'AAAAAAAAAEAAAIBAAADAQAAAAEEAACBBAABAQQAAYEEAAIBBAACQQQ=='
    b'eJxjYAADBwaGBiA+AMQMjgwMCkDsAMQJQNwAxBMcAVbKBVc='

Note how encoding the data with 64 bit precision results in an output string of
length :math:`108` characters that is about twice as long compared to encoding the data
with 32 bit precision which is of length :math:`56` characters.  However, this
difference disappears when zlib compression is used and the resulting string is
shorter still.

Numpress Encoding
-----------------

We can do even better, using the :index:`numpress <numpress>` compression. The numpress algorithm
uses lossy compression, similar to jpeg compression, which is capable of
compressing data even further but at the cost of not being able to recover the
original input data exactly:

.. code-block:: python
    :linenos:

    data = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0 + 1e-8]
    print(data)
    r = []

    c = oms.NumpressConfig()
    c.np_compression = oms.MSNumpressCoder.NumpressCompression.LINEAR
    res = oms.String()
    oms.MSNumpressCoder().encodeNP(data, res, False, c)
    print(res)

    oms.MSNumpressCoder().decodeNP(res, r, False, c)
    print(r)

    c.np_compression = oms.MSNumpressCoder.NumpressCompression.PIC
    oms.MSNumpressCoder().encodeNP(data, res, False, c)
    print(res)

    oms.MSNumpressCoder().decodeNP(res, r, False, c)
    print(r)

.. code-block:: output

    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.00000001]
    b'Qc////+AAAAAAAAA/v//f4iIiIew'
    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.00000001024455]
    b'hydHZ4enx+YBYhA='
    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]

Note how the lossy numpress compression leads to even shorter data, with 16
characters for PIC compression and 28 characters for linear compression. This
makes the encoding much more efficient than lossless encoding that we have
discussed above, however this is at the price of accuracy. 

Different numpress compression schemes result in different accuracy, the LINEAR
compression scheme introduced an inaccuracy of 10e-10 while the PIC (positive
integer compression) can only store positive integers and results in greater
loss of accuracy. 
