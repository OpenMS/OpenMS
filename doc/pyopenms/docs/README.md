pyOpenMS Documentation
======================
The source code for the pyOpenMS documentation lies here. We use sphinx to
build documentation.

Preparation

Install Sphinx (which is a Python package) and some of its modules/plugins.
We recommend doing this in a [python venv](https://docs.python.org/3/library/venv.html).

    # create it:
    python -m venv /path/to/myenv
    
    # activate it, e.g.
    
    # Linux:
    source <venv>/bin/activate
    
    # Windows:
    c:\path\to\myenv\Scripts\activate.bat

Once the environment is active, you can install all required python packages using

    cd <pyOpenMS_dir>/docs
    pip install -r requirements.txt
    python install_pyopenms.py


To build the docs run the line below (works on all operating systems), but make sure that your working dir is `<pyOpenMS_dir>/docs` (see above):

    make html
    
and check validity of links using

    make linkcheck
    



