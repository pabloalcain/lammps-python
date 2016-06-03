Installation
============

First build the analysis library:

.. code-block:: bash

   cd analysis/
   make

Then run, in the root directory

.. code-block:: bash

   pip install -r requirements.txt
   python setup.py install

Testing
-------

To run the tests, go to the test directory:

.. code-block:: bash

   cd test
   nosetests . --with-coverage --rednose --cover-package=analysis,postprocess,pylammps
