Example file
============

The example file in ``examples/neutronstar.py`` is a typical setup for
Neutron Star Matter calculation. In this file, we set a proton
fraction of :math:`x=0.5` and a density of :math:`\rho=0.05`. Starting
from a relatively high temperature of :math:`T=2.0`, we go down to
:math:`T=0.5`, thermalizing each temperature we set. In it we can see
the use of computes in the analyzer and the logger that is
re-instantiated each time we change the system parameters.

.. literalinclude:: ../../examples/neutronstar.py
   :language: python
   :linenos:
