.. _usage-plugins-xrayDiffraction:

X-ray diffraction
-----------------

This plugin calculates X-ray diffraction from particle positions.

.. math::

   I({\bf q}) &= \frac{F({\bf q})F({\bf q})^*}{Np} \\
   F({\bf q}) &= \sum^{N}_{j=1} f_j w_j \exp(i {\bf q \cdot \bf r_j})

============================== ================================================================================
Variable                       Meaning
============================== ================================================================================
:math:`\bf r_j`                The position of macro particle *j* .
:math:`\bf q`                  Scattering vector, which is defined as :math:`|{\bf q}| = \frac{4 \pi \sin \theta}{\lambda}`.
:math:`f_j`                    The form factor of macro particle *j*.
:math:`w_j`                    The weighting of macro particle *j*.
:math:`Np`                     Number of real particles (total weighting of macro particles).
:math:`N`                      Number of macro particles.
============================== ================================================================================

This is based on the kinematic model of scattering.
It can calculate scattering intensity from all the particle positions for each time step.
Right now, only the plane wave condition is considered.
Multiple Scattering CAN NOT be handled in this model.

Temporal and transversal X-ray envelope will be taken into account in the future.

For free electrons, :math:`f_j = 1`.
For ions, further computations are required, which is not implemented yet.

.cfg file
^^^^^^^^^

For a specific charged species ``<species>``, e.g. ``e``, the scattering can be computed by the following commands.

========================================= ==============================================================================================================================
Command line option                       Description
========================================= ==============================================================================================================================
``--<species>_xrayDiffraction.period``    Gives the number of time steps between which the scattering intensity should be calculated.
                                          Default is ``0``, which means that the scattering intensity is never calculated.
                                          Using ``1`` calculates the scattering intensity constantly.
``--<species>_xrayDiffraction.qx_max``    Upper bound of reciprocal space range in qx direction. The unit is :math:`Å^{-1}`
                                          Default is ``5``.
``--<species>_xrayDiffraction.qy_max``    Upper bound of reciprocal space range in qy direction. The unit is :math:`Å^{-1}`
                                          Default is ``5``.
``--<species>_xrayDiffraction.qz_max``    Upper bound of reciprocal space range in qz direction. The unit is :math:`Å^{-1}`
                                          Default is ``0``.
``--<species>_xrayDiffraction.qx_min``    Lower bound of reciprocal space range in qx direction. The unit is :math:`Å^{-1}`
                                          Default is ``5``.
``--<species>_xrayDiffraction.qy_min``    Lower bound of reciprocal space range in qy direction. The unit is :math:`Å^{-1}`
                                          Default is ``5``.
``--<species>_xrayDiffraction.qz_min``    Lower bound of reciprocal space range in qz direction. The unit is :math:`Å^{-1}`
                                          Default is ``0``.
``--<species>_xrayDiffraction.n_qx``      Number of scattering vectors in qx direction.
                                          Default is ``100``.
``--<species>_xrayDiffraction.n_qy``      Number of scattering vectors in qy direction.
                                          Default is ``100``.
``--<species>_xrayDiffraction.n_qz``      Number of scattering vectors in qz direction.
                                          Default is ``1``.
========================================= ==============================================================================================================================


Output
^^^^^^

``<species>_xrayDiffraction.[timestep].dat``

An *ASCII* file that contains scattering intensity for each scattering vector defined by the reciprocal space range. The first 3 columns are the 3 components of scattering vector, and the 4th column is the corresponding scattering intensity.

``<species>_xrayDiffraction.[timestep].log``

An *ASCII* file that contains the number of real particles and the number of macro particles.

Memory Complexity
^^^^^^^^^^^^^^^^^

Accelerator
"""""""""""

One complex number per each scattering vector.

Host
""""

Two complex numbers per each scattering vector.

Known Limitations
^^^^^^^^^^^^^^^^^

- the plugin correctly handles only electrons, not ions
- may be time-consuming for large amount of macro particles

References
^^^^^^^^^^

.. [E2018]
    J.C. E, L. Wang, S. Chen, Y.Y. Zhang, and S.N. Luo.
	*GAPD: a GPU-accelerated atom-based polychromatic diffraction simulation code*,
	Journal of Synchrotron Radiation 25, pp. 604-611 (2018),
	`DOI:10.1107/S1600577517016733 <https://doi.org/10.1107/S1600577517016733>`_
