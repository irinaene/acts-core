Run the telescope tracking example
===============================

Prerequisites
-------------

Acts must be build with activated examples support
(``ACTS_BUILD_EXAMPLES=on``). ``<build>``
is used to identify the path to the build directory.

Run the telescope tracking
----------------------

.. code-block:: console

   $ <build>/bin/ActsRecTelescopeTracks \
       --input-dir=alpide_input \
       --output-dir=alpide_output \
       --output-obj=1 \
       -n 1

The ``--output-obj`` should be turned on when we want to write the trajectories to per-event obj file.

Track visualization
----------------------

With the ``--output-obj`` turned on as above, obj files named ``eventXXXXX-TelescopeTrack.obj`` will be produced in the ``output-dir``.
The obj file could be opened with e.g. meshlab by ``meshlab eventXXXXX-TelescopeTrack.obj``.

Tracking results
----------------------

The reconstructed tracks could be written out into root file (named ``telescope_tracks.root``) or csv file (not implemented yet) in the ``output-dir``.

The output root file includes a tree named with ``tracks`` with one entry representing one trajectory. From this file, one could check the information of every measurement track state on the trajectory,
i.e. the residual, pull and chi2 etc.

For example, to look at the residual of (fitted) smoothed track parameters w.r.t. measurement for local x coordinate on first plane, one could use:

.. code-block:: console

   $root <output-dir>/telescope_tracks.root 
   $tracks->Draw("res_eLOC0_smt[0]>>h(100,-0.5,0.5)")

Look at the telescope tracking performance
----------------------

The reconstruction will generate a root file (named ``performace_telescope_tracking.root``) containing tracking performance plots in the ``output-dir``. Currently, this is still quite empty since how to define the performance is still to be determined.
