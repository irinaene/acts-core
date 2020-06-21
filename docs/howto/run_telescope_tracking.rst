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
       --input-dir=alpide_input_data \
       --output-dir=alpide_output \
       -n 1

Look at the telescope tracking performance
----------------------

The reconstruction will generate three root files (the name of those root files are currently not configurable via the command line) in the ``output-dir``:

*   ``telescope_tracks.root``
This includes a tree with one entry representing one trajectory. From this file, one could check the information of every measurement track state on the trajectory,
i.e. the residual, pull and chi2 etc.

*  ``performace_telescope_tracking.root``
This includes plots showing performance of the telescope tracking. Currently, it's basically empty.
