Run the telescope tracking example
===============================

Prerequisites
-------------

Acts must be build with activated examples and Pythia8 support
(``ACTS_BUILD_EXAMPLES=on``) to enable the fast simulation. ``<build>``
is used to identify the path to the build directory.

Run the telescope tracking
----------------------

.. code-block:: console

   $ <build>/bin/ActsRecTelescopeTracks \
       --input-dir=alpide_input_data \
       --output-dir=alpide_output

The magnetic field setup should be consistent between simulation and truth tracking. 

Look at the telescope tracking performance
----------------------

The reconstruction will generate three root files (the name of those root files are currently not configurable via the command line) in the ``output-dir``:

*   ``telescope_tracks.root``
This includes a tree with one entry representing one trajectory. From this file, one could check the information of every measurement track state on the trajectory.

*  ``performace_telescope_tracking.root``
This includes a tree showing performance of the truth track finding.
