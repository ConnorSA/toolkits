In here:

Collection of random scripts that might be useful.

CASTEP_readers:
scripts to parse castep output files for use in the python environment.

PhononFromCastep: Gets the frequencies and q-points from a phonon castep job,
can be used to rescale q-point path from 0 to 1 scale if a predefined path is
defined. This is reading from the job.castep output file.

Note: Assumes you are using the (THz) units when running your castep job.


BandsFromCastep: Gets the energy eigenvalues for k-points along a specified
path. rescales the x axis to 0 to 1, with specified symmetry points being at
equispaced positions (i.e if 4 symmetry points are in the path, the turning
points are: 0, 0.33, 0.66, 1). This is reading from the .bands file.

reads the fermi-energy with BandsFromCastep.get_fermi_level().

Units = eV.
