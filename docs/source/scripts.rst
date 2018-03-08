
User Scripts
============

In the ``./scripts`` directory of the repository you can find several
scripts for various tasks.  Since they are using the ``mascado``
module, you have to make that available in your python path, or run
the scripts from the repository root using::

  $ python3 -m scripts.SCRIPTNAME ARGS...


Zemax Grid Distortion Data
--------------------------

In Zemax export distortions by clicking `Analysis` in the menu.  Then
select `Miscellaneous` and `Grid Distortion`.  After a right click
into the window, a dialog with setup options pops up.  Click on `Text`
in the menu to get a text file with the data points and save that
somewhere.  The exported TXT file is Latin-1 encoded and should look
like::

  Listing of Grid Distortion Data
  
  File : ...
  Title: EELT Optical System Specification
  Date : ...
  
  
  Units are Millimeters.
  Field units are degrees
  Wavelength: 1.00000 µm
  Reference Coordinates: Xref = 0.00000E+000, Yref = 0.00000E+000
  
     i    j       X-Field       Y-Field       R-Field   Predicted X   Predicted Y        Real X        Real Y     Distortion
    -6   -6 -7.50000E-003 -7.50000E-003  1.06066E-002 -8.95381E+001 -8.95381E+001 -8.95362E+001 -8.95362E+001     -0.00
  ...
     6    6  7.50000E-003  7.50000E-003  1.06066E-002  8.95381E+001  8.95381E+001  8.95362E+001  8.95362E+001     -0.002150%
  
  Maximum distortion: -0.0021%
  Predicted coordinate ABCD matrix:
  A =    6.84e+005
  B =            0
  C =            0
  D =    6.84e+005
  
  SMIA TV Distortion -0.0011%

The ``scripts.analyze_grid`` and ``scripts.compare_grids`` scripts
read these files and use the ``X-Field``, ``Y-Field``, ``Real X``, and
``Real Y`` columns to calculate properties of the distortion pattern.
The plate scale from Zemax is not used, but an affine transformation,
defined be the least-squares solution of ``Real`` to ``Field``
coordinates translates the distorted pattern back onto the sky.  For
that reason the pattern displayed to you looks a bit different than
the pattern displayed in Zemax.  Polynomial fits are done in the
transformed coordinate set, i.e. sky coordinates.

The ``Field`` units are assumed to be degrees.  If that is not the
case, use the ``--scale`` command line argument to define a plate
scale relating to degrees.  Otherwise the displayed units are wrong.
In any case, try the ``--help`` argument for a description of script
options.

An example invocation would be::

  $ python -m scripts.analyze_grid --maxorder 6 "../Zemax Grids/ELT_MC20.TXT" --saveplot ../tmp/MC20.png
