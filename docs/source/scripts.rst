
User Scripts
============

In the ``./scripts`` directory of the repository you can find several
scripts for various tasks.  Since they are using the ``mascado``
module, you have to make that available in your python path, or run
the scripts from the repository root using::

  $ python3 -m scripts.SCRIPTNAME ARGS...


Zemax Grid Distortion Data
--------------------------

Usage
^^^^^

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

If your exported Zemax file has a different encoding than Latin-1 (ISO
8859-1), use the ``--encoding`` argument of the scripts.  Although
UTF-8 and Latin-1 are indistinguishable for ASCII characters.


Explanation
^^^^^^^^^^^

The output figure has four panels.  The upper left panel shows the
distortions in sky coordinates after the affine transform.  The
reference point in Zemax is irrelevant, because a least-squares
solution is used to relate the detector to on-sky coordinates.  The
upper right panel shows the residuals of the highest order fit and
their marginalized distributions.  The arrow key is, similar to the
first panel, twice the RMS of the residuals.  With simulation data the
residuals often show small scale systematic patterns that may come
from the accuracy of Zemax's simulations.

The steps performed to obtain the polynomial fit are as follows:

#. Read on-sky grid :math:`x_\text{ref}` and detector grid
   :math:`x_\text{detector}` from data file.
#. Compute affine transform :math:`A` from detector coordinates to sky
   coordinates and apply: :math:`\hat x = A x_\text{detector}`.  (The
   ``compare_grids`` script applies the same trafo to both input
   grids.)
#. The distortions are defined as :math:`d = \hat x - x_\text{ref}`.
   The affine trafo was chosen in the previous step, such that the
   distortions are minimal in the least-squares sense.  (The drift in
   ``compare_grids`` is calculated by subtracting the two distortions
   grids: :math:`d = d_2 - d_1`.)
#. Normalize the reference positions :math:`x_\text{norm} =
   x_\text{ref} / \text{scale} - \text{shift}` by a shift and a scale
   into the domain :math:`[-1, 1]\times[-1, 1]`.  This normalization
   is transparent in the script (you don't have to care about it) but
   an important thing to remember when reusing the code.  The offset
   of the normalization is not restored, such that the resulting field
   is always centered on :math:`(0, 0)`.  The magnitude of the
   distortions is not scaled during calculations.
#. Fit a polynomial vector field :math:`P(x_\text{norm}) \approx d`.

The lower two panels are scatter plots showing properties of the
polynomial fit at different orders.  *Since the polynomial is used to
fit the distortions only, the zeroth and first order terms are
incomplete, because the affine transform was already removed.* The
drift in ``compare_grids`` compares two distortion solutions after the
same affine trafo, so the zeroth and first order terms have a
significant interpretation that is only influenced by the scale of the
affine trafo.

The lower left panel shows the RMS of residuals :math:`(r_{x,i},
r_{y,i})`, :math:`i=1...n`, calculated as

.. math::
   \operatorname{rms} r = \sqrt{\frac{\sum_i^n \left(r_{x,i}^2 + r_{y,i}^2\right)}{2 n}}
     = \operatorname{rms}\{\operatorname{rms} r_x, \operatorname{rms} r_y\}

The lower right panel shows the RMS of the distortions encoded by the
terms of a specific order in the highest order fit.  It is calculated
by setting the coefficients for all other terms to zero and evaluating
the model on a regular :math:`100\times100` point grid.


Input Format Variants
^^^^^^^^^^^^^^^^^^^^^

Scripts and macros in Zemax might produce different outputs, so some
slightly different formats might be implemented over time.  Currently
there is one variant additional to the default format.

The format variant can be chosen using the ``--format LETTER``
argument.

Variant B::

  Executing PATH
  start
  A -7.44068E+005
  B 1.70508E-003
  C -3.78466E-005
  D -7.50941E+005
  EFFL  -6.31490E+005   mm
         Npoint         Input_X deg     Input_Y         Distorted_X mm          Distorted_Y
  1.00000E+000 -7.49999E-003  -7.49995E-003  9.44980E+001  9.66661E+001
  2.00000E+000 -7.49993E-003  -6.74997E-003  9.47787E+001  8.72730E+001
  ...

In most cases the exact content of non-tabular lines is ignored,
because only a specific number of lines is skipped in the beginning
and end.
