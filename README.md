
# MASCADO

Package for analyzing geometric distortions of optical systems.  This
software is meant for quick-look inspection of distortion data sets.

Dependencies are

  - Python3.5 or higher
    - 3.5 or higher needed for matrix multiplication operator (``@``)
  - Matplotlib: e.g. `pip install matplotlib`
  - NumPy >= 1.14: e.g. `pip install numpy=1.14`
    - 1.14 or higher needed for ``encoding`` keyword to ``numpy.genfromtxt()``
  - SciPy: e.g. `pip install scipy`
  - pandas: e.g. `pip install pandas`

Make sure that you are using the Python3 version of pip (explicitly
`pip3 ...` or `python3 -m pip ...`) or `conda`, because otherwise the
packages are not visible to the scripts.

The documentation is available at https://tychons.net/mascado/  
If you want to build the documentation yourself, get Sphinx and run

    $ make docupdate
    $ make docbuild

Then open `./docs/build/html/index.html`.

There are a few unit tests for the package, though not nearly
extensive.  Run with

    $ python -m unittest discover --start-directory tests

or

    $ make test


## License

Copyright 2018 Hannes Riechert <hannes@tychons.net>  
at Max-Planck-Institute for Astronomy, Heidelberg.  
Licensed under GPL-3.0-or-later.  
You can find the full text of the license in COPYING.

If you plan to publish the obtained results, please contact the
authors for permission and latest bug fixes:

> Hannes Riechert <hannes@tychons.net> or <riechert@mpia.de>  
> Jörg-Uwe Pott <jpott@mpia.de>
