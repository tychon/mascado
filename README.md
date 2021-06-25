
# MASCADO

Package for analyzing geometric distortions of optical systems.  This
software is meant for quick-look inspection of distortion data sets.


## Install

Use ``setuptools`` via

    python setup.py install [--user] [--dry-run]

in the project root directory.  This makes the ``mascado`` package
available as import and puts the executable scripts into an
appropriate directory for executables.  In the output, ``setuptools``
tells you where exactly.

Dependencies are

  - Python3.6 or higher
    - 3.5 or higher needed for matrix multiplication operator (``@``)
    - 3.6 for f-strings
  - Matplotlib
  - NumPy >= 1.14
    - 1.14 or higher needed for ``encoding`` keyword to ``numpy.genfromtxt()``
  - SciPy
  - pandas

Make sure that you are using the Python3 version of pip (explicitly
`pip3 ...` or `python3 -m pip ...`) or `conda`, because otherwise the
packages are not visible to the scripts.


## Documentation

The documentation is available at https://tychons.net/mascado/  
If you want to build the documentation yourself, get Sphinx and run

    $ make docupdate
    $ make docbuild

Then open `./docs/build/html/index.html`.


## Tests

There are a few unit tests for the package, though not nearly
extensive.  Run with

    $ python -m unittest discover --start-directory tests

or

    $ python setup.py test

or

    $ make test

## Contributing

- Codestyle: PEP8 with some exceptions: 
  - Maximum Line length 
  - 

## License

Copyright 2018 Hannes Riechert <hannes@tychons.net>  
at Max-Planck-Institute for Astronomy, Heidelberg.  
Licensed under GPL-3.0-or-later.  
You can find the full text of the license in COPYING.

If you plan to publish the obtained results, please contact the
authors:

> Hannes Riechert <hannes@tychons.net> or <riechert@mpia.de>  
> Jörg-Uwe Pott <jpott@mpia.de>
