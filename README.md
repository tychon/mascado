
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
All available inspections and tests can be run via tox.

    $ pip/conda install tox  # in case its not installed 
    $ tox

There are a few unit tests for the package, though not nearly
extensive. Tests can be run with

    $ python -m unittest discover --start-directory tests

or

    $ python setup.py test

or

    $ make test

## Development workflow/Contributing

- Ensure changes conform to coding style and don't break anything by running the test cases, pylint, mypy and flake8. 
  This can be done in one go in an isolated environment. Just run `tox`.
  Add `-p auto` to parallelize for all test-environments, `-e py<yourversion>` to only test against a single version of python.
  
- Do not merge user-visible changes without documenting them in the changelog. Also create/request a release from the
  old version if it does not exist yet.

- Before merging a pull-request, make sure that at least one other person has reviewed or at least
  sanity checked the changes.
  
- Please use meaningful commit-messages, describe what and if necessary how you changed things. Not just "Fixed Stuff".
  
- When opening a pull request, please make sure all changes are self-contained and complete. In case you made some 
  incidental fixes put them in a separate branch 
  (e.g. by cherry-picking/ interactive rebasing the desired changes onto new branches. 
  
- Unless traffic in this repo increases dramatically, you don't squash changes before merging. 
  Merged branches can be deleted

### Codestyle

- Style follows PEP8, with the exception of line-length of 120 characters

- Docstrings, especially for functions, follow to the numpy guidelines.

- Prefer meaningful variable names over abbreviations. E.g. `magnitude, flux, flux_error` over `m, f, df` unless
  its very clear from context (`for i in range(...):`) or established notation.

- Write at least basic smoke-tests (i.e. "if I run this, does it explode right away?") for new functionality and when
  refactoring existing code.

## Changelog

###v1.0
Initial release

### unreleased Changes

- (example) Testdata available as package resources: `mascado.resources.example_files`

## License

Copyright 2018 Hannes Riechert <hannes@tychons.net>  
at Max-Planck-Institute for Astronomy, Heidelberg.  
Licensed under GPL-3.0-or-later.  
You can find the full text of the license in COPYING.

If you plan to publish the obtained results, please contact the
authors:

> Hannes Riechert <hannes@tychons.net> or <riechert@mpia.de>  
> Jörg-Uwe Pott <jpott@mpia.de>
