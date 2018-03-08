#!/usr/bin/python3

# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.


import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

import maskastrometry.utility.affine as affine
import maskastrometry.utility.zemax as zemax
import maskastrometry.utility.plotting as plotting


def main():
    parser = argparse.ArgumentParser(
        description="Analyze properties of distorted pattern compared"
                    " to on-sky grid.")
    parser.add_argument(
        "zemaxfile1", metavar="GRIDDATA1",
        help="Path to Grid Distortion Data exported from Zemax as TXT file.")
    parser.add_argument(
        "zemaxfile2", metavar="GRIDDATA2",
        help="Path to Grid Distortion Data exported from Zemax as TXT file.")
    parser.add_argument(
        "--scale", type=float, default=1, metavar="SCALE",
        help="Plate scale for undistorted grid.")
    parser.add_argument(
        "--maxorder", type=int, default=5, metavar="MAXORDER")
    parser.add_argument(
        "--saveplot", metavar="PLOTFILE",
        help="Save plot to file and don't show on screen.")
    args = parser.parse_args()

    df1 = zemax.load_grid_data(args.zemaxfile1)
    df2 = zemax.load_grid_data(args.zemaxfile2)
    assert len(df1)**0.5 % 1 == 0, "Not a square number of rows in file 1?"
    assert len(df2)**0.5 % 1 == 0, "Not a square number of rows in file 2?"
    if not zemax.have_same_grid([df1, df2]):
        print("Both catalogs are required to have the same reference grid")
        print("(x-field, y-field columns).")
        sys.exit(1)

    # reference grid
    refpos = df1.loc[:, ['x-field', 'y-field']].as_matrix()  # degree
    refpos = refpos * args.scale
    refpos = refpos * 3600  # degree to arcsec
    # distorted grid in focal plane
    dpos1 = df1.loc[:, ['x-real', 'y-real']].as_matrix()  # mm
    dpos2 = df2.loc[:, ['x-real', 'y-real']].as_matrix()  # mm

    # remove same affine trafo in both catalogs
    atrafo = affine.affine_lstsq(dpos1, refpos)
    dpos1 = affine.affine_trafo(dpos1, atrafo)
    dpos2 = affine.affine_trafo(dpos2, atrafo)
    # distortions
    dis1 = dpos1 - refpos
    dis2 = dpos2 - refpos
    # drift
    drift = dis2 - dis1
    # now everything is in arcseconds
    print("Affine transform: distorted (mm) -> reference grid (arcseconds):")
    print(atrafo)
    print("Is applied to both catalogs GD1 and GD2.")
    print("Drift = GD2 - GD1")

    # normalize positions
    posmin, posmax = np.min(refpos), np.max(refpos)
    posscale = (posmax - posmin) / 2
    posshift = (posmax + posmin) / 2
    positions = refpos / posscale - posshift

    print()
    plotting.set_fontsize(medium=9)
    fig = plt.figure(figsize=(8*1.3, 6*1.3))
    plotting.make_grid_analysis(fig, positions, drift,
                                posscale, args.maxorder, name="drift")
    if args.saveplot:
        print()
        print("Writing plot to", args.saveplot)
        plt.savefig(args.saveplot, dpi=250)
    else:
        plt.show()


if __name__ == '__main__':
    main()
