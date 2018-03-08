#!/usr/bin/python3

# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.


import argparse
import numpy as np
import matplotlib.pyplot as plt

import mascado.utility.affine as affine
import mascado.utility.zemax as zemax
import mascado.utility.plotting as plotting


def main():
    parser = argparse.ArgumentParser(
        description="Analyze properties of distorted pattern compared"
                    " to on-sky grid.")
    parser.add_argument(
        "zemaxfile", metavar="ZEMAXGRIDDATA",
        help="Path to Grid Distortion Data exported from Zemax as TXT file.")
    parser.add_argument(
        "--scale", type=float, default=1, metavar="SCALE",
        help="Plate scale for undistorted grid.")
    parser.add_argument(
        "--maxorder", type=int, default=5, metavar="MAXORDER")
    parser.add_argument(
        "--saveplot", metavar="PLOTFILE")
    args = parser.parse_args()

    df = zemax.load_grid_data(args.zemaxfile)
    assert len(df)**0.5 % 1 == 0, "Not a square number of rows?"

    # reference grid
    refpos = df.loc[:, ['x-field', 'y-field']].as_matrix()  # degree
    refpos = refpos * args.scale
    refpos = refpos * 3600  # degree to arcsec
    # distorted grid in focal plane
    dpos = df.loc[:, ['x-real', 'y-real']].as_matrix()  # mm

    # remove affine part
    atrafo = affine.affine_lstsq(dpos, refpos)
    dpos = affine.affine_trafo(dpos, atrafo)
    # now everything is in arcseconds
    distortions = dpos - refpos
    print("Affine transform: distorted (mm) -> reference grid (arcseconds):")
    print(atrafo)

    # normalize positions
    posmin, posmax = np.min(refpos), np.max(refpos)
    posscale = (posmax - posmin) / 2
    posshift = (posmax + posmin) / 2
    positions = refpos / posscale - posshift

    print()
    plotting.set_fontsize(medium=9)
    fig = plt.figure(figsize=(8*1.3, 6*1.3))
    plotting.make_grid_analysis(fig, positions, distortions,
                                posscale, args.maxorder)
    if args.saveplot:
        print()
        print("Writing plot to", args.saveplot)
        plt.savefig(args.saveplot, dpi=250)
    else:
        plt.show()


if __name__ == '__main__':
    main()
