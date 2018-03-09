#!/usr/bin/python3

# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.


import sys
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
    atrafo, posnormscale, positions, (dis1, dis2) = \
        zemax.distortions_on_sky([df1, df2], scale=args.scale)
    print("Affine transform: focal plane (mm) -> reference grid (arcseconds):")
    print(atrafo)
    print("Is applied to both catalogs GD1 and GD2.")
    print("Drift = GD2 - GD1")
    print()

    drift = dis2 - dis1

    plotting.set_fontsize(medium=9)
    fig = plt.figure(figsize=(8*1.3, 6*1.3))
    plotting.make_grid_analysis(fig, positions, drift,
                                posnormscale, args.maxorder, name="drift")
    if args.saveplot:
        print()
        print("Writing plot to", args.saveplot)
        plt.savefig(args.saveplot, dpi=250)
    else:
        plt.show()
    plt.close()


if __name__ == '__main__':
    main()
