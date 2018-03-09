#!/usr/bin/python3

# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.


import argparse
import matplotlib.pyplot as plt

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
        help="Additional scale for undistorted grid.")
    parser.add_argument(
        "--maxorder", type=int, default=5, metavar="MAXORDER")
    parser.add_argument(
        "--saveplot", metavar="PLOTFILE")
    args = parser.parse_args()

    df = zemax.load_grid_data(args.zemaxfile)
    atrafo, posnormscale, positions, (distortions,) = \
        zemax.distortions_on_sky([df], scale=args.scale)
    print("Affine transform: focal plane (mm) -> reference grid (arcseconds):")
    print(atrafo)
    print()

    plotting.set_fontsize(medium=9)
    fig = plt.figure(figsize=(8*1.3, 6*1.3))
    plotting.make_grid_analysis(fig, positions, distortions,
                                posnormscale, args.maxorder)
    if args.saveplot:
        print()
        print("Writing plot to", args.saveplot)
        plt.savefig(args.saveplot, dpi=250)
    else:
        plt.show()
    plt.close()


if __name__ == '__main__':
    main()
