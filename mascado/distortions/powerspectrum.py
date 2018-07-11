
import numpy as np


def psd_spectrums(vf, posscale, grid=None):
    r"""Calculate power spectrum in both directions of vector field.

    Parameters
    ----------
    vf : :class:`mascado.distortions.polynomials.PolyVectorField`
        2D vector field.
    posscale : float
        Scale of vector field from [-1, 1] to arcsec.
        ``2*posscale`` = FOV.
    grid : ((N, N), (N, N)) 2-tuple of square float arrays
        Regular quadratic  meshgrid (xx, yy) with ij-indexing covering
        :math:`[-1,1]^2` to evaluate vector field on.  The side length
        has to be an size.  If no grid is given, one is created
        with side length of vector field degree times 8 + 1.

    Returns
    -------
    xpsd : (N, N) float array
    ypsd : (N, N) float array
        Power spectra of x- and y-component of vector field.  The zeroth
        frequency / offset is in the center.  Every pixel to the left, right,
        top, or bottom are increasing in frequency by :math:`1/\text{FOV}`.

    Notes
    -----
    The power spectrum :math:`P` in spatial frequencies :math:`k` of
    the distortion pattern :math:`D` is given by

    .. math:: P[D](k) = \left| \mathcal F[D](k) \right|^2

    Where the Fourier transform is

    .. math:: \mathcal F[D](k) = \int_{-\text{FOV}/2}^{\text{FOV}/2} D(x)\,
        e^{2 \pi i k x} \,\mathrm d x

    With the substitution :math:`x^\prime=\frac{x}{\text{FOV}/2}` we
    get :math:`D` in normalized coordinates

    .. math:: \mathcal F[D](k) = \frac{1}{2}\text{FOV} \int_{-1}^1
              \underbrace{D(x^\prime\,\text{FOV}/2)}_{D_\text{norm}(x^\prime)}\,
              e^{2 \pi i k x^\prime\, \text{FOV}/2}
              \,\mathrm d x^\prime

    The input to NumPy's discrete Fourier transform is an image and it
    corresponds to integration bounds from 0 to 1, so another
    substitution :math:`x^{\prime\prime}=\frac{x^\prime+1}{2}` is
    needed:

    .. math:: \mathcal F[D](k) &= \text{FOV} \int_0^1
        \underbrace{D(\text{FOV}\, x^{\prime\prime} - \text{FOV}/2)}
          _{D_\text{norm}^\text{img}(x^{\prime\prime})}\,
        e^{2 \pi i k (\text{FOV}\, x^{\prime\prime} - \text{FOV}/2)}
        \,\mathrm d x^{\prime\prime} \\[1em]
        &= \text{FOV} \cdot \mathcal F
             [D_\text{norm}^\text{img}](\text{FOV}\cdot k)
           \cdot \underbrace{e^{-2 \pi i k\, \text{FOV}/2}}
             _{\text{neglect phase}}

    Here we see, that the :math:`n`-th frequency in the transformed
    result is the spatial frequency in units of :math:`n/\text{FOV}`.
    """
    # make grid
    if grid is None:
        degree = max(vf.xpoly.degree, vf.ypoly.degree)
        n = degree * 8 + 1
        # don't evaluate at boundary but centered in 'pixels'
        xs = np.linspace(-1, 1, n+1)[:-1] + 1/n
        xx, yy = np.meshgrid(xs, xs, indexing='ij')
    else:
        assert xx.shape[0] == xx.shape[1] == yy.shape[0] == yy.shape[1]
        xx, yy = grid
        n = xx.shape[0]
    if n % 2 != 1:
        raise ValueError("Side length of grid has to be odd.")
    # evaluate on grid
    points = np.stack([xx.ravel(), yy.ravel()], axis=1)
    model = vf.model(points)
    # get power spectrum by fft
    mxx, myy = model[:, 0].reshape((n, n)), model[:, 1].reshape((n, n))
    sxx = np.fft.fftshift(np.absolute(np.fft.fft2(mxx)*(2*posscale))**2)
    syy = np.fft.fftshift(np.absolute(np.fft.fft2(myy)*(2*posscale))**2)
    return sxx, syy


def psd_histogram(spectrum):
    """Bin all pixels of spectrum with same frequency.

    The pixels on the border of a square centered on the spectrum's
    center are considered to have the same frequency.  Then Nyquist
    frequency corresponds to all outermost pixels.

    These pixels of the same frequency are summed up into one bin.

    Parameters
    ----------
    spectrum : (N, N)-shaped float array
        Power spectrum

    Returns
    -------
    bins : (N//2+1)-shaped float array
        Binned spectrum.
    """
    # m: number of frequencies excluding zero
    m = spectrum.shape[0] // 2
    lastsum = 0
    bins = np.zeros(m+1)
    for i in range(m+1):
        subsum = np.sum(spectrum[m-i:m+i+1, m-i:m+i+1])
        bins[i] = subsum - lastsum
        lastsum = subsum
    return bins


def psd_histogram_cumulative(spectrum):
    """Bin all pixels of spectrum up to some frequency.

    Cumulative version of ``psd_histogram``.  The pixels inside and on
    the border of a square are summed up into one bin.

    Parameters
    ----------
    spectrum : (N, N)-shaped float array
        Power spectrum

    Returns
    -------
    bins : (N//2+1)-shaped float array
        Binned spectrum.
    """
    m = spectrum.shape[0] // 2
    bins = np.zeros(m+1)
    for i in range(m+1):
        subsum = np.sum(spectrum[m-i:m+i+1, m-i:m+i+1])
        bins[i] = subsum
    return bins
