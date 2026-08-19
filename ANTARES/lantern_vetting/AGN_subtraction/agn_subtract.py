import numpy as np
from scipy.ndimage import shift
from scipy.optimize import least_squares
from scipy.ndimage import fourier_shift


def normalize_psf(psf):
    """
    Normalize a PSF image to unit total flux.
    """
    psf = np.asarray(psf, dtype=float)

    if psf.ndim != 2:
        raise ValueError("PSF must be a 2D array.")

    if psf.shape[0] % 2 == 0 or psf.shape[1] % 2 == 0:
        raise ValueError(
            "For this simple version, use a PSF stamp with odd dimensions."
        )

    total = np.sum(psf)

    if total <= 0:
        raise ValueError("PSF must have positive total flux.")

    return psf / total


def add_psf(model, psf, x, y, flux):
    """
    Add one subpixel-shifted point source to model.

    x = column
    y = row
    """

    ny, nx = model.shape
    py, px = psf.shape

    # nearest integer pixel
    ix = int(np.round(x))
    iy = int(np.round(y))

    # fractional shift
    dx = x - ix
    dy = y - iy

    # --------------------------------------------------
    # Fourier shift
    #
    # Pad first to avoid wrap-around contamination.
    # --------------------------------------------------

    pad_y = py
    pad_x = px

    psf_pad = np.pad(
        psf,
        ((pad_y, pad_y), (pad_x, pad_x)),
        mode="constant",
    )

    shifted_pad = np.fft.ifftn(
        fourier_shift(
            np.fft.fftn(psf_pad),
            shift=(dy, dx),
        )
    ).real

    # crop back to original PSF size
    shifted = shifted_pad[
        pad_y:pad_y + py,
        pad_x:pad_x + px
    ]

    # preserve normalization
    s = shifted.sum()

    if s != 0:
        shifted /= s

    # --------------------------------------------------
    # Insert into science image
    # --------------------------------------------------

    hy = py // 2
    hx = px // 2

    y0 = iy - hy
    y1 = y0 + py

    x0 = ix - hx
    x1 = x0 + px

    sy0 = max(y0, 0)
    sy1 = min(y1, ny)

    sx0 = max(x0, 0)
    sx1 = min(x1, nx)

    if sy0 >= sy1 or sx0 >= sx1:
        return

    py0 = sy0 - y0
    py1 = py0 + (sy1 - sy0)

    px0 = sx0 - x0
    px1 = px0 + (sx1 - sx0)

    model[sy0:sy1, sx0:sx1] += (
        flux * shifted[py0:py1, px0:px1]
    )


def make_agn_model(shape, psf, positions, fluxes):
    """
    Construct the summed AGN PSF model.
    """

    model = np.zeros(shape, dtype=float)

    for (x, y), flux in zip(positions, fluxes):
        add_psf(
            model=model,
            psf=psf,
            x=x,
            y=y,
            flux=flux,
        )

    return model


def _initial_fluxes(image, variance, psf, positions, valid):
    """
    Estimate initial source fluxes by weighted linear least squares
    while keeping the input positions fixed.
    """

    sigma = np.sqrt(variance[valid])

    basis = []

    for x, y in positions:

        unit_model = np.zeros_like(image, dtype=float)

        add_psf(
            model=unit_model,
            psf=psf,
            x=x,
            y=y,
            flux=1.0,
        )

        basis.append(unit_model[valid] / sigma)

    A = np.column_stack(basis)
    b = image[valid] / sigma

    fluxes, *_ = np.linalg.lstsq(A, b, rcond=None)

    # Avoid negative starting values
    fluxes = np.maximum(fluxes, 1e-10)

    return fluxes


def subtract_agn(
    image,
    variance,
    psf,
    positions,
    bad_mask=None,
    max_position_shift=2.0,
):
    """
    Simultaneously fit and subtract multiple point sources.

    Parameters
    ----------
    image : 2D ndarray
        Science image.

    variance : 2D ndarray
        Pixel variance image.

    psf : 2D ndarray
        Centered PSF stamp.

    positions : array-like, shape (N, 2)
        Approximate source positions:

            [(x1, y1),
             (x2, y2),
             ...]

        Coordinates follow the usual image convention:
        x = column, y = row.

    bad_mask : 2D bool ndarray, optional
        True for pixels that should NOT be used.

    max_position_shift : float
        Maximum allowed shift, in pixels, away from each
        initial source position.

    Returns
    -------
    result : dict
        Dictionary containing fitted positions, fluxes,
        AGN model, residual image, and scipy optimizer result.
    """

    image = np.asarray(image, dtype=float)
    variance = np.asarray(variance, dtype=float)
    positions0 = np.asarray(positions, dtype=float)

    if image.shape != variance.shape:
        raise ValueError("image and variance must have the same shape.")

    if positions0.ndim != 2 or positions0.shape[1] != 2:
        raise ValueError("positions must have shape (N_sources, 2).")

    psf = normalize_psf(psf)

    # Valid pixels
    valid = (
        np.isfinite(image)
        & np.isfinite(variance)
        & (variance > 0)
    )

    if bad_mask is not None:
        valid &= ~np.asarray(bad_mask, dtype=bool)

    # Initial flux estimates
    fluxes0 = _initial_fluxes(
        image=image,
        variance=variance,
        psf=psf,
        positions=positions0,
        valid=valid,
    )

    nsrc = len(positions0)

    # Parameter order:
    #
    # x1, y1, flux1,
    # x2, y2, flux2,
    # ...

    p0 = []

    lower = []
    upper = []

    for (x, y), flux in zip(positions0, fluxes0):

        p0.extend([x, y, flux])

        lower.extend([
            x - max_position_shift,
            y - max_position_shift,
            0.0,
        ])

        upper.extend([
            x + max_position_shift,
            y + max_position_shift,
            np.inf,
        ])

    p0 = np.asarray(p0)
    lower = np.asarray(lower)
    upper = np.asarray(upper)

    sigma = np.sqrt(variance[valid])

    def unpack(params):

        pos = np.zeros((nsrc, 2))
        flux = np.zeros(nsrc)

        for i in range(nsrc):
            pos[i, 0] = params[3 * i]
            pos[i, 1] = params[3 * i + 1]
            flux[i] = params[3 * i + 2]

        return pos, flux

    def residual_function(params):

        pos, flux = unpack(params)

        model = make_agn_model(
            shape=image.shape,
            psf=psf,
            positions=pos,
            fluxes=flux,
        )

        return (image[valid] - model[valid]) / sigma

    opt = least_squares(
        residual_function,
        p0,
        bounds=(lower, upper),
        method="trf",
        x_scale="jac",
    )

    positions_fit, fluxes_fit = unpack(opt.x)

    agn_model = make_agn_model(
        shape=image.shape,
        psf=psf,
        positions=positions_fit,
        fluxes=fluxes_fit,
    )

    residual = image - agn_model

    return {
        "positions": positions_fit,
        "fluxes": fluxes_fit,
        "model": agn_model,
        "residual": residual,
        "optimizer": opt,
    }


def recenter_psf(psf):
    """
    Recenter the PSF so that its flux centroid is exactly
    at the central pixel.

    Assumes an odd-sized PSF stamp.
    """

    psf = np.asarray(psf, dtype=float)
    psf = psf / psf.sum()

    yy, xx = np.indices(psf.shape)

    xcen = np.sum(xx * psf)
    ycen = np.sum(yy * psf)

    x_target = (psf.shape[1] - 1) / 2
    y_target = (psf.shape[0] - 1) / 2

    dx = x_target - xcen
    dy = y_target - ycen

    print(f"PSF centroid before: x={xcen:.5f}, y={ycen:.5f}")
    print(f"Shifting PSF by:      dx={dx:.5f}, dy={dy:.5f}")

    # Pad to avoid Fourier wrap-around
    py, px = psf.shape

    psf_pad = np.pad(
        psf,
        ((py, py), (px, px)),
        mode="constant",
    )

    psf_shifted_pad = np.fft.ifftn(
        fourier_shift(
            np.fft.fftn(psf_pad),
            shift=(dy, dx),
        )
    ).real

    psf_shifted = psf_shifted_pad[
        py:2 * py,
        px:2 * px,
    ]

    psf_shifted /= psf_shifted.sum()

    return psf_shifted

#################################################
# residual deconvolution
#################################################
def _psf_to_otf(psf, shape):
    """
    Convert a centered PSF stamp into an optical transfer function
    with the requested image shape.
    """

    psf = normalize_psf(psf)

    py, px = psf.shape
    ny, nx = shape

    if py > ny or px > nx:
        raise ValueError("PSF cannot be larger than the image.")

    kernel = np.zeros(shape, dtype=float)

    kernel[:py, :px] = psf

    # Move the PSF center to the FFT origin
    kernel = np.roll(kernel, -(py // 2), axis=0)
    kernel = np.roll(kernel, -(px // 2), axis=1)

    return np.fft.fft2(kernel)


def wiener_deconvolve(
    image,
    psf,
    variance=None,
    balance=None,
    pad=True,
):
    """
    Simple Wiener deconvolution.

    Parameters
    ----------
    image : 2D ndarray
        Image to deconvolve. In our application this should
        normally be the AGN-subtracted residual.

    psf : 2D ndarray
        PSF stamp.

    variance : 2D ndarray, optional
        Variance image. Used to estimate `balance` automatically
        if balance=None.

    balance : float, optional
        Wiener regularization strength.

        If None and variance is supplied, a rough value is
        estimated automatically.

    pad : bool
        Reflect-pad the image before FFT deconvolution to reduce
        edge artifacts.

    Returns
    -------
    deconvolved : 2D ndarray
    """

    image = np.asarray(image, dtype=float)

    if balance is None:

        if variance is None:
            balance = 1e-3

        else:
            variance = np.asarray(variance, dtype=float)

            good = (
                np.isfinite(image)
                & np.isfinite(variance)
                & (variance > 0)
            )

            noise_var = np.median(variance[good])
            observed_var = np.var(image[good])

            signal_var = max(
                observed_var - noise_var,
                0.1 * noise_var,
            )

            balance = noise_var / signal_var

    if pad:

        py = psf.shape[0] // 2
        px = psf.shape[1] // 2

        work = np.pad(
            image,
            ((py, py), (px, px)),
            mode="reflect",
        )

    else:
        work = image

    H = _psf_to_otf(psf, work.shape)
    G = np.fft.fft2(work)

    # Wiener filter
    F = (
        np.conj(H)
        /
        (np.abs(H) ** 2 + balance)
        * G
    )

    result = np.real(np.fft.ifft2(F))

    if pad:
        result = result[
            py:py + image.shape[0],
            px:px + image.shape[1],
        ]

    return result