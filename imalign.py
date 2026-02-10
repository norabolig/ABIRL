import numpy as np
import argparse
from astropy.io import fits
from scipy.ndimage import affine_transform
from skimage.registration import phase_cross_correlation
from skimage.transform import warp_polar

def discover_transform(reference, target):
    """
    Finds the rotation and translation (dx, dy) to align target to reference.
    Uses Fourier-Mellin approach.
    """
    # 1. Estimate Rotation using Log-Polar Transform
    # We transform both images to log-polar coordinates
    radius = np.hypot(*reference.shape) / 2
    ref_polar = warp_polar(reference, radius=radius, scaling='log')
    tgt_polar = warp_polar(target, radius=radius, scaling='log')

    # Phase correlation in polar space gives us the rotation angle
    shifts_polar, error, phasediff = phase_cross_correlation(ref_polar, tgt_polar, upsample_factor=10)
    
    # The 'y' shift in polar space corresponds to the rotation angle
    # (assuming 360 degrees mapped across the height of the polar image)
    detected_rotation = shifts_polar[0] * (360.0 / ref_polar.shape[0])

    # 2. Correct rotation temporarily to find translation
    # (Using a simple rotation for the discovery phase)
    # Note: We use -detected_rotation to 'undo' the difference
    from scipy.ndimage import rotate
    temp_target = rotate(target, -detected_rotation, reshape=False, order=1)

    # 3. Estimate Translation (dx, dy)
    shifts, error, phasediff = phase_cross_correlation(reference, temp_target, upsample_factor=10)
    dy, dx = shifts # Scipy/Numpy use (y, x)

    return dx, dy, -detected_rotation

def apply_alignment(data, dx, dy, angle_deg, fill_val=0.0):
    """
    The flux-preserving warp from the previous implementation.
    """
    angle_rad = np.radians(angle_deg)
    ny, nx = data.shape
    center = np.array([nx / 2.0, ny / 2.0])
    
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    matrix = np.array([[cos_a, sin_a], [-sin_a, cos_a]])
    
    shift_vec = np.array([dy, dx])
    offset = center[::-1] - np.dot(matrix, center[::-1]) - shift_vec

    return affine_transform(data, matrix, offset=offset, order=3, cval=fill_val)

def main():
    parser = argparse.ArgumentParser(description="Auto-align a target FITS to a reference FITS.")
    parser.add_argument("reference", help="Reference FITS file")
    parser.add_argument("target", help="Target FITS file to be moved")
    parser.add_argument("output", help="Output FITS file")
    parser.add_argument("--fill", type=float, default=0.0, help="Background fill value")
    
    args = parser.parse_args()

    with fits.open(args.reference) as h_ref, fits.open(args.target) as h_tgt:
        ref_data = h_ref[0].data.astype(np.float32)
        tgt_data = h_tgt[0].data.astype(np.float32)
        header = h_tgt[0].header

    print("Analyzing image offsets and rotation...")
    
    # Discovery
    dx, dy, rot = discover_transform(ref_data, tgt_data)
    
    print(f"Detected: dx={dx:.2f}, dy={dy:.2f}, rotation={rot:.2f}°")

    # Alignment
    aligned = apply_alignment(tgt_data, dx, dy, rot, fill_val=args.fill)

    # Header update
    header['HISTORY'] = f"Auto-Aligned to {args.reference}"
    header['ALIG_DX'] = (dx, "Horizontal shift detected")
    header['ALIG_DY'] = (dy, "Vertical shift detected")
    header['ALIG_ROT'] = (rot, "Rotation detected")

    fits.PrimaryHDU(aligned, header=header).writeto(args.output, overwrite=True)
    print(f"Alignment complete. Saved to {args.output}")

if __name__ == "__main__":
    main()
