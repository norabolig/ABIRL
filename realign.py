import numpy as np
import argparse
from astropy.io import fits
from scipy.ndimage import affine_transform
from scipy.ndimage import rotate as scipy_rotate

def align_image(data, dx, dy, angle_deg, center=None, mode='constant', cval=0.0):
    """
    Performs a flux-preserving transformation (shift and rotation).
    """
    # Convert angle to radians for matrix math
    angle_rad = np.radians(angle_deg)
    
    # Define the center of rotation (default to image center)
    ny, nx = data.shape
    if center is None:
        center = np.array([nx / 2.0, ny / 2.0])

    # The transformation matrix (Rotation + Shift)
    # We define the inverse mapping: from output pixel to input pixel
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)
    
    # Matrix for rotation
    # [[cos, sin], [-sin, cos]]
    matrix = np.array([[cos_a, sin_a], [-sin_a, cos_a]])
    
    # Calculate the offset to keep rotation centered and add the user shifts
    # offset = center - matrix dot center - shifts
    # Note: scipy uses (y, x) indexing, so we flip dx, dy
    shift_vec = np.array([dy, dx])
    offset = center[::-1] - np.dot(matrix, center[::-1]) - shift_vec

    # Apply transformation
    # order=3 is cubic spline (good balance of speed and flux conservation)
    aligned_data = affine_transform(
        data, 
        matrix, 
        offset=offset, 
        order=3, 
        mode=mode, 
        cval=cval
    )
    
    return aligned_data

def main():
    parser = argparse.ArgumentParser(description="Replicate IRAF imalign functionality.")
    parser.add_argument("input", help="Input FITS file")
    parser.add_argument("output", help="Output FITS file")
    parser.add_argument("--dx", type=float, default=0.0, help="Horizontal shift (pixels)")
    parser.add_argument("--dy", type=float, default=0.0, help="Vertical shift (pixels)")
    parser.add_argument("--rot", type=float, default=0.0, help="Rotation angle (degrees, counter-clockwise)")
    parser.add_argument("--fill", type=float, default=0.0, help="Value for pixels outside the frame")
    parser.add_argument("--crop", action="store_true", help="Crop the image to remove padding (not implemented in this simple version)")

    args = parser.parse_args()

    # Load data
    with fits.open(args.input) as hdul:
        header = hdul[0].header
        data = hdul[0].data.astype(np.float32)

    # Process alignment
    aligned = align_image(data, args.dx, args.dy, args.rot, cval=args.fill)

    # Update header to reflect changes
    header['HISTORY'] = f"Aligned: dx={args.dx}, dy={args.dy}, rot={args.rot}"
    
    # Save result
    hdu = fits.PrimaryHDU(aligned, header=header)
    hdu.writeto(args.output, overwrite=True)
    print(f"Successfully saved aligned image to {args.output}")

if __name__ == "__main__":
    main()
