import numpy as np
import argparse
import csv
import sys
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import affine_transform, map_coordinates
from skimage.feature import ORB, match_descriptors
from skimage.measure import ransac
from skimage.transform import SimilarityTransform

def detect_and_match(ref, tgt, threshold=0.05, max_src=200):
    """Identifies stars in both images and finds the geometric transformation."""
    def norm(img):
        img_min, img_max = np.min(img), np.max(img)
        return (img - img_min) / (img_max - img_min + 1e-6)

    detector = ORB(n_keypoints=max_src, fast_threshold=threshold)
    
    detector.detect_and_extract(norm(ref))
    kp_ref, desc_ref = detector.keypoints, detector.descriptors
    
    detector.detect_and_extract(norm(tgt))
    kp_tgt, desc_tgt = detector.keypoints, detector.descriptors

    matches = match_descriptors(desc_ref, desc_tgt, cross_check=True)
    
    if len(matches) < 4:
        raise ValueError(f"Found only {len(matches)} matches. Lower --thresh or check image overlap.")

    src = kp_tgt[matches[:, 1]][:, ::-1]
    dst = kp_ref[matches[:, 0]][:, ::-1]

    model, inliers = ransac((src, dst), SimilarityTransform, min_samples=3,
                               residual_threshold=2, max_trials=1000)
    
    matched_stars = {
        'ref_xy': dst[inliers],
        'tgt_xy': src[inliers]
    }
    
    dx, dy = model.translation
    rot = np.degrees(model.rotation)
    
    return dx, dy, rot, matched_stars

def align_by_wcs(ref_header, tgt_header, tgt_data, fill_val=0.0):
    """
    Aligns the target image to the reference grid using FITS WCS keywords 
    (CRVAL, CRPIX, CD matrix, etc.).
    """
    # Initialize WCS objects from headers
    # Astropy natively reads CRVAL, CRPIX, and CD1_1 etc., from headers here
    w_ref = WCS(ref_header)
    w_tgt = WCS(tgt_header)

    # 1. Create a grid of pixel coordinates for the output (reference image shape)
    ny, nx = ref_header['NAXIS2'], ref_header['NAXIS1']
    y, x = np.mgrid[0:ny, 0:nx]

    # 2. Convert reference pixel coordinates to world coordinates (RA, Dec)
    # origin=0 because numpy/scipy use 0-based indexing
    ra, dec = w_ref.wcs_pix2world(x, y, 0)

    # 3. Convert those world coordinates back to target image pixel coordinates
    x_tgt, y_tgt = w_tgt.wcs_world2pix(ra, dec, 0)

    # 4. Map coordinates expects input in (y, x) order for interpolation
    coords = np.array([y_tgt, x_tgt])

    print("Reprojecting target frame pixels using WCS coordinates...")
    # order=3 is bicubic spline (flux-preserving for rigid transforms)
    aligned = map_coordinates(tgt_data, coords, order=3, mode='constant', cval=fill_val)
    
    return aligned

def save_catalogue(filename, star_data):
    """Writes the sub-pixel positions of matched stars to a CSV."""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ref_x', 'ref_y', 'target_x', 'target_y'])
        for (ref, tgt) in zip(star_data['ref_xy'], star_data['tgt_xy']):
            writer.writerow([f"{ref[0]:.3f}", f"{ref[1]:.3f}", 
                             f"{tgt[0]:.3f}", f"{tgt[1]:.3f}"])

def main():
    parser = argparse.ArgumentParser(description="Professional Image Alignment (WCS or Star-Matching)")
    parser.add_argument("reference", help="Reference FITS image")
    parser.add_argument("target", help="Target FITS image")
    parser.add_argument("output", help="Output aligned FITS filename")
    parser.add_argument("--wcs", action="store_true", help="Use WCS header keywords to align instead of stars")
    parser.add_argument("--catalogue", help="Optional: CSV filename to save star matches (only works if not using --wcs)")
    parser.add_argument("--thresh", type=float, default=0.05, help="Star detection sensitivity")
    parser.add_argument("--max_src", type=int, default=200, help="Max stars to detect")
    parser.add_argument("--fill", type=float, default=0.0, help="Value for empty pixels")

    args = parser.parse_args()

    # Load Data
    try:
        with fits.open(args.reference) as h_ref, fits.open(args.target) as h_tgt:
            ref_data = h_ref[0].data.astype(np.float32)
            tgt_data = h_tgt[0].data.astype(np.float32)
            
            ref_hdr = h_ref[0].header
            tgt_hdr = h_tgt[0].header
    except Exception as e:
        print(f"File Error: {e}")
        sys.exit(1)

    # Process alignment based on selected mode
    if args.wcs:
        print(f"Aligning using WCS data present in headers...")
        try:
            aligned = align_by_wcs(ref_hdr, tgt_hdr, tgt_data, fill_val=args.fill)
            
            # Update the output header to adopt the reference's WCS
            # Since the image is now locked to the reference grid, it shares its WCS!
            out_hdr = tgt_hdr.copy()
            for key in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
                if key in ref_hdr:
                    out_hdr[key] = ref_hdr[key]
                    
            out_hdr['HISTORY'] = f"WCS-Aligned to {args.reference} using py_align_final"
            
        except Exception as e:
            print(f"WCS Alignment Error: {e}")
            print("Ensure both headers have valid astrometric keywords.")
            sys.exit(1)
            
    else:
        # Fall back to Star Matching Mode
        print(f"Searching for stars to match images...")
        try:
            dx, dy, rot, star_data = detect_and_match(ref_data, tgt_data, args.thresh, args.max_src)
        except ValueError as e:
            print(f"Alignment Error: {e}")
            sys.exit(1)

        print(f"Matches Verified: {len(star_data['ref_xy'])} stars")
        print(f"Calculated: ΔX={dx:.3f}, ΔY={dy:.3f}, Rotation={rot:.3f}°")

        if args.catalogue:
            save_catalogue(args.catalogue, star_data)
            print(f"Catalog saved to: {args.catalogue}")

        # Geometric Warp
        ny, nx = tgt_data.shape
        center = np.array([nx/2.0, ny/2.0])
        angle_rad = np.radians(rot)
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
        matrix = np.array([[cos_a, sin_a], [-sin_a, cos_a]])
        shift_vec = np.array([dy, dx])
        offset = center[::-1] - np.dot(matrix, center[::-1]) - shift_vec
        
        aligned = affine_transform(tgt_data, matrix, offset=offset, order=3, mode='constant', cval=args.fill)

        out_hdr = tgt_hdr.copy()
        out_hdr['HISTORY'] = f"Feature-Aligned to {args.reference} using py_align_final"
        out_hdr['ALIG_DX'] = (dx, "Horizontal pixel shift")
        out_hdr['ALIG_DY'] = (dy, "Vertical pixel shift")
        out_hdr['ALIG_ROT'] = (rot, "Rotation angle in degrees")
        out_hdr['ALIG_N'] = (len(star_data['ref_xy']), "Number of stars matched")

    # Save Output
    fits.PrimaryHDU(aligned, header=out_hdr).writeto(args.output, overwrite=True)
    print(f"Success! Aligned FITS saved to: {args.output}")

if __name__ == "__main__":
    main()
