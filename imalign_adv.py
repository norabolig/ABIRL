import numpy as np
import argparse
import csv
import sys
from astropy.io import fits
from scipy.ndimage import affine_transform
from skimage.feature import ORB, match_descriptors
from skimage.measure import ransac
from skimage.transform import SimilarityTransform

def detect_and_match(ref, tgt, threshold=0.05, max_src=200):
    """
    Identifies stars in both images and finds the geometric 
    transformation using RANSAC for outlier rejection.
    """
    # Normalize for feature detection (0.0 to 1.0)
    def norm(img):
        img_min, img_max = np.min(img), np.max(img)
        return (img - img_min) / (img_max - img_min + 1e-6)

    # Initialize ORB detector
    detector = ORB(n_keypoints=max_src, fast_threshold=threshold)
    
    # Extract features from Reference
    detector.detect_and_extract(norm(ref))
    kp_ref, desc_ref = detector.keypoints, detector.descriptors
    
    # Extract features from Target
    detector.detect_and_extract(norm(tgt))
    kp_tgt, desc_tgt = detector.keypoints, detector.descriptors

    # Match features based on descriptors
    matches = match_descriptors(desc_ref, desc_tgt, cross_check=True)
    
    if len(matches) < 4:
        raise ValueError(f"Found only {len(matches)} matches. Lower --thresh or check image overlap.")

    # Convert coordinates (y,x) to (x,y) for the transform model
    src = kp_tgt[matches[:, 1]][:, ::-1]
    dst = kp_ref[matches[:, 0]][:, ::-1]

    # Robustly estimate the transform (Rotation + Translation + Scaling=1)
    # RANSAC ignores 'matches' that don't fit the global geometric model
    model, inliers = ransac((src, dst), SimilarityTransform, min_samples=3,
                               residual_threshold=2, max_trials=1000)
    
    # Package data for the catalog
    matched_stars = {
        'ref_xy': dst[inliers],
        'tgt_xy': src[inliers]
    }
    
    dx, dy = model.translation
    rot = np.degrees(model.rotation)
    
    return dx, dy, rot, matched_stars

def save_catalog(filename, star_data):
    """Writes the sub-pixel positions of matched stars to a CSV."""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ref_x', 'ref_y', 'target_x', 'target_y'])
        for (ref, tgt) in zip(star_data['ref_xy'], star_data['tgt_xy']):
            writer.writerow([f"{ref[0]:.3f}", f"{ref[1]:.3f}", 
                             f"{tgt[0]:.3f}", f"{tgt[1]:.3f}"])

def main():
    parser = argparse.ArgumentParser(description="Professional Star-Matching Image Alignment")
    parser.add_argument("reference", help="Reference FITS image (the 'anchor')")
    parser.add_argument("target", help="Target FITS image (to be moved)")
    parser.add_argument("output", help="Output aligned FITS filename")
    parser.add_argument("--catalog", help="Optional: CSV filename to save star matches")
    parser.add_argument("--thresh", type=float, default=0.05, help="Detection threshold (lower = more sensitive)")
    parser.add_argument("--max_src", type=int, default=200, help="Max stars to detect")
    parser.add_argument("--fill", type=float, default=0.0, help="Value for empty pixels")

    args = parser.parse_args()

    # 1. Load Data
    try:
        with fits.open(args.reference) as h_ref, fits.open(args.target) as h_tgt:
            ref_data = h_ref[0].data.astype(np.float32)
            tgt_data = h_tgt[0].data.astype(np.float32)
            header = h_tgt[0].header
    except Exception as e:
        print(f"File Error: {e}")
        sys.exit(1)

    # 2. Find Transformation
    print(f"Searching for stars in {args.target}...")
    try:
        dx, dy, rot, star_data = detect_and_match(ref_data, tgt_data, args.thresh, args.max_src)
    except ValueError as e:
        print(f"Alignment Error: {e}")
        sys.exit(1)

    print(f"Matches Verified: {len(star_data['ref_xy'])} stars")
    print(f"Calculated: ΔX={dx:.3f}, ΔY={dy:.3f}, Rotation={rot:.3f}°")

    # 3. Export Catalog
    if args.catalog:
        save_catalog(args.catalog, star_data)
        print(f"Catalog saved to: {args.catalog}")

    # 4. Apply Flux-Preserving Transformation
    # We rotate around the physical center of the image
    ny, nx = tgt_data.shape
    center = np.array([nx/2.0, ny/2.0])
    
    angle_rad = np.radians(rot)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    
    # Construct the affine matrix (inverse mapping)
    matrix = np.array([[cos_a, sin_a], [-sin_a, cos_a]])
    
    # Calculate offset to maintain rotation center and apply shifts
    # scipy uses (y, x) indexing for the transform
    shift_vec = np.array([dy, dx])
    offset = center[::-1] - np.dot(matrix, center[::-1]) - shift_vec
    
    # order=3 (Bicubic Spline) is used for flux preservation
    aligned = affine_transform(
        tgt_data, 
        matrix, 
        offset=offset, 
        order=3, 
        mode='constant', 
        cval=args.fill
    )

    # 5. Save Output
    header['HISTORY'] = f"Aligned to {args.reference} using py_align_final"
    header['ALIG_DX'] = (dx, "Horizontal pixel shift")
    header['ALIG_DY'] = (dy, "Vertical pixel shift")
    header['ALIG_ROT'] = (rot, "Rotation angle in degrees")
    header['ALIG_N'] = (len(star_data['ref_xy']), "Number of stars matched")

    fits.PrimaryHDU(aligned, header=header).writeto(args.output, overwrite=True)
    print(f"Success! Aligned FITS saved to: {args.output}")

if __name__ == "__main__":
    main()
