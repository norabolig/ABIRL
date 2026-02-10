import argparse
import csv
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import maximum_filter, label, center_of_mass

def identify_sources(filename, threshold_sigma, max_sources, output_csv=None):
    # 1. Load the FITS file
    with fits.open(filename) as hdul:
        data = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim == 2:
                data = hdu.data.astype(float)
                break
        
        if data is None:
            print("Error: No 2D image data found in FITS file.")
            return

    # 2. Estimate background and noise
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    threshold_value = median + (threshold_sigma * std)
    
    # 3. Detect peaks
    neighborhood_size = 5
    data_max = maximum_filter(data, size=neighborhood_size)
    peak_mask = (data == data_max) & (data > threshold_value)
    
    labeled_array, num_features = label(peak_mask)
    
    peaks = []
    for i in range(1, num_features + 1):
        y, x = np.where(labeled_array == i)
        peaks.append((y[0], x[0], data[y[0], x[0]]))
    
    peaks = sorted(peaks, key=lambda x: x[2], reverse=True)[:max_sources]
    
    # Prepare results list
    results = []
    print(f"\n{'ID':<4} | {'X Pixel':<10} | {'Y Pixel':<10} | {'FWHM (px)':<10}")
    print("-" * 45)

    box_size = 10 
    
    for i, (py, px, val) in enumerate(peaks):
        y0, y1 = max(0, py-box_size), min(data.shape[0], py+box_size+1)
        x0, x1 = max(0, px-box_size), min(data.shape[1], px+box_size+1)
        cutout = data[y0:y1, x0:x1] - median
        cutout = np.maximum(cutout, 0)
        
        # 4. Refine Centroid
        cy_rel, cx_rel = center_of_mass(cutout)
        cx = x0 + cx_rel
        cy = y0 + cy_rel
        
        # 5. Calculate FWHM
        y_idx, x_idx = np.indices(cutout.shape)
        total_sum = np.sum(cutout)
        
        if total_sum > 0:
            var_x = np.sum(cutout * (x_idx - cx_rel)**2) / total_sum
            var_y = np.sum(cutout * (y_idx - cy_rel)**2) / total_sum
            sigma_avg = np.sqrt(np.sqrt(var_x * var_y))
            fwhm = 2.355 * sigma_avg
        else:
            fwhm = 0.0

        results.append({'id': i+1, 'x': cx, 'y': cy, 'fwhm': fwhm})
        print(f"{i+1:<4} | {cx:<10.2f} | {cy:<10.2f} | {fwhm:<10.2f}")

    # 6. Save to CSV if requested
    if output_csv:
        keys = ['id', 'x', 'y', 'fwhm']
        try:
            with open(output_csv, 'w', newline='') as f:
                dict_writer = csv.DictWriter(f, fieldnames=keys)
                dict_writer.writeheader()
                dict_writer.writerows(results)
            print(f"\nSuccess: Results saved to {output_csv}")
        except Exception as e:
            print(f"\nError saving CSV: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify point sources in a FITS image.")
    parser.add_argument("filename", help="Path to the FITS file")
    parser.add_argument("-t", "--threshold", type=float, default=5.0, 
                        help="Detection threshold in sigmas above background")
    parser.add_argument("-m", "--max", type=int, default=50, 
                        help="Maximum number of sources to report")
    parser.add_argument("-o", "--output", help="Optional path to save results as a CSV file")
    
    args = parser.parse_args()
    identify_sources(args.filename, args.threshold, args.max, args.output)
