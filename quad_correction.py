import numpy
import logging

from analysis.pixel_detector import _cmc as cmc

def fix_upper_right_quadrant(img_quad, msk_quad, signal_threshold=100, min_nr_pixels_per_median=50, rescaling_factors_asics=None):
    asics = [img_quad[:,512/4*i:512/4*(i+1)] for i in range(4)]
    msk_asics = [msk_quad[:,512/4*i:512/4*(i+1)] for i in range(4)]
    
    for asic, msk_asic in zip(asics, msk_asics):
        # Vertical common mode correction column by column
        cmc(asic, msk_asic, axis=0, signal_threshold=signal_threshold, min_nr_pixels_per_median=min_nr_pixels_per_median)               
        # Horizontal common mode correction row by row
        cmc(asic, msk_asic, axis=1, signal_threshold=signal_threshold, min_nr_pixels_per_median=min_nr_pixels_per_median)

    if rescaling_factors_asics is not None:
        for f,asic in zip(rescaling_factors_asics, asics):
            if abs(f-1.0) > 1E-6:
                asic *= f

