# Usage and outputs of `input_output_comp.py`

The script `input_output_comp.py` reads an input PhoSim instance catalog through
the option `--truth-table`. This catalog acts as a "truth table". Then we specify
the path to the L2 catalog that we want to match with this truth table using the
option `--input-catalog`, and the path to the calibrated exposure corresponding
the that catalog with `--calexp`. Finally we pass the path to the magnitude truth
table with `--mag-truth` (in a FITS table).

The script matches the PhoSim inputs and L2 outputs and produces some useful plots
comparing the input and output a,b (semi-major and semi-minor axis of the ellipses
that represent each galaxy) in `bias_ellipse_*.png`. It also tests the performance
of `base_ClassificationExtendedness_value` as a star/galaxy separator producing
the plot `SG_classification.png`.

This script also produces a histogram representing the redshift distribution
`redshift_hist.png`, magnitude distribution of the detected galaxies, `mag_hist_galaxies.png` 
and a scatter plot comparing the input and output magnitudes, `mag_scatter_one_chip.png`.

It also tests the number of detected objects around bright stars, producing a plot of
the pseudo-derivative of the number of galaxies detected respect to the angular distance
to the star for stars in different flux percentiles, `dn_dtheta_gal_star.png`.

Finally it produces a plot showing the input objects (blue X), detected galaxies
(cyan +), and detected stars (red +) in a 200 x 200 pixel cutout around a bright star
and a faint object, `bright_star_example.png`, and `faint_example.png`, respectively.

 
