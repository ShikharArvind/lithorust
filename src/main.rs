mod fft_helpers;
use fft_helpers::{abs_of_complex_array, fftfreq, fftshift1d, ifftshift1d};
use ndarray::prelude::*;
use ndrustfft::{ndfft_par, ndifft_par, FftHandler};
use num_complex::Complex;
use plotpy::{Curve, Plot};
fn main() {
    println!("Simple binary mask transmission function");
    let dot_size = 250.0; //Chromium dot radius in nm for a transmission mask
                          // Generate an x range from -1000 nm to 1000 nm with 128 points. Since we need to do FFT of this 2^n points are very efficient
    let no_of_points: f64 = 128.0;
    let x_support = Array::linspace(-1000.0, 1000.0, no_of_points as usize);
    let dx: f64 = 2000.0 / no_of_points as f64; // Pixel size
                                                // Define the binary mask transmission, 0 for x outside of -125 to 125 and 1 for within that range. We use complex floats as it is easy to do FFT on them.
    let mask = x_support.map(|x_ref: &f64| {
        if (*x_ref).abs() < dot_size / 2.0 {
            Complex::new(1.0, 0.0)
        } else {
            Complex::new(0.0, 0.0)
        }
    });

    // Plotting mask transmission function
    let mut curve = Curve::new();
    curve.set_label("Mask transmission curve");
    curve.draw(&x_support.to_vec(), &abs_of_complex_array(&mask).to_vec());
    let mut plot = Plot::new();
    plot.add(&curve)
        .set_title("Mask Transmission")
        .grid_labels_legend("x positon", "Mask transmisson");
    plot.show("mask_transmision");

    println!("FFT of the mask transmission");
    // Mask
    let mut fft_handler = FftHandler::new(no_of_points as usize);
    let mut mask_fft = Array1::<Complex<f64>>::zeros(no_of_points as usize);
    ndfft_par(&mask, &mut mask_fft, &mut fft_handler, 0);
    //println!("{:?}", mask_fft);
    mask_fft = fftshift1d(&mask_fft);

    // Freq support
    let mut freq_support = fftfreq(no_of_points as usize, dx);
    freq_support = fftshift1d(&freq_support);

    // Plotting mask transmission function
    println!("Plotting");
    let mut curve = Curve::new();
    curve.set_label("Mask FT curve");
    curve.draw(
        &freq_support.to_vec(),
        &abs_of_complex_array(&mask_fft).to_vec(),
    );
    let mut plot = Plot::new();
    plot.add(&curve)
        .set_title("Mask FT")
        .grid_labels_legend("Frequency [1/nm]", "Mask spectrum");
    plot.show("mask_ft");

    // Abbe implementation
    fn compute_abbe_1d(
        sigma: f64,
        NA: f64,
        wavelength: f64,
        mask_ft: &Array1<Complex<f64>>,
        freq: &Array1<f64>,
    ) -> Array1<f64> {
        // Define source points (this is one approach, can also create a mask using mapv on freq and then do select on)
        let source_points: Vec<f64> = freq
            .iter()
            .filter(|&x| x.abs() <= (sigma * NA) / wavelength)
            .cloned()
            .collect();
        //println!("{:?}", source_points);
        // Initialize aerial image array
        let mut aerial_image: Array1<f64> = Array1::zeros(mask_ft.len());

        for &freq_src in &source_points {
            // Shift frequency support relative to the current source point
            let freq_mask_shift = freq.map(|x| (*x) - freq_src);

            // Shifted transfer function of the projection lens
            let pupil_shifted = freq_mask_shift.map(|x| {
                if (*x).abs() <= NA / wavelength {
                    Complex::new(1.0,0.0)
                } else {
                    Complex::new(0.0,0.0)
                }
            });

            // Shifted transfer function applied  to mask spectra
            let mask_lpf = mask_ft * pupil_shifted;


            // IFFT of mask_lpf 
            let mut ifft_handler = FftHandler::new(mask_lpf.len());
            let mut mask_lpf_ifft = Array1::<Complex<f64>>::zeros(mask_lpf.len());
            ndifft_par(&mask_lpf, &mut mask_lpf_ifft, &mut ifft_handler, 0);

            aerial_image = aerial_image + abs_of_complex_array(&mask_lpf_ifft).pow2() ;

        }

        // Normalization with number of source points
        aerial_image.map_inplace(|x| *x /= source_points.len() as f64);
        return  aerial_image;

    }
    let aerial_image_1d = compute_abbe_1d(0.7, 0.57, 365.0, &mask_fft, &freq_support);
    
    // Plotting 1D aerial image
    println!("Plotting");
    let mut curve = Curve::new();
    curve.set_label("Aerial image");
    curve.draw(
        &x_support.to_vec(),
        &aerial_image_1d.to_vec(),
    );
    let mut plot = Plot::new();
    plot.add(&curve)
        .set_title("Aerial image")
        .grid_labels_legend("x [nm]", "Intensity (au)");
    plot.show("aerial_image_1d");

}
