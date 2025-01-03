use ndarray::prelude::*;
use plotpy::{Curve, Plot, StrError};

fn main() {
    println!("Simple binary mask transmission function");
    let dot_size = 250.0; //Chromium dot radius in nm for a transmission mask
                          // Generate an x range from -1000 nm to 1000 nm with 128 point. Since we need to do FFT of this 2^n points are very efficient
    let x_support = Array::linspace(-1000.0, 1000.0, 128);
    let dx: f64 = 2000.0 / 128 as f64; // Pixel size
                                       // Define the binary mask transmission, 0 for outside of -250 to 250, 1 for inside and including
    let mask = x_support.map(|x_ref: &f64| {
        if (*x_ref).abs() < dot_size / 2.0 {
            1.0
        } else {
            0.0
        }
    });
    println!("{:?}", mask);
    // Plotting
    let mut curve = Curve::new();
    curve.set_label("Mask transmission curve");
    curve.draw(&x_support.to_vec(), &mask.to_vec());
    let mut plot = Plot::new();
    plot.add(&curve)
        .set_title("Mask Transmission")
        .grid_labels_legend("x positon", "Mask transmisson");
    plot.show("mask_transmision");
}
