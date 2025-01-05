use ndarray::{concatenate, s, Array1, Axis};
use num_complex::Complex;

// One dimensional fftshit (swap left and right havles of the array)
pub fn fftshift1d<T: Clone>(x: &Array1<T>) -> Array1<T> {
    let len = (*x).len();
    let mut middle_value = len / 2; // Default middle value, works if the array len is even
    if len % 2 != 0 {
        middle_value = (len + 1) / 2; // Ceil of middle values, for odd lengthed arrays
    }
    let left_side = (*x).slice(s![0..middle_value]);
    let right_side = (*x).slice(s![middle_value..]);
    return  concatenate![Axis(0), right_side, left_side]; // Swap the right side to left and vice versa
}

// Generate  Discrete Fourier Transform sample frequencies 1D 
pub fn fftfreq(n: usize , d: f64) -> Array1<f64> {
    let spacing = 1.0 / ((n as f64 * d) );
    let mut spatial_freq = Array1::zeros(n);
    let N = ((n - 1) / 2) + 1;
    let p1 = Array1::range(0.0, N as f64, 1.0);
    spatial_freq.slice_mut(s![..N as usize]).assign(&p1); 
    let p2 = Array1::range(-1.0*((n/2) as f64), 0.0, 1.0);
    spatial_freq.slice_mut(s![N as usize..]).assign(&p2); 
    spatial_freq.mapv_inplace(|x| x*spacing);
    spatial_freq

}
pub fn abs_of_complex_array (x: &Array1<Complex<f64>>) -> Array1<f64> {
    (*x).mapv(|z| z.norm())
}