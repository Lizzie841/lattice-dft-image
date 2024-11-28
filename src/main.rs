use image::{EncodableLayout, GenericImage, GenericImageView, Pixel};
use rustfft;

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
	let img = image::ImageReader::open("enzyme_xcrystal.png")?.decode()?;
	let mut v = Vec::<rustfft::num_complex::Complex64>::new();
	//
	println!("#parsing image");
	let img_size = (img.width() * img.height()) as usize;
	for x in 0..img.width() {
		for y in 0..img.height() {
			let norm_luma = (img.get_pixel(x, y).to_luma().0[0] as f64) / 255.0;
			v.push(rustfft::num_complex::Complex{
				re: norm_luma,
				im: 0.0,
			});
		}
	}
	//
	println!("#calc dft");
	let mut planner = rustfft::FftPlanner::new();
	let fft = planner.plan_fft_forward(4);
	fft.process(&mut v);
	//
	let mut out_v = Vec::<rustfft::num_complex::Complex<f64>>::new();
	for x in 0..img.width() {
		for y in 0..img.height() {
			let mut result = rustfft::num_complex::Complex{ re: 0.0, im: 0.0 };
			let n = 10;
			for _ in 0..n {
				let norm_luma = (img.get_pixel(x, y).to_luma().0[0] as f64) / 255.0;
				let fj = rustfft::num_complex::Complex{
					re: norm_luma,
					im: 0.0,
				};
				let xj = (x as f64) / img.width() as f64;
				let yj = (y as f64) / img.height() as f64;
				let zj = norm_luma;
				// Miller indices
				let h = (xj + xj).sqrt();
				let k = (yj + yj).sqrt();
				let l = (zj + zj).sqrt();
				result += fj * rustfft::num_complex::Complex{
					re: f64::cos(h * xj + k * yj + l * zj),
					im: -f64::sin(h * xj + k * yj + l * zj)
				};
			}
			out_v.push(result);
		}
	}
	//
	let mut max = rustfft::num_complex::Complex{ re: 1.0, im: 1.0 };
	for e in out_v.iter() {
		if e.im > max.im {
			max.im = e.im;
		} else if e.re > max.re {
			max.re = e.re;
		}
	}
	//
	println!("#write final");
	let mut img2 = image::DynamicImage::new(img.width(), img.height(), image::ColorType::Rgb8);
	let bytes = img2.as_mut_rgb8().unwrap();
	for i in 0..img_size {
		let pos = (i as u32 % img.width(), i as u32 / img.width());
		//
		let norm_luma = (img.get_pixel(pos.0, pos.1).to_luma().0[0] as f64) / 255.0;
		let norm_values = [
			out_v[i].re / max.re,
			out_v[i].im / max.im,
			norm_luma
		];
		let mut rgb: image::Rgb<u8> = image::Rgb::from([0u8; 3]);
		rgb.0[0] = ((norm_values[0] * 255.0) as u8);
		rgb.0[1] = ((norm_values[1] * 255.0) as u8);
		rgb.0[2] = ((norm_values[2] * 255.0) as u8);
		bytes[pos] = rgb;
	}
	img2.save_with_format("lattice.png", image::ImageFormat::Png)?;
	//
	Ok(())
}
