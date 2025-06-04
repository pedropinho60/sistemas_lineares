use nalgebra::DMatrix;
use std::io::{self, Write};
use std::{fs::File, path::Path};

fn main() -> io::Result<()> {
    let mut file = File::create("matrizes/b3.txt")?;

    let m1: DMatrix<f64> = nalgebra::io::cs_matrix_from_matrix_market(Path::new("matrizes/a3.mtx"))
        .unwrap()
        .into();

    let b: Vec<f64> = (0..m1.nrows())
        .map(|_| rand::random_range(1.0..=50.0))
        .collect();

    for value in b {
        writeln!(file, "{}", value)?;
    }

    Ok(())
}
