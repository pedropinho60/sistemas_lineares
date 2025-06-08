use std::{path::Path, time::Instant};

use nalgebra::{DMatrix, io};
use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> std::io::Result<()> {
    let m1: DMatrix<f64> = io::cs_matrix_from_matrix_market(Path::new("matrizes/a3.mtx"))
        .unwrap()
        .into();

    let b1 = from_file("b3.txt")?;

    // Mede tempo de execução de cada método
    let inicio1 = Instant::now();
    run_gauss(m1.clone(), b1.clone(), 1e-6);
    let tempo1 = inicio1.elapsed();

    let inicio2 = Instant::now();
    run_lu(m1.clone(), b1.clone(), 1e-6);
    let tempo2 = inicio2.elapsed();

    let inicio3 = Instant::now();
    run_gauss_jacobi(m1.clone(), b1.clone(), 1e-6, 1000);
    let tempo3 = inicio3.elapsed();

    let inicio4 = Instant::now();
    run_gauss_seidel(m1.clone(), b1.clone(), 1e-6, 1000);
    let tempo4 = inicio4.elapsed();

    let total = inicio1.elapsed();

    println!("Tempo elim. Gauss: {}μs", tempo1.as_micros());
    println!("Tempo fat. LU: {}μs", tempo2.as_micros());
    println!("Tempo Gauss-Jacobi: {}μs", tempo3.as_micros());
    println!("Tempo Gauss-Seidel: {}μs", tempo4.as_micros());
    println!("Tempo total: {}μs", total.as_micros());

    Ok(())
}
