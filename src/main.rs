use std::path::Path;

use nalgebra::{DMatrix, DVector, io};
use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> std::io::Result<()> {
    let m1: DMatrix<f64> = io::cs_matrix_from_matrix_market(Path::new("matrizes/a3.mtx"))
        .unwrap()
        .into();

    // let m1 =
    //     DMatrix::from_vec(3, 3, vec![-9.0, 5.0, 6.0, 2.0, 3.0, 1.0, -1.0, 1.0, -3.0]).transpose();

    let b1 = from_file("b3.txt")?;

    // let b1 = DVector::from_vec(vec![11.0, 4.0, -2.0]);

    run_gauss(m1.clone(), b1.clone(), 1e-6);
    run_lu(m1.clone(), b1.clone(), 1e-6);
    run_gauss_jacobi(m1.clone(), b1.clone(), 1e-6, 1000);
    run_gauss_seidel(m1.clone(), b1.clone(), 1e-6, 1000);

    Ok(())
}
