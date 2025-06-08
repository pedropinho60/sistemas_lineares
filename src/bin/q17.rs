use nalgebra::{DMatrix, DVector};
use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> std::io::Result<()> {
    let m1 = DMatrix::from_vec(
        5,
        5,
        vec![
            1.0, 1.0, -1.0, 2.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 4.0, 0.0,
            0.0, 16.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0,
        ],
    )
    .transpose();
    // Transpose feito porque a matriz é lida por colunas

    println!("Matriz A: {m1}");

    let b1 = DVector::from_vec(vec![2.0, 2.0, 2.0, 20.0, 4.0]);

    println!("Vetor B: {b1}");

    run_gauss(m1.clone(), b1.clone(), 1e-6);
    run_gauss_seidel(m1.clone(), b1.clone(), 1e-3, 10);

    Ok(())
}
