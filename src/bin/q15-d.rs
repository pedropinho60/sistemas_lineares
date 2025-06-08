use nalgebra::{DMatrix, DVector};
use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> std::io::Result<()> {
    let m1 = DMatrix::from_vec(
        4,
        4,
        vec![
            3.0, -2.0, 5.0, 1.0, -6.0, 4.0, -8.0, 1.0, 9.0, -6.0, 19.0, 1.0, 6.0, -4.0, -6.0, 15.0,
        ],
    )
    .transpose();
    // Transpose feito porque a matriz é lida por colunas

    println!("Matriz A: {m1}");

    let b1 = DVector::from_vec(vec![7.0, -9.0, 23.0, 11.0]);

    println!("Vetor B: {b1}");

    run_gauss(m1.clone(), b1.clone(), 1e-6);
    run_lu(m1.clone(), b1.clone(), 1e-6);

    Ok(())
}
