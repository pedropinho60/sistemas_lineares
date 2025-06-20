use nalgebra::{DMatrix, DVector};
use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> std::io::Result<()> {
    let m1 = DMatrix::from_vec(
        4,
        4,
        vec![
            10.0, -1.0, 2.0, 0.0, -1.0, 11.0, -1.0, 3.0, 2.0, -1.0, 10.0, -1.0, 0.0, 3.0, -1.0, 8.0,
        ],
    )
    .transpose();
    // Transpose feito porque a matriz é lida por colunas

    println!("Matriz A: {m1}");

    let b1 = DVector::from_vec(vec![6.0, 25.0, -11.0, 15.0]);

    println!("Vetor B: {b1}");

    run_gauss_jacobi(m1.clone(), b1.clone(), 1e-3, 5);
    run_gauss_seidel(m1.clone(), b1.clone(), 1e-3, 5);

    Ok(())
}
