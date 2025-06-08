use nalgebra::{DMatrix, DVector};

use crate::gaussian_elimination;

pub fn calcula_resultado(coeficientes: DVector<f64>, x: f64) -> f64 {
    let mut res = 0.0;

    for (i, c) in coeficientes.iter().enumerate() {
        res += x.powi(i as i32) * c;
    }

    res
}

/// Função que recebe pontos (x, y) e retorna coeficientes do polinômio interpolador
pub fn interpolacao_sistema_linear(xs: &[f64], ys: &[f64]) -> Result<impl Fn(f64) -> f64, String> {
    let n = xs.len();
    assert_eq!(n, ys.len(), "xs e ys devem ter o mesmo tamanho");

    // Montar matriz Vandermonde (n x n)
    let mut vander = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            vander[(i, j)] = xs[i].powi(j as i32);
        }
    }

    // Vetor y
    let y = DVector::from_vec(ys.to_vec());

    // Resolver V * a = y
    let a = gaussian_elimination(vander, y, 10e-6)?;
    println!("A = {}", a);

    Ok(move |x: f64| {
        let mut res = 0.0;

        for (i, c) in a.iter().enumerate() {
            res += x.powi(i as i32) * c;
        }

        res
    })
}


pub fn intepolacao_lagrange(xs: &[f64], ys: &[f64]) -> impl Fn(f64) -> f64 {
    let n = xs.len();
    assert_eq!(n, ys.len(), "xs e ys devem ter o mesmo tamanho");

    move |x| {
        let mut result = 0.0;

        for i in 0..n {
            // Calcular o polinômio base L_i(x)
            let mut li = 1.0;
            for j in 0..n {
                if j != i {
                    li *= (x - xs[j]) / (xs[i] - xs[j]);
                }
            }
            result += ys[i] * li;
        }

        result
    }
}
/// Calcula os coeficientes do polinômio de Newton via diferenças divididas
pub fn newton_diferencas_divididas(xs: &[f64], ys: &[f64]) -> impl Fn(f64) -> f64 {
    let n = xs.len();
    assert_eq!(n, ys.len(), "xs e ys devem ter o mesmo tamanho");

    // Inicializa uma cópia do vetor ys para ser modificado em-place
    let mut coef = ys.to_vec();

    // Construção da tabela triangular (diferenças divididas)
    for j in 1..n {
        for i in (j..n).rev() {
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j]);
        }
    }

    move |x| {
        let n = coef.len();
        let mut result = coef[n - 1];

        for i in (0..n - 1).rev() {
            result = result * (x - xs[i]) + coef[i];
        }

        result
    }
}
