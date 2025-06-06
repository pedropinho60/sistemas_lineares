use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use nalgebra::{DMatrix, DVector};

pub mod interpolation;
pub mod methods;

pub use methods::*;

pub fn from_file(name: &str) -> std::io::Result<DVector<f64>> {
    let file = File::open(format!("matrizes/{}", name))?;
    let reader = BufReader::new(file);

    let mut b1: Vec<f64> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        b1.push(line.parse().unwrap());
    }

    Ok(DVector::from_vec(b1))
}

pub fn residual(a: DMatrix<f64>, b: DVector<f64>, x: DVector<f64>) -> DVector<f64> {
    let b2 = a * x;

    b2.iter()
        .zip(b.iter())
        .map(|(a, b)| (a - b).abs())
        .collect::<Vec<_>>()
        .into()
}

pub fn run_gauss(a: DMatrix<f64>, b: DVector<f64>, eps: f64) {
    let x = gaussian_elimination(a.clone(), b.clone(), eps);
    println!("Eliminação de Gauss: ");
    match x {
        Ok(x) => {
            println!("{}", x);

            let res = residual(a.clone(), b.clone(), x.clone());
            println!("Resíduo: {}", res);
        }

        Err(e) => println!("Erro: {}\n", e),
    }
}

pub fn run_lu(a: DMatrix<f64>, b: DVector<f64>, eps: f64) {
    let x = lower_upper_decomposition(a.clone(), b.clone(), eps);
    println!("Fatoração LU: ");
    match x {
        Ok(x) => {
            println!("{}", x);

            let res = residual(a.clone(), b.clone(), x.clone());
            println!("Resíduo: {}", res);
        }

        Err(e) => println!("Erro: {}\n", e),
    }
}

pub fn run_gauss_jacobi(a: DMatrix<f64>, b: DVector<f64>, eps: f64, lim_iter: u64) {
    let x = gauss_jacobi(a.clone(), b.clone(), eps, lim_iter);
    println!("Gauss-Jacobi: ");
    match x {
        Ok(x) => {
            println!("{}", x);

            let res = residual(a.clone(), b.clone(), x.clone());
            println!("Resíduo: {}", res);
        }

        Err(e) => println!("Erro: {}\n", e),
    }
}

pub fn run_gauss_seidel(a: DMatrix<f64>, b: DVector<f64>, eps: f64, lim_iter: u64) {
    let x = gauss_seidel(a.clone(), b.clone(), eps, lim_iter);

    println!("Gauss-Seidel:");
    match x {
        Ok(x) => {
            println!("{}", x);

            let res = residual(a.clone(), b.clone(), x.clone());
            println!("Resíduo: {}", res);
        }
        Err(e) => println!("Erro: {}\n", e),
    }
}
