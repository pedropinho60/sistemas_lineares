use nalgebra::{DMatrix, DVector};

pub fn gaussian_elimination(mut a: DMatrix<f64>, mut b: DVector<f64>, eps: f64) -> DVector<f64> {
    let n = a.nrows();

    for k in 0..n {
        if a[(k, k)].abs() < eps {
            let mut max_linha = k;
            for i in k + 1..n {
                if a[(i, k)].abs() > a[(max_linha, k)] {
                    max_linha = i;
                }
            }
            a.swap_rows(k, max_linha);
            b.swap_rows(k, max_linha);
        }

        for i in k + 1..n {
            let mult = a[(i, k)] / a[(k, k)];
            for j in 0..n {
                a[(i, j)] -= mult * a[(k, j)];
            }
            b[i] -= mult * b[k];
        }
    }

    let mut x = DVector::zeros(n);

    for i in (0..n).rev() {
        let mut soma = 0.0;
        for j in i + 1..n {
            soma += a[(i, j)] * x[j];
        }
        x[i] = (b[i] - soma) / a[(i, i)];
    }

    x
}

pub fn lower_upper_decomposition(a: DMatrix<f64>, mut b: DVector<f64>, eps: f64) -> DVector<f64> {
    let n = a.nrows();

    let mut l = DMatrix::identity(n, n);

    let mut u = a.clone();

    for k in 0..n {
        if u[(0, 0)].abs() < eps {
            let mut max_linha = k;
            for i in k + 1..n {
                if u[(i, k)].abs() > u[(max_linha, k)] {
                    max_linha = i;
                }
            }
            u.swap_rows(k, max_linha);
            b.swap_rows(k, max_linha);
        }
        for i in k + 1..n {
            let mult = u[(i, k)] / u[(k, k)];
            l[(i, k)] = mult;
            for j in 0..n {
                u[(i, j)] -= mult * u[(k, j)];
            }
        }
    }

    let mut y = DVector::zeros(n);

    for i in 0..n {
        let mut soma = 0.0;
        for j in 0..i {
            soma += l[(i, j)] * y[j];
        }
        y[i] = b[i] - soma;
    }

    let mut x = DVector::zeros(n);

    for i in (0..n).rev() {
        let mut soma = 0.0;
        for j in i + 1..n {
            soma += u[(i, j)] * x[j];
        }
        x[i] = (y[i] - soma) / u[(i, i)];
    }

    x
}

pub fn gauss_jacobi(
    a: DMatrix<f64>,
    b: DVector<f64>,
    eps: f64,
    lim_iter: u64,
) -> Result<DVector<f64>, String> {
    let n = a.nrows();
    let mut x_ant = DVector::zeros(n);

    for i in 0..n {
        x_ant[i] = b[i] / a[(i, i)];
    }

    let mut x = x_ant.clone();

    for _ in 0..lim_iter {
        for i in 0..n {
            let mut soma = 0.0;
            for j in 0..n {
                if j != i {
                    soma += a[(i, j)] * x_ant[j];
                }
            }
            x[i] = (b[i] - soma) / a[(i, i)];
        }

        let dist = x
            .iter()
            .zip(x_ant.iter())
            .map(|(a, b)| (a - b).abs())
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();

        let x_max = x
            .iter()
            .map(|a| a.abs())
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();

        if dist / x_max < eps {
            return Ok(x);
        }

        x_ant = x.clone();
    }

    Err(String::from("Limite de iterações atingido"))
}

pub fn gauss_seidel(
    a: DMatrix<f64>,
    b: DVector<f64>,
    eps: f64,
    lim_iter: u64,
) -> Result<DVector<f64>, String> {
    let n = a.nrows();
    let mut x_ant = DVector::zeros(n);

    for i in 0..n {
        x_ant[i] = b[i] / a[(i, i)];
    }

    let mut x = x_ant.clone();

    for _ in 0..lim_iter {
        for i in 0..n {
            let mut soma = 0.0;
            for j in 0..i {
                soma += a[(i, j)] * x[j];
            }
            for j in i + 1..n {
                soma += a[(i, j)] * x_ant[j];
            }
            x[i] = (b[i] - soma) / a[(i, i)];
        }

        let dist = x
            .iter()
            .zip(x_ant.iter())
            .map(|(a, b)| (a - b).abs())
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();

        let x_max = x
            .iter()
            .map(|a| a.abs())
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();

        if dist / x_max < eps {
            return Ok(x);
        }

        x_ant = x.clone();
    }

    Err(String::from("Limite de iterações atingido"))
}
