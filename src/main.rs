fn main() {
    let a = vec![
        vec![10.0, 2.0, 1.0],
        vec![1.0, 5.0, 1.0],
        vec![2.0, 3.0, 10.0],
    ];

    let b = vec![7.0, -8.0, 6.0];

    let x = gaussian_elimination(a.clone(), b.clone());
    let x2 = lower_upper_decomposition(a.clone(), b.clone());
    let x3 = gauss_jacobi(a.clone(), b.clone());
    let x4 = gauss_seidel(a.clone(), b.clone());

    println!("{:?}", x);
    println!("{:?}", x2);
    println!("{:?}", x3);
    println!("{:?}", x4);
}

fn gaussian_elimination(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = a.len();

    for k in 0..n {
        if a[k][k] == 0.0 {
            let mut max_linha = k;
            for i in k + 1..n {
                if a[i][k].abs() > a[max_linha][k] {
                    max_linha = i;
                }
            }
            a.swap(k, max_linha);
            b.swap(k, max_linha);
        }

        for i in k + 1..n {
            let mult = a[i][k] / a[k][k];
            for j in 0..n {
                a[i][j] -= mult * a[k][j];
            }
            b[i] -= mult * b[k];
        }
    }

    let mut x = vec![0.0; n];

    for i in (0..n).rev() {
        let mut soma = 0.0;
        for j in i + 1..n {
            soma += a[i][j] * x[j];
        }
        x[i] = (b[i] - soma) / a[i][i];
    }

    x
}

fn lower_upper_decomposition(a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = a.len();

    let mut l = vec![vec![0.0; n]; n];

    for k in 0..n {
        l[k][k] = 1.0;
    }

    let mut u = a.clone();

    for k in 0..n {
        if u[0][0] == 0.0 {
            let mut max_linha = k;
            for i in k + 1..n {
                if u[i][k].abs() > u[max_linha][k] {
                    max_linha = i;
                }
            }
            u.swap(k, max_linha);
            b.swap(k, max_linha);
        }
        for i in k + 1..n {
            let mult = u[i][k] / u[k][k];
            l[i][k] = mult;
            for j in 0..n {
                u[i][j] -= mult * u[k][j];
            }
        }
    }

    let mut y = vec![0.0; n];

    for i in 0..n {
        let mut soma = 0.0;
        for j in 0..i {
            soma += l[i][j] * y[j];
        }
        y[i] = b[i] - soma;
    }

    let mut x = vec![0.0; n];

    for i in (0..n).rev() {
        let mut soma = 0.0;
        for j in i + 1..n {
            soma += u[i][j] * x[j];
        }
        x[i] = (y[i] - soma) / u[i][i];
    }

    x
}

fn gauss_jacobi(a: Vec<Vec<f64>>, b: Vec<f64>) -> Result<Vec<f64>, String> {
    let n = a.len();
    let mut x_ant = vec![0.0; n];

    for i in 0..n {
        x_ant[i] = b[i] / a[i][i];
    }

    let mut x = x_ant.clone();

    let prec = 1e-8;

    for _ in 0..500 {
        for i in 0..n {
            let mut soma = 0.0;
            for j in 0..n {
                if j != i {
                    soma += a[i][j] * x_ant[j];
                }
            }
            x[i] = (b[i] - soma) / a[i][i];
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

        if dist / x_max < prec {
            return Ok(x);
        }

        x_ant = x.clone();
    }

    Err(String::from("Limite de iterações atingido"))
}

fn gauss_seidel(a: Vec<Vec<f64>>, b: Vec<f64>) -> Result<Vec<f64>, String> {
    let n = a.len();
    let mut x_ant = vec![0.0; n];

    for i in 0..n {
        x_ant[i] = b[i] / a[i][i];
    }

    let mut x = x_ant.clone();

    let prec = 1e-8;

    for _ in 0..500 {
        for i in 0..n {
            let mut soma = 0.0;
            for j in 0..i {
                soma += a[i][j] * x[j];
            }
            for j in i + 1..n {
                soma += a[i][j] * x_ant[j];
            }
            x[i] = (b[i] - soma) / a[i][i];
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

        if dist / x_max < prec {
            return Ok(x);
        }

        x_ant = x.clone();
    }

    Err(String::from("Limite de iterações atingido"))
}
