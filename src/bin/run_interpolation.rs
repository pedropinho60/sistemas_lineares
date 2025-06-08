use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> Result<(), String> {
    // Função usada: 1 / (1 + 25x²)
    let f = |x: f64| 1.0 / (1.0 + 25.0 * x.powi(2));

    // let n = 5;
    // let xs = vec![-1.0 , -0.6, -0.2, 0.2, 0.6, 1.0];
    // let ys = vec![0.03846154, 0.1, 0.5, 0.5, 0.1, 0.03846154];

    let n = 10;
    let xs = vec![-1.0, -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.];
    let ys = vec![
        0.03846154, 0.05882353, 0.1, 0.2, 0.5, 1.0, 0.5, 0.2, 0.1, 0.05882353, 0.03846154,
    ];

    // let n = 15;
    // let xs = vec![
    //     -1.0, -0.86666667, -0.73333333, -0.6, -0.46666667,
    //     -0.33333333, -0.2, -0.06666667, 0.06666667, 0.2,
    //     0.33333333, 0.46666667, 0.6, 0.73333333, 0.86666667, 1.0
    // ];
    //
    // let ys = vec![
    //     0.03846154, 0.0505618, 0.06923077, 0.1, 0.15517241,
    //     0.26470588, 0.5, 0.9, 0.9, 0.5,
    //     0.26470588, 0.15517241, 0.1, 0.06923077, 0.0505618, 0.03846154
    // ];

    // ---------------- Test Vectors ------------------------------
    let xtest = (-20..=20).map(|x| x as f64 / 20.0).collect::<Vec<_>>();

    let ytest = xtest.iter().map(|&x| f(x)).collect::<Vec<_>>();

    let (_, p1) = interpolacao_sistema_linear(&xs, &ys)?;
    // let p2 = intepolacao_lagrange(&xs, &ys);
    // let p3 = newton_diferencas_divididas(&xs, &ys);

    println!("n = {}", n);
    println!("{:<8} {:<12} {:<12}", "x", "y", "p(x)");
    println!("{:-<34}", "");

    for i in 0..xtest.len() {
        println!(
            "{:<8} {:<12.6} {:<12.6}",
            format!("{: >+6.2}", xtest[i]),
            ytest[i],
            p1(xtest[i])
        );
    }

    Ok(())
}
