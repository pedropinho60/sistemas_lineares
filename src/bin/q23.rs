use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> Result<(), String> {
    let xs = vec![0.0, 6.0, 10.0];
    let ys = vec![6.67, 17.33, 42.67];

    let (coefs, p) = interpolacao_sistema_linear(&xs, &ys)?;
    println!("Peso no dia 7: {}", p(7.0));

    let a = coefs[2];
    let b = coefs[1];
    let c = coefs[0] - 10.0;

    // ax^2 + bx + c - 10 = 0

    let delta = b.powi(2) - 4.0 * a * c;

    let x1 = (-b + delta.sqrt()) / (2.0 * a);
    let x2 = (-b - delta.sqrt()) / (2.0 * a);

    println!("Dias onde f(x) = 10: {} e {}", x1, x2);

    Ok(())
}
