use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> Result<(), String> {
    let xs = vec![2008.0, 2011.0, 2014.0];
    let ys = vec![323.0, 1430.0, 7310.0];

    let (coefs, p) = interpolacao_sistema_linear(&xs, &ys)?;
    println!("Capacidade estimada em 2010: {}MW", p(2010.0));
    println!("Capacidade estimada em 2015: {}MW", p(2015.0));

    let a = coefs[2];
    let b = coefs[1];
    let c = coefs[0] - 5000.0;

    // ax^2 + bx + c - 5000 = 0

    let delta = b.powi(2) - 4.0 * a * c;

    let x1 = (-b + delta.sqrt()) / (2.0 * a);
    let x2 = (-b - delta.sqrt()) / (2.0 * a);

    println!("Anos em que f(x) = 5000: {} e {}", x1 as i32, x2 as i32);

    Ok(())
}
