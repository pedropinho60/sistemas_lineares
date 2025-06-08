use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> Result<(), String> {
    let xs2 = vec![1990.0, 2000.0, 2010.0];
    let ys2 = vec![335.0, 370.0, 388.0];

    let xs3 = vec![1980.0, 1990.0, 2000.0, 2010.0];
    let ys3 = vec![337.0, 335.0, 370.0, 388.0];

    let (coefs, p2) = interpolacao_sistema_linear(&xs2, &ys2)?;
    let (_, p3) = interpolacao_sistema_linear(&xs3, &ys3)?;
    println!("Concentração real em 2008: 381ppm");
    println!("Concentração estimada em 2008 (2° grau): {}ppm", p2(2008.0));
    println!("Concentração estimada em 2008 (3° grau): {}ppm", p3(2008.0));

    let a = coefs[2];
    let b = coefs[1];
    let c = coefs[0] - 350.0;

    // ax^2 + bx + c - 350 = 0

    let delta = b.powi(2) - 4.0 * a * c;

    let x1 = (-b + delta.sqrt()) / (2.0 * a);
    let x2 = (-b - delta.sqrt()) / (2.0 * a);

    println!("Anos em que f(x) = 350: {} e {}", x1 as i32, x2 as i32);

    Ok(())
}
