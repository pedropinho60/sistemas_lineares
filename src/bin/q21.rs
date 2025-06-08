use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> Result<(), String> {
    let xs1 = vec![0.6, 0.8];
    let ys1 = vec![2.0333, 2.6965];

    let xs2 = vec![0.6, 0.8, 1.0];
    let ys2 = vec![2.0333, 2.6965, 3.7183];

    let p1 = intepolacao_lagrange(&xs1, &ys1);
    let p2 = intepolacao_lagrange(&xs2, &ys2);
    println!("f(0.75) estimado com grau 1: {}", p1(0.75));
    println!("f(0.75) estimado com grau 2: {}", p2(0.75));
    Ok(())
}
