use sistemas_lineares_e_interpolacao_polinomial::interpolation::{
    calcula_resultado, intepolacao_lagrange, interpolacao_sistema_linear,
    newton_diferencas_divididas,
};

fn main() -> Result<(), String> {
    // Exemplo: 3 pontos
    let xs = vec![1.0, 2.0, 3.0];
    let ys = vec![2.0, 3.0, 5.0];

    let p1 = interpolacao_sistema_linear(&xs, &ys)?;
    let p2 = intepolacao_lagrange(&xs, &ys);
    let p3 = newton_diferencas_divididas(&xs, &ys);

    let x = 10.0;
    println!("f({x}) = {}", p1(x));
    println!("f({x}) = {}", p2(x));
    println!("f({x}) = {}", p3(x));

    Ok(())
}
