use sistemas_lineares_e_interpolacao_polinomial::*;

fn main() -> Result<(), String> {
    let xs = vec![0.0, 10.0, 20.0, 30.0];
    let ys = vec![0.0, 20.56, 30.67, 67.78];

    let p = newton_diferencas_divididas(&xs, &ys);
    println!("Distância após 15.6 minutos: {}", p(15.6));

    Ok(())
}
