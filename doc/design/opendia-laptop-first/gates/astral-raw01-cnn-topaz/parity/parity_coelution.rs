// Numeric parity harness: TOPAZ coelution + MS1/MS2 agreement on fixed inputs (zero-weight learned parts).
use candle_core::{DType, Device, Tensor};
use candle_nn::VarBuilder;
use redeem_topaz::building_blocks::coelution::CoelutionHead;
use redeem_topaz::building_blocks::conv_encoder::TraceEncoder;
use redeem_topaz::model::topaz::TopazConfig;

fn input(n: usize, c: usize, l: usize, seed: u64) -> Vec<f32> {
    let mut s = seed;
    (0..n * c * l)
        .map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 33) % 1000) as f32
        })
        .collect()
}

fn main() -> candle_core::Result<()> {
    let dev = Device::Cpu;
    let vb = VarBuilder::zeros(DType::F32, &dev);
    let (n, l) = (3usize, 32usize);
    let ms2 = Tensor::from_vec(input(n, 6, l, 7), (n, 6, l), &dev)?;
    let ms1 = Tensor::from_vec(input(n, 3, l, 11), (n, 3, l), &dev)?;
    let head = CoelutionHead::new(vb.pp("c6"), 6, 10.0, 2, 8)?;
    let head3 = CoelutionHead::new(vb.pp("c3"), 3, 10.0, 2, 8)?;
    let cfg = TopazConfig { ms1_cmax: 3, ms2_cmax: 6, l, use_coelution_head: true, ..Default::default() };
    let enc = TraceEncoder::new(vb.pp("enc"), &cfg)?;
    let x = Tensor::cat(&[ms1.clone(), ms2.clone()], 1)?;
    let (_emb, _coe, ms12) = enc.forward_components(&x)?;
    for (name, t) in [("coe6", head.forward(&ms2)?), ("coe3", head3.forward(&ms1)?), ("ms12", ms12)] {
        for row in t.to_vec2::<f32>()? {
            println!("{name}\t{}", row.iter().map(|v| format!("{v:.6}")).collect::<Vec<_>>().join("\t"));
        }
    }
    Ok(())
}
