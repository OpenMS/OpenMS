use candle_core::{DType, Device};
use candle_nn::{VarBuilder, VarMap};
use redeem_topaz::building_blocks::conv_encoder::ConvBranch;
use redeem_topaz::building_blocks::coelution::CoelutionHead;
use redeem_topaz::building_blocks::trace_input::TraceInputMode;
fn count(vm: &VarMap) -> usize { vm.all_vars().iter().map(|v| v.elem_count()).sum() }
fn main() -> candle_core::Result<()> {
    for c in [6usize, 3] {
        let vm = VarMap::new();
        let vb = VarBuilder::from_varmap(&vm, DType::F32, &Device::Cpu);
        let _b = ConvBranch::new(vb.pp("b"), c, 32, 64, TraceInputMode::Dual, false, 16)?;
        let vm2 = VarMap::new();
        let vb2 = VarBuilder::from_varmap(&vm2, DType::F32, &Device::Cpu);
        let _h = CoelutionHead::new(vb2.pp("h"), c, 10.0, 2, 8)?;
        println!("branch_c{c}\t{}\tcoelution_c{c}\t{}", count(&vm), count(&vm2));
    }
    Ok(())
}
