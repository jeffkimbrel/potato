# Potato Validation Summary

All 10 potatoes now pass validation with standardized schema.

## Fixes Applied

### 1. Schema Standardization
- Replaced `is_marker` with `marker` in:
  - bhac.json (5 instances)
  - test_glycolysis.json (3 instances)
  - nitrogen_fixation.json (removed - see below)

### 2. Removed Duplicate Pathway
- **Deleted:** nitrogen_fixation.json
- **Reason:** Superseded by nitrogen_fixation_mo.json and nitrogen_fixation_v.json
- The original had sequential steps but Mo and V systems are alternative pathways

### 3. Fixed Input/Output Targets
- **nitrogen_fixation_mo.json:** Fixed targets to reference actual node `nifHDK_1`
- **nitrogen_fixation_v.json:** Fixed targets to reference actual node `vnfDKGH_1`

### 4. Added Input/Output Specifications
All potatoes now have input/output specs:

| Potato | Input | Output |
|--------|-------|--------|
| bhac | glyoxylate (C00048) | aspartate (C00049) |
| entner_doudoroff_classic | D-glucose (C00031) | pyruvate (C00022) |
| entner_doudoroff_np | D-gluconate (C00257) | glyceraldehyde-3P (C00118) |
| entner_doudoroff_semi_phos | D-gluconate (C00257) | glycerate-3P (C00197) |
| entner_doudoroff_semi_phos_alt | D-gluconate (C00257) | glyceraldehyde-3P (C00118) |
| glyoxylate_cycle | acetyl-CoA (C00024) | oxaloacetate (C00036) |
| microcystin_degradation | microcystin-LR (C05371) | ADDA + peptides |
| nitrogen_fixation_mo | N2 (C00697) | 2 NH3 (C00014) |
| nitrogen_fixation_v | N2 (C00697) | 2 NH3 (C00014) |
| test_glycolysis | glyceraldehyde-3P (C00118) | 2-phosphoglycerate (C00631) |

## Current Status

✅ All 10 potatoes pass `validate_potato()` with no errors or warnings
✅ Standardized schema (marker field)
✅ Complete input/output specifications
✅ Correct node references in input/output targets
✅ `verified: false` on all potatoes (ready for biological review)

## Next Steps (Biological Accuracy Review)

Each pathway needs individual biological accuracy review:

1. **bhac** - Beta-hydroxyaspartate cycle
2. **entner_doudoroff_classic** - Classic ED pathway
3. **entner_doudoroff_np** - Non-phosphorylative ED
4. **entner_doudoroff_semi_phos** - Semi-phosphorylative ED
5. **entner_doudoroff_semi_phos_alt** - Alt semi-phosphorylative ED
6. **glyoxylate_cycle** - Glyoxylate shunt
7. **microcystin_degradation** - mlrABCD pathway
8. **nitrogen_fixation_mo** - Mo-nitrogenase
9. **nitrogen_fixation_v** - V-nitrogenase
10. **test_glycolysis** - Test pathway (3 steps)

For each, verify:
- Gene assignments are correct
- OR branches represent true biochemical alternatives
- Required vs. optional designations are accurate
- Marker genes are appropriately identified
- Detection methods (KO/BLAST/HMM) are comprehensive
- Edge connectivity reflects actual biochemistry
- min_fraction thresholds are reasonable
