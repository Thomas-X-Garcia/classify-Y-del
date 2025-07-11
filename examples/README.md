# Example Files for classify-Y-del

This directory contains example input files demonstrating various Y-chromosomal deletion patterns.

## File Descriptions

### no_deletion.tsv
- Normal male pattern with all markers present
- Expected output: `NO_DELETION_DETECTED`

### complete_azfc.tsv
- Complete AZFc deletion (b2/b4 type)
- sY254 and sY255 absent, sY160 present
- Expected output: `COMPLETE_AZFC_DELETION_B2/B4`
- TESE prognosis: ~50% success rate

### complete_azfb.tsv
- Complete AZFb deletion (P5/proximal P1)
- sY127, sY134, sY121, and sY1192 absent
- Expected output: `COMPLETE_AZFB_DELETION_P5/PROXIMAL_P1`
- TESE prognosis: Virtually impossible

### partial_azfb.tsv
- Partial AZFb deletion with sY1192 retained
- sY127 and sY134 absent, but sY1192 present
- Expected output: `PARTIAL_AZFB_DELETION`
- TESE prognosis: Possible

### grgr_deletion.tsv
- gr/gr partial AZFc deletion
- sY1291 absent, sY1191 present
- Expected output: `PARTIAL_AZFC_GR/GR_DELETION`
- Clinical significance: Population-specific risk factor

### 46xx_male.tsv
- 46,XX male pattern
- All Y-specific markers absent, ZFX/ZFY present
- Expected output: `46,XX_MALE_OR_COMPLETE_Y_CHROMOSOME_ABSENCE`
- TESE prognosis: Not possible

## Usage

To test any example file:

```bash
# Basic classification
python3 ../classify-Y-del.py no_deletion.tsv

# Verbose mode with clinical recommendations
python3 ../classify-Y-del.py complete_azfc.tsv --verbose

# Validate marker completeness
python3 ../classify-Y-del.py partial_azfb.tsv --validate-only
```

## Notes

- All example files contain the complete set of recommended markers from the EAA/EMQN 2023 guidelines
- Files use tab-separated format with header row
- Marker status is case-insensitive ("present" or "absent")