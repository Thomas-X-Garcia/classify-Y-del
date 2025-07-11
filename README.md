# classify-Y-del: Y-Chromosomal Microdeletion Classifier

A Python implementation of the EAA/EMQN best practice guidelines for molecular diagnosis of Y-chromosomal microdeletions (Krausz et al., 2023).

**Author:** Thomas X. Garcia, PhD, HCLD  
**Version:** 1.0  
**License:** MIT

## Table of Contents
- [Overview](#overview)
- [Background](#background)
- [Installation](#installation)
- [Usage](#usage)
- [Input Format](#input-format)
- [Output Classifications](#output-classifications)
- [Clinical Significance](#clinical-significance)
- [Marker Requirements](#marker-requirements)
- [Examples](#examples)
- [Citation](#citation)
- [References](#references)

## Overview

`classify-Y-del.py` is a clinical-grade Python script that classifies Y-chromosomal microdeletions according to the 2023 European Academy of Andrology (EAA) and European Molecular Genetics Quality Network (EMQN) best practice guidelines. The script provides accurate classification of AZF (AZoospermia Factor) deletions, which are the second most common genetic cause of male infertility after Klinefelter syndrome.

### Key Features
- Full compliance with EAA/EMQN 2023 guidelines
- Two-step diagnostic process (basic + extension analysis)
- Comprehensive validation of control markers
- TESE (Testicular Sperm Extraction) prognosis for clinical decision-making
- Detection of partial deletions including gr/gr
- Robust error handling and quality control checks
- Detailed clinical recommendations in verbose mode

## Background

### Y-Chromosomal Microdeletions

Y-chromosomal microdeletions occur in approximately 1 in 4,000 men in the general population, but their frequency is significantly increased in infertile men:
- **Azoospermic men**: 2-10% (or higher)
- **Severely oligozoospermic men** (<2 × 10⁶ sperm/mL): Lower frequency than azoospermic men

The deletions occur through Non-Allelic Homologous Recombination (NAHR) between repetitive sequences in the AZF regions of the Y chromosome.

### AZF Regions

Three regions on the Y chromosome long arm are critical for spermatogenesis:

1. **AZFa** (792 kb)
   - Contains: USP9Y and DDX3Y genes
   - Deletion mechanism: NAHR between HERVyq1 and HERVyq2 retroviral sequences
   - Phenotype: Sertoli Cell-Only (SCO), azoospermia

2. **AZFb** (6.2 Mb)
   - Contains: 32 transcription units
   - Deletion mechanism: NAHR between P5/proximal P1 palindromes
   - Phenotype: Azoospermia with SCO or spermatogenic arrest

3. **AZFc** (3.5 Mb)
   - Contains: 21 transcription units including DAZ gene family
   - Deletion mechanism: NAHR between b2/b4 amplicons
   - Phenotype: Azoospermia or severe oligozoospermia

### Clinical Relevance

Testing for AZF deletions is indicated for:
- All azoospermic men
- Men with severe oligozoospermia (<2 × 10⁶ sperm/mL)
- Pre-TESE evaluation
- Before varicocelectomy in azoo/oligozoospermic men

## Installation

### Requirements
- Python 3.6 or higher
- No external dependencies required (uses only Python standard library)

### Download
```bash
git clone https://github.com/Thomas-X-Garcia/classify-Y-del.git
cd classify-Y-del
chmod +x classify-Y-del.py
```

## Usage

### Basic Usage
```bash
python3 classify-Y-del.py <marker_file.tsv>
```

### Verbose Mode (Recommended for Clinical Use)
```bash
python3 classify-Y-del.py <marker_file.tsv> --verbose
```

### Validate Marker Completeness
```bash
python3 classify-Y-del.py <marker_file.tsv> --validate-only
```

### Command Line Options
- `-v, --verbose`: Generate detailed report with clinical recommendations
- `--validate-only`: Check if all required markers are present without classification
- `-h, --help`: Show help message with marker requirements

## Input Format

The script expects a tab-separated values (TSV) file with two columns:
1. **marker_name**: The STS marker identifier
2. **status**: Either "present" or "absent" (case-insensitive)

### Example Input File
```tsv
marker	status
sY14	present
ZFX/ZFY	present
sY84	absent
sY86	absent
sY127	present
sY134	present
sY254	present
sY255	present
```

**Note:** The script automatically detects and skips header rows.

## Output Classifications

### Complete Deletions

#### AZFa Deletions
- **COMPLETE_AZFA_DELETION**: Complete loss of AZFa region
  - TESE prognosis: Virtually impossible
  - All cases result in azoospermia with SCO

#### AZFb Deletions
- **COMPLETE_AZFB_DELETION_P5/PROXIMAL_P1**: Complete P5/proximal P1 deletion
  - TESE prognosis: Virtually impossible
  - Key marker: sY1192 absent
- **PARTIAL_AZFB_DELETION**: Retention of sY1192
  - TESE prognosis: Possible
  - May have residual spermatogenesis

#### AZFc Deletions
- **COMPLETE_AZFC_DELETION_B2/B4**: Non-terminal complete deletion
  - TESE prognosis: Approximately 50% success rate
  - sY160 present
- **TERMINAL_AZFC_DELETION**: Terminal deletion
  - sY160 absent
  - Requires karyotype analysis for 46,XY/45,X mosaicism

#### Combined Deletions
- **COMPLETE_AZFBC_DELETION_P5/DISTAL_P1**: Larger deletion
  - TESE prognosis: Virtually impossible
- **COMPLETE_AZFBC_DELETION_P4/DISTAL_P1**: Smaller deletion
  - TESE prognosis: May be positive
- **COMPLETE_AZFABC_DELETION**: All three regions deleted
  - Often associated with karyotype abnormalities

### Partial Deletions
- **PARTIAL_AZFC_GR/GR_DELETION**: gr/gr partial deletion
  - Population-specific risk factor
  - Increased risk for testicular germ cell tumors
  - Compatible with varying sperm counts

### Special Cases
- **46,XX_MALE_OR_COMPLETE_Y_CHROMOSOME_ABSENCE**: sY14 absent, ZFX/ZFY present
- **METHODOLOGICAL_ERROR**: sY254/sY255 discordance
- **NO_DELETION_DETECTED**: All markers present as expected

## Clinical Significance

### TESE Prognosis by Deletion Type

| Deletion Type | TESE Success Rate | Clinical Recommendation |
|--------------|-------------------|------------------------|
| Complete AZFa | 0% | TESE not recommended |
| Complete AZFb (P5/proximal P1) | ~0% | TESE not recommended |
| Partial AZFb (sY1192+) | Possible | TESE may be attempted |
| Complete AZFc (b2/b4) | ~50% | TESE recommended with counseling |
| Complete AZFbc (P5/distal P1) | ~0% | TESE not recommended |
| Complete AZFbc (P4/distal P1) | Possible | TESE may be attempted |

### Genetic Counseling Considerations

1. **Inheritance**: All Y-chromosomal deletions are transmitted to male offspring
2. **Phenotype prediction**: Son's fertility cannot be precisely predicted due to:
   - Different genetic background
   - Environmental factors
   - Potential deletion expansion
3. **Family screening**: Indicated only for:
   - Complete AZFc deletions
   - Partial AZFa/AZFb deletions
   - gr/gr deletions

## Marker Requirements

### Control Markers (Always Required)
- **sY14**: SRY gene marker (testis determining factor)
- **ZFX/ZFY** (or ZFX/Y): Internal PCR control

### Basic Analysis Markers
According to EAA/EMQN guidelines Appendix B:

| Region | Markers | Expected Result in Deletion |
|--------|---------|----------------------------|
| AZFa | sY84, sY86 | Absent |
| AZFb | sY127, sY134 | Absent |
| AZFc | sY254, sY255 | Absent |

### Extension Analysis Markers
Required for confirmed deletions to determine subtype and prognosis:

#### AZFa Extension
- **Proximal boundary**: sY82 (present), sY1064 (absent)
- **Distal boundary**: sY1065 or sY1182 (absent), sY88 (present)

#### AZFb Extension
- **Proximal boundary**: sY105 (present), sY121 (absent), sY1224 (variable)
- **Distal boundary**: sY1192 (absent for complete), sY153 (present)

#### AZFc Extension
- **sY160**: Present (b2/b4) or Absent (terminal)

#### gr/gr Partial Deletion
- **sY1291**: Absent
- **sY1191**: Present

### Critical Updates in 2023 Guidelines

1. **sY1064 replaces sY83** for AZFa proximal boundary
   - Some complete AZFa deletions retain sY83
   
2. **sY1192 replaces sY143** for AZFb distal boundary
   - Only sY1192 distinguishes complete from partial deletions
   - Critical for TESE prognosis

3. **sY84 primer update** for Asian populations
   - SNP (rs72609647) requires modified reverse primer

## Examples

### Example 1: Complete AZFc Deletion
```bash
python3 classify-Y-del.py examples/complete_azfc.tsv --verbose
```

Output:
```
=== Y-CHROMOSOMAL MICRODELETION ANALYSIS REPORT ===

CLASSIFICATION: COMPLETE_AZFC_DELETION_B2/B4 (TESE: approximately_50_percent)

MARKER SUMMARY:
Control markers:
  sY14: present
  ZFX/ZFY: present

Basic AZF markers:
  AZFa:
    sY84: present
    sY86: present
  AZFb:
    sY127: present
    sY134: present
  AZFc:
    sY254: absent
    sY255: absent

CLINICAL RECOMMENDATIONS:
- TESE may be attempted (approximately 50% success rate)
- Genetic counseling mandatory (deletion will be transmitted to male offspring)
- Consider karyotype analysis to rule out 46,XY/45,X mosaicism
```

### Example 2: Partial AZFb Deletion
Input file with sY1192 present indicates partial deletion compatible with spermatogenesis.

### Example 3: gr/gr Deletion
Population-specific risk factor that may or may not affect fertility.

## Citation

If you use this software in your research or clinical practice, please cite:

```
Garcia, T.X. (2025). classify-Y-del: Y-Chromosomal Microdeletion Classifier v1.0. 
GitHub repository: https://github.com/Thomas-X-Garcia/classify-Y-del
```

And the guidelines paper:
```
Krausz, C., Navarro-Costa, P., Wilke, M., & Tüttelmann, F. (2023). 
EAA/EMQN best practice guidelines for molecular diagnosis of Y-chromosomal 
microdeletions: State of the art 2023. Andrology. 
https://doi.org/10.1111/andr.13514
```

## References

1. Krausz C, Navarro-Costa P, Wilke M, Tüttelmann F. EAA/EMQN best practice guidelines for molecular diagnosis of Y-chromosomal microdeletions: State of the art 2023. Andrology. 2023. doi:10.1111/andr.13514

2. Tiepolo L, Zuffardi O. Localization of factors controlling spermatogenesis in the nonfluorescent portion of the human Y chromosome long arm. Hum Genet. 1976;34(2):119-124.

3. Vogt PH, Edelmann A, Kirsch S, et al. Human Y chromosome azoospermia factors (AZF) mapped to different subregions in Yq11. Hum Mol Genet. 1996;5(7):933-943.

## Support

For questions, bug reports, or feature requests, please open an issue on the [GitHub repository](https://github.com/Thomas-X-Garcia/classify-Y-del).

## Disclaimer

This software is intended for research and clinical use by qualified professionals. Always correlate results with clinical findings and follow local guidelines for genetic testing and counseling.
