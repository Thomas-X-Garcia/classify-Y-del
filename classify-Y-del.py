#!/usr/bin/env python3
"""
Y-Chromosomal Microdeletion Classifier

This script implements the EAA/EMQN best practice guidelines (2023) for molecular diagnosis 
of Y-chromosomal microdeletions as described in Krausz et al., 2023.

Key features:
- Full compliance with Appendix B marker requirements
- Two-step analysis process (basic + extension)
- Comprehensive error checking and validation
- Handling of partial deletions and special cases
- Detailed and guideline-compliant output messages
"""

import sys
import csv
import argparse
from typing import Dict, List, Tuple, Optional

# Define marker sets according to EAA/EMQN 2023 guidelines (Appendix B)
CONTROL_MARKERS = ['sY14', 'ZFX/ZFY']  # Also written as ZFX/Y in some contexts

# Basic analysis markers (first step)
BASIC_MARKERS = {
    'AZFa': ['sY84', 'sY86'],
    'AZFb': ['sY127', 'sY134'],
    'AZFc': ['sY254', 'sY255']
}

# Extension analysis markers (second step) - for confirmed deletions
EXTENSION_MARKERS = {
    'AZFa': {
        'proximal': {'present': ['sY82'], 'absent': ['sY1064']},
        'distal': {'present': ['sY88'], 'absent': ['sY1065', 'sY1182']}  # either sY1065 or sY1182
    },
    'AZFb': {
        'proximal': {'present': ['sY105'], 'absent': ['sY121'], 'variable': ['sY1224']},
        'distal': {'present': ['sY153'], 'absent': ['sY1192']}
    },
    'AZFc': {
        'terminal': {'sY160': 'absent'},  # terminal deletion
        'non_terminal': {'sY160': 'present'}  # b2/b4 deletion
    },
    'gr/gr': {
        'markers': {'sY1291': 'absent', 'sY1191': 'present'}
    }
}

# All markers that should be tested according to guidelines
ALL_REQUIRED_MARKERS = (
    CONTROL_MARKERS + 
    list(BASIC_MARKERS['AZFa']) + 
    list(BASIC_MARKERS['AZFb']) + 
    list(BASIC_MARKERS['AZFc']) +
    ['sY82', 'sY1064', 'sY1065', 'sY1182', 'sY88',  # AZFa extension
     'sY105', 'sY121', 'sY1224', 'sY1192', 'sY153',  # AZFb extension
     'sY160',  # AZFc extension
     'sY1291', 'sY1191']  # gr/gr partial
)

# Markers that SHOULD be present in healthy males (according to Appendix B expected results)
EXPECTED_PRESENT_MARKERS = ['sY14', 'sY82', 'sY88', 'sY105', 'sY153', 'sY160', 'sY1191', 'ZFX/ZFY']

# Markers that SHOULD be absent in complete deletions (according to the paper)
EXPECTED_ABSENT_IN_DELETIONS = {
    'AZFa': ['sY84', 'sY86', 'sY1064', 'sY1065', 'sY1182'],
    'AZFb': ['sY127', 'sY134', 'sY121', 'sY1192'],
    'AZFc': ['sY254', 'sY255'],
    'AZFbc': ['sY127', 'sY134', 'sY254', 'sY255', 'sY121', 'sY1192'],
    'gr/gr': ['sY1291']
}


def parse_marker_file(filepath: str) -> Optional[Dict[str, str]]:
    """
    Parses a two-column TSV file into a dictionary of marker statuses.
    
    Args:
        filepath: Path to the input TSV file
        
    Returns:
        Dictionary mapping marker names to their status ('present'/'absent')
        Returns None if file cannot be read
    """
    markers = {}
    try:
        with open(filepath, 'r', newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            
            # Process first row to check if it's a header
            try:
                first_row = next(reader)
                if len(first_row) >= 2:
                    # Check if second column looks like data (present/absent)
                    if first_row[1].lower().strip() in ['present', 'absent']:
                        # It's data, process it
                        markers[first_row[0].strip()] = first_row[1].lower().strip()
                    # Otherwise skip as header
                    
                # Process remaining rows
                for row in reader:
                    if len(row) >= 2:
                        marker_name = row[0].strip()
                        marker_status = row[1].lower().strip()
                        if marker_status in ['present', 'absent']:
                            markers[marker_name] = marker_status
                        else:
                            print(f"Warning: Invalid status '{marker_status}' for marker {marker_name}", 
                                  file=sys.stderr)
                                  
            except StopIteration:
                print("Error: Empty file", file=sys.stderr)
                return None
                
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error: An error occurred while reading the file: {e}", file=sys.stderr)
        return None
        
    return markers


def check_control_markers(markers: Dict[str, str]) -> Tuple[bool, str]:
    """
    Check control markers (sY14/SRY and ZFX/ZFY) for quality control.
    
    Returns:
        Tuple of (is_valid, message)
    """
    # Check ZFX/ZFY (should always be present as internal control)
    # Handle both ZFX/ZFY and ZFX/Y variants
    zfxy_status = markers.get('ZFX/ZFY') or markers.get('ZFX/Y')
    if not zfxy_status:
        return False, "CRITICAL: ZFX/ZFY control marker not tested - cannot validate results"
    if zfxy_status != 'present':
        return False, "CRITICAL: ZFX/ZFY control marker absent - technical failure or sample issue"
    
    # Check sY14 (SRY)
    if 'sY14' not in markers:
        return False, "CRITICAL: sY14 (SRY) marker not tested - cannot determine sex chromosome status"
    
    return True, "Control markers valid"


def check_basic_deletion(markers: Dict[str, str], region: str) -> bool:
    """Check if basic markers indicate a deletion in the specified region."""
    if region not in BASIC_MARKERS:
        return False
    
    region_markers = BASIC_MARKERS[region]
    # All markers in the region must be tested and absent for a deletion
    for marker in region_markers:
        if marker not in markers:
            return False  # Marker not tested
        if markers[marker] != 'absent':
            return False  # Marker present
    
    return True


def check_extension_markers(markers: Dict[str, str], region: str) -> Dict[str, any]:
    """
    Perform extension analysis for a confirmed deletion.
    
    Returns dictionary with deletion subtype information.
    """
    result = {'complete': False, 'subtype': 'unknown', 'details': {}}
    
    if region == 'AZFa':
        # Check proximal boundary
        if markers.get('sY82') == 'present' and markers.get('sY1064') == 'absent':
            result['details']['proximal'] = 'correct'
        else:
            result['details']['proximal'] = 'atypical'
            
        # Check distal boundary (need either sY1065 or sY1182 absent)
        distal_absent = (markers.get('sY1065') == 'absent' or 
                        markers.get('sY1182') == 'absent')
        if markers.get('sY88') == 'present' and distal_absent:
            result['details']['distal'] = 'correct'
        else:
            result['details']['distal'] = 'atypical'
            
        # Complete deletion if both boundaries are correct
        if (result['details']['proximal'] == 'correct' and 
            result['details']['distal'] == 'correct'):
            result['complete'] = True
            result['subtype'] = 'complete'
        else:
            result['subtype'] = 'partial_or_atypical'
            
    elif region == 'AZFb':
        # Check for sY1192 - critical for TESE prognosis
        has_sY1192 = markers.get('sY1192') == 'present'
        
        # Check boundaries
        if (markers.get('sY105') == 'present' and 
            markers.get('sY121') == 'absent' and
            markers.get('sY1192') == 'absent' and
            markers.get('sY153') == 'present'):
            result['complete'] = True
            result['subtype'] = 'P5/proximal_P1'
            result['details']['TESE_prognosis'] = 'virtually_impossible'
        elif has_sY1192:
            result['complete'] = False
            result['subtype'] = 'partial'
            result['details']['sY1192'] = 'present'
            result['details']['TESE_prognosis'] = 'possible'
        else:
            result['subtype'] = 'atypical'
            
        # Note sY1224 status (can be variable)
        result['details']['sY1224'] = markers.get('sY1224', 'not_tested')
        
    elif region == 'AZFc':
        # Check sY160 to distinguish terminal from non-terminal
        if 'sY160' not in markers:
            result['subtype'] = 'unknown_terminal_status'
        elif markers['sY160'] == 'present':
            result['complete'] = True
            result['subtype'] = 'b2/b4'
            result['details']['terminal'] = False
            result['details']['TESE_prognosis'] = 'approximately_50_percent'
        else:  # sY160 absent
            result['complete'] = True
            result['subtype'] = 'terminal'
            result['details']['terminal'] = True
            result['details']['karyotype_recommended'] = True
            
    return result


def check_azfbc_deletion(markers: Dict[str, str]) -> Tuple[bool, str]:
    """Check for combined AZFbc deletion."""
    # All four markers must be absent
    required_absent = ['sY127', 'sY134', 'sY254', 'sY255']
    
    for marker in required_absent:
        if marker not in markers or markers[marker] != 'absent':
            return False, ""
    
    # Check for P4/P5 distinction using sY116
    if 'sY116' in markers:
        if markers['sY116'] == 'absent':
            subtype = "P5/distal_P1"
            tese = "virtually_impossible"
        else:
            subtype = "P4/distal_P1"
            tese = "may_be_positive"
    else:
        subtype = "undetermined_P4_P5"
        tese = "unknown"
        
    # Check terminal status
    terminal = ""
    if 'sY160' in markers:
        if markers['sY160'] == 'absent':
            terminal = "_TERMINAL"
            
    return True, f"COMPLETE_AZFBC_DELETION_{subtype}{terminal} (TESE: {tese})"


def check_grgr_deletion(markers: Dict[str, str]) -> bool:
    """Check for gr/gr partial AZFc deletion."""
    return (markers.get('sY1291') == 'absent' and 
            markers.get('sY1191') == 'present')


def classify_y_chromosome_state(markers: Dict[str, str]) -> str:
    """
    Main classification function following EAA/EMQN 2023 guidelines.
    
    Args:
        markers: Dictionary of marker statuses
        
    Returns:
        Classification string
    """
    # Step 1: Validate control markers
    control_valid, control_msg = check_control_markers(markers)
    if not control_valid:
        return control_msg
    
    # Check for 46,XX male
    if markers.get('sY14') == 'absent' and markers.get('ZFX/ZFY', markers.get('ZFX/Y')) == 'present':
        # Check if all Y markers are absent (suggesting XX male)
        y_markers_absent = all(
            markers.get(m) == 'absent' 
            for m in ['sY84', 'sY86', 'sY127', 'sY134', 'sY254', 'sY255']
            if m in markers
        )
        if y_markers_absent:
            return "46,XX_MALE_OR_COMPLETE_Y_CHROMOSOME_ABSENCE"
    
    # Step 2: Basic deletion analysis
    azfa_deleted = check_basic_deletion(markers, 'AZFa')
    azfb_deleted = check_basic_deletion(markers, 'AZFb')
    azfc_deleted = check_basic_deletion(markers, 'AZFc')
    
    # Check for methodological errors
    if 'sY254' in markers and 'sY255' in markers:
        if markers['sY254'] != markers['sY255']:
            return "METHODOLOGICAL_ERROR: sY254 and sY255 results discordant (both should have same status)"
    
    # Step 3: Classify based on deletion patterns
    
    # Check for AZFabc deletion (all regions deleted)
    if azfa_deleted and azfb_deleted and azfc_deleted:
        # With sY14 absent, suggests 46,XX male or Y absence
        if markers.get('sY14') == 'absent':
            return "AZFABC_DELETION_WITH_46,XX_MALE_OR_Y_CHROMOSOME_ABSENCE"
        return "COMPLETE_AZFABC_DELETION"
    
    # Check for AZFbc deletion
    if azfb_deleted and azfc_deleted and not azfa_deleted:
        is_azfbc, azfbc_msg = check_azfbc_deletion(markers)
        if is_azfbc:
            return azfbc_msg
    
    # Check individual deletions with extension analysis
    if azfa_deleted and not azfb_deleted and not azfc_deleted:
        ext_result = check_extension_markers(markers, 'AZFa')
        if ext_result['complete']:
            return "COMPLETE_AZFA_DELETION (TESE: virtually_impossible)"
        else:
            return f"PARTIAL_OR_ATYPICAL_AZFA_DELETION ({ext_result['subtype']})"
    
    if azfb_deleted and not azfa_deleted and not azfc_deleted:
        ext_result = check_extension_markers(markers, 'AZFb')
        if ext_result['subtype'] == 'P5/proximal_P1':
            return f"COMPLETE_AZFB_DELETION_P5/PROXIMAL_P1 (TESE: {ext_result['details']['TESE_prognosis']})"
        elif ext_result['subtype'] == 'partial':
            return f"PARTIAL_AZFB_DELETION (sY1192 present, TESE: {ext_result['details']['TESE_prognosis']})"
        else:
            return f"AZFB_DELETION_ATYPICAL_PATTERN"
    
    if azfc_deleted and not azfa_deleted and not azfb_deleted:
        ext_result = check_extension_markers(markers, 'AZFc')
        if ext_result['subtype'] == 'b2/b4':
            return f"COMPLETE_AZFC_DELETION_B2/B4 (TESE: {ext_result['details']['TESE_prognosis']})"
        elif ext_result['subtype'] == 'terminal':
            return "TERMINAL_AZFC_DELETION (karyotype analysis recommended)"
        else:
            return "AZFC_DELETION_UNDETERMINED_TYPE (sY160 not tested)"
    
    # Check for gr/gr partial deletion (can coexist with no complete deletions)
    if check_grgr_deletion(markers):
        if not (azfa_deleted or azfb_deleted or azfc_deleted):
            return "PARTIAL_AZFC_GR/GR_DELETION (population-specific risk factor)"
        else:
            # gr/gr with other deletions
            return f"MULTIPLE_DELETIONS_INCLUDING_GR/GR"
    
    # No deletions detected
    if not (azfa_deleted or azfb_deleted or azfc_deleted):
        # Verify expected present markers are actually present
        missing_expected = []
        for marker in EXPECTED_PRESENT_MARKERS:
            if marker in markers and markers[marker] != 'present':
                missing_expected.append(marker)
        
        if missing_expected:
            return f"INCONSISTENT_DATA: Expected present markers are absent: {', '.join(missing_expected)}"
        
        return "NO_DELETION_DETECTED"
    
    # Unclassified pattern
    return "UNCLASSIFIED_DELETION_PATTERN (manual review required)"


def generate_report(markers: Dict[str, str], classification: str) -> str:
    """Generate a detailed report with recommendations."""
    report = []
    report.append("=== Y-CHROMOSOMAL MICRODELETION ANALYSIS REPORT ===")
    report.append(f"\nCLASSIFICATION: {classification}")
    
    # Add marker summary
    report.append("\nMARKER SUMMARY:")
    report.append("Control markers:")
    for marker in CONTROL_MARKERS:
        if marker == 'ZFX/ZFY':
            # Handle both variants
            status = markers.get('ZFX/ZFY') or markers.get('ZFX/Y') or 'not_tested'
        else:
            status = markers.get(marker, 'not_tested')
        report.append(f"  {marker}: {status}")
    
    report.append("\nBasic AZF markers:")
    for region, region_markers in BASIC_MARKERS.items():
        report.append(f"  {region}:")
        for marker in region_markers:
            status = markers.get(marker, 'not_tested')
            report.append(f"    {marker}: {status}")
    
    # Add clinical recommendations based on classification
    report.append("\nCLINICAL RECOMMENDATIONS:")
    
    if "AZFA_DELETION" in classification and "COMPLETE" in classification:
        report.append("- TESE not recommended (virtually impossible to retrieve sperm)")
        report.append("- Consider sperm donation or adoption")
    elif "AZFB_DELETION" in classification and "P5/PROXIMAL_P1" in classification:
        report.append("- TESE not recommended (virtually impossible to retrieve sperm)")
        report.append("- Consider sperm donation or adoption")
    elif "AZFC_DELETION" in classification and "B2/B4" in classification:
        report.append("- TESE may be attempted (approximately 50% success rate)")
        report.append("- Genetic counseling mandatory (deletion will be transmitted to male offspring)")
        report.append("- Consider karyotype analysis to rule out 46,XY/45,X mosaicism")
    elif "TERMINAL" in classification:
        report.append("- Karyotype analysis strongly recommended")
        report.append("- Check for 46,XY/45,X mosaicism")
    elif "GR/GR" in classification:
        report.append("- Population-specific risk factor for impaired spermatogenesis")
        report.append("- Increased risk for testicular germ cell tumors")
        report.append("- Will be transmitted to male offspring")
    elif "46,XX_MALE" in classification:
        report.append("- Karyotype confirmation recommended")
        report.append("- TESE not possible")
        report.append("- Consider genetic counseling")
    elif "NO_DELETION" in classification:
        report.append("- No Y-chromosomal microdeletion detected")
        report.append("- Consider other causes of infertility")
    
    # Add warnings for missing markers
    missing_critical = []
    for marker in ALL_REQUIRED_MARKERS:
        if marker == 'ZFX/ZFY':
            # Check both variants
            if 'ZFX/ZFY' not in markers and 'ZFX/Y' not in markers:
                missing_critical.append(marker)
        elif marker not in markers:
            missing_critical.append(marker)
    
    if missing_critical:
        report.append("\nWARNING: The following markers were not tested:")
        for marker in missing_critical:
            report.append(f"  - {marker}")
        report.append("Complete testing is recommended for accurate diagnosis.")
    
    report.append("\n" + "=" * 50)
    
    return "\n".join(report)


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Classify Y-chromosomal microdeletions according to EAA/EMQN 2023 guidelines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Expected input format (TSV):
  marker_name<TAB>status
  
Where status is 'present' or 'absent' (case-insensitive).

Example:
  sY14	present
  sY84	absent
  sY86	absent
  ...

For complete analysis, all markers from Appendix B should be tested:
- Control: sY14, ZFX/ZFY
- Basic: sY84, sY86, sY127, sY134, sY254, sY255
- Extension: sY82, sY1064, sY1065/sY1182, sY88, sY105, sY121, sY1224, sY1192, sY153, sY160
- gr/gr: sY1291, sY1191
        """
    )
    
    parser.add_argument('marker_file', help='Path to TSV file containing marker results')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Generate detailed report with recommendations')
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate that all required markers are present')
    
    args = parser.parse_args()
    
    # Parse marker file
    marker_data = parse_marker_file(args.marker_file)
    if not marker_data:
        sys.exit(1)
    
    # Validate mode - just check if all markers are present
    if args.validate_only:
        missing = [m for m in ALL_REQUIRED_MARKERS if m not in marker_data]
        if missing:
            print(f"Missing {len(missing)} required markers:")
            for m in missing:
                print(f"  - {m}")
            sys.exit(1)
        else:
            print("All required markers are present.")
            sys.exit(0)
    
    # Classify
    classification = classify_y_chromosome_state(marker_data)
    
    if args.verbose:
        # Generate and print detailed report
        report = generate_report(marker_data, classification)
        print(report)
    else:
        # Just print classification
        print(classification)


if __name__ == "__main__":
    main()