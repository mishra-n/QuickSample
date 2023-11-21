from helpers import *
from header import *
from spectral import *
from selector import *
from select_absorbers import *
from upper_limits import *


def genIdentificationTable(path, field_name):
    """
    Generate an identification table from a FITS file.

    Parameters:
    path (str): The path to the FITS file.
    field_name (str): The name of the field.

    Returns:
    None
    """
    known_lines = QTable.read(path+field_name+'_final_abscal_identifications.fits')
    known_lines['wave_r'] = known_lines['wave_r'] * u.Angstrom
    known_lines['wave_o'] = known_lines['wave_o'] * u.Angstrom
    known_lines['Wr'] = known_lines['Wr'] * u.Angstrom
    known_lines['WrErr'] = known_lines['WrErr'] * u.Angstrom
    known_lines['b'] = known_lines['b'] * u.km / u.second
    known_lines['bErr'] = known_lines['bErr'] * u.km / u.second

    return known_lines


def genQuickLookTableColumns():
    """
    Generate the columns for the quick look table.

    Returns:
    QTable: The quick look table columns.
    """
    new_table = QTable(names=['name',
                              'z',
                              'line',
                              'ion',
                              'out_of_range',
                              'upper_limit',
                              'default_upper_limit',
                              'default_SNR',
                              'SNR',
                              'default_Wr',
                              'default_WrErr',
                              'default_AOD',
                              'default_AODErr',
                              'v',
                              'vmin',
                              'vmax',
                              'Wr',
                              'WrErr',
                              'AOD',
                              'AODErr',
                              'N',
                              'NErr',
                              'UL_N',
                              'default_UL_N'],
                      dtype=['str','f4', 'S2', 'S2', 
                             'bool', 'bool', 'bool',
                             'f4', 'f4',
                             'f4', 'f4', 'f4', 'f4', 
                             'f4', 'f4', 'f4', 
                             'f4', 'f4', 'f4', 'f4',
                             'f4', 'f4', 'f4', 'f4'],
                      units={'default_Wr':'Angstrom',
                             'default_WrErr':'Angstrom',
                             'default_AOD':'Angstrom',
                             'default_AODErr':'Angstrom',
                             'v':'km/s',
                             'vmin':'km/s',
                             'vmax':'km/s',
                             'Wr':'Angstrom',
                             'WrErr':'Angstrom',
                             'AOD':'Angstrom',
                             'AODErr':'Angstrom',
                             'N': '1/cm**2',
                             'NErr': '1/cm**2',
                             'UL_N': '1/cm**2',
                             'default_UL_N':'1/cm**2'})
    return new_table


def addIonRows(z, name, selected_lines=np.array(['HI','OVI', 'CIII']), new_table=None):
    """
    Add rows for specific ions to the quick look table.

    Parameters:
    z (float): The redshift value.
    name (str): The name of the galaxy.
    selected_lines (np.array): The selected lines to add.
    new_table (QTable): The quick look table.

    Returns:
    QTable: The updated quick look table.
    """
    if np.array_equal(new_table, None):
        new_table = genQuickLookTableColumns()
    for line in SEARCH_LINES:
        if (line['name'] == selected_lines).sum() > 0:
            new_table.add_row(vals={'name':name, 'z':z, 'line': line['tempname'], 'ion': line['name']})
    return new_table


def genQuickLookTable(galaxy_sample_table, identification_table, spectrum, inner_search_dv, outer_search_dv, upper_limit_dv):
    """
    Generate the quick look table.

    Parameters:
    galaxy_sample_table (QTable): The galaxy sample table.
    identification_table (QTable): The identification table.
    spectrum (Spectrum): The spectrum data.
    inner_search_dv (float): The inner search velocity range.
    outer_search_dv (float): The outer search velocity range.
    upper_limit_dv (float): The upper limit velocity range.

    Returns:
    QTable: The quick look table.
    """
    new_table = genQuickLookTableColumns()
    for galaxy in galaxy_sample_table:
        z = galaxy['REDSHIFT']
        name = galaxy['SHORTNAME']
        new_table = addIonRows(z, name, new_table=new_table)
    
    for row in new_table:
        line = SEARCH_LINES[SEARCH_LINES['tempname']==row['line']]
        if line['wave'][0]*(1+z) < 1100 or line['wave'][0]*(1+z) > 1790:
            row['out_of_range'] = True
            continue
        
        inrange_trough_indices, inrange_troughs = find_inrange_troughs(spectrum, line, z, dv=outer_search_dv)
        sig_uk_inrange_peak_indices, sig_uk_inrange_peaks = find_significant_unknown_inrange_peaks(spectrum, line, z, identification_table, dv=outer_search_dv)

        mask = mask_known_line_range(spectrum, line, z, identification_table, dv=outer_search_dv)
        if len(sig_uk_inrange_peaks) !=0:
            print(line['wave'][0])
            left_peak_bound, right_peak_bound, center_peak, left_peak_bound_index, right_peak_bound_index, center_peak_index = find_uk_peak_range(spectrum, line, z, identification_table, dv=outer_search_dv)
            row['vmin'] = left_peak_bound
            row['vmax'] = right_peak_bound
            row['v'] = center_peak
            AOD, AODErr = calc_AOD(spectrum, line, z, v=center_peak, dv=[left_peak_bound, right_peak_bound])
            Wr, WrErr = calc_Wr(spectrum, line, z, v=center_peak, dv=[left_peak_bound, right_peak_bound])
            row['Wr'] = Wr
            row['WrErr'] = WrErr
            row['SNR'] = Wr/WrErr
            row['AOD'] = AOD
            row['AODErr'] = AODErr

        free = default_range_free(spectrum, line, z, identification_table, dv=outer_search_dv, upper_limit_dv=upper_limit_dv)
        if free:
            print('free')
            AOD, AODErr = calc_AOD(spectrum, line, z, v=0, dv=[-60 * u.km/u.second, 60 * u.km/u.second])
            Wr, WrErr = calc_Wr(spectrum, line, z, v=0, dv=[-60 * u.km/u.second, 60 * u.km/u.second])
            row['default_Wr'] = Wr
            row['default_WrErr'] = WrErr
            row['default_SNR'] = Wr/WrErr
            row['default_AOD'] = AOD
            row['default_AODErr'] = AODErr
            row['default_upper_limit'] = True
            Wr, WrErr = calc_Wr(spectrum, line, z, v=0, dv=[-60 * u.km/u.second, 60 * u.km/u.second], mask=mask)
            AOD, AODErr = calc_AOD(spectrum, line, z, v=0, dv=[-60 * u.km/u.second, 60 * u.km/u.second], mask=mask)
            N, NErr = AOD_to_N(AOD, AODErr, line)
            row['default_UL_N'] = 3 * NErr
            
        max_SNR = new_table[new_table['name'] ==  name]['SNR'].max()
        selected_row = new_table[(new_table['SNR'] == max_SNR)]
        vmin = selected_row['vmin'][0]
        vmax = selected_row['vmax'][0]
        v = selected_row['v'][0]

        for row in new_table:
            if row['name'] != name:
                continue
            if max_SNR < 3:
                continue
            line = SEARCH_LINES[SEARCH_LINES['tempname'] == row['line']]
            if line['wave']*(1+z) < 1115 or line['wave']*(1+z) > 1790:
                row['out_of_range'] = True
                continue
            
            inrange_trough_indices, inrange_troughs = find_inrange_troughs(spectrum, line, z, dv=outer_search_dv)
            sig_uk_inrange_peak_indices, sig_uk_inrange_peaks = find_significant_unknown_inrange_peaks(spectrum, line, z, identification_table, dv=outer_search_dv)

            mask = mask_known_line_range(spectrum, line, z, identification_table, dv=outer_search_dv)
            
            row['vmin'] = vmin
            row['vmax'] = vmax
            row['v'] = v
            AOD, AODErr = calc_AOD(spectrum, line, z, v=v, dv=[vmin, vmax], mask=mask)
            Wr, WrErr = calc_Wr(spectrum, line, z, v=v, dv=[vmin, vmax], mask=mask)
            row['Wr'] = Wr
            row['WrErr'] = WrErr

            row['Wr'] = Wr
            row['SNR'] = Wr/WrErr
            row['AOD'] = AOD
            row['AODErr'] = AODErr

            N, NErr = AOD_to_N(AOD, AODErr, line)
            row['N'] = N
            row['NErr'] = NErr
            row['UL_N'] = 3 * NErr
    return new_table
