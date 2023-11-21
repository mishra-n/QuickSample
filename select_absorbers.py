from header import *
from helpers import *

SIGMA = 10  # Standard deviation for Gaussian filter
FLUX_CUT = 0.8  # Flux threshold for significant peaks
KNOWN_DIFF = 100 * u.km/u.second  # Known velocity difference
DV = [-400 * u.km/u.second, 400 * u.km/u.second]  # Velocity range
DV2 = [-250 * u.km/u.second, 250 * u.km/u.second]  # Velocity range
default_DV = [-60 * u.km/u.second, 60 * u.km/u.second]

filtered_flux = None

def find_critical_points(spectrum, line, z):
    """
    Finds the critical points (zero-crossings) in the spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        
    Returns:
        numpy.ndarray: Indices of the critical points.
        numpy.ndarray: Values of the critical points.
    """
    line_wave = line['wave'][0]
    global filtered_flux
    if filtered_flux is None:
        filtered_flux = gaussian_filter1d(spectrum.flux, sigma=SIGMA)
    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)
    gradient = np.gradient(filtered_flux)
    crit_indices, crit_points = find_zero_crossings(velocity_axis, gradient)
    return crit_indices, crit_points

def find_peaks(spectrum, line, z):
    """
    Finds the peaks in the spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        
    Returns:
        numpy.ndarray: Indices of the peaks.
        numpy.ndarray: Values of the peaks.
    """
    line_wave = line['wave'][0]
    global filtered_flux
    if filtered_flux is None:
        filtered_flux = gaussian_filter1d(spectrum.flux, sigma=SIGMA)

    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)
    crit_indices, crit_points = find_critical_points(spectrum, line, z)
    second_grad = np.gradient(np.gradient(filtered_flux))
    is_peak_bool_arr = np.sign(second_grad[crit_indices]) == 1
    peaks = crit_points[is_peak_bool_arr]
    peak_indices = crit_indices[is_peak_bool_arr]
    return peak_indices, peaks

def find_troughs(spectrum, line, z):
    """
    Finds the troughs in the spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        
    Returns:
        numpy.ndarray: Indices of the troughs.
        numpy.ndarray: Values of the troughs.
    """
    line_wave = line['wave'][0]
    global filtered_flux
    if filtered_flux is None:
        filtered_flux = gaussian_filter1d(spectrum.flux, sigma=SIGMA)

    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)
    crit_indices, crit_points = find_critical_points(spectrum, line, z)
    second_grad = np.gradient(np.gradient(filtered_flux))
    is_trough_bool_arr = np.sign(second_grad[crit_indices]) == -1
    troughs = crit_points[is_trough_bool_arr]
    trough_indices = crit_indices[is_trough_bool_arr]
    return trough_indices, troughs    

def find_significant_peaks(spectrum, line, z):
    """
    Finds the significant peaks in the spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        
    Returns:
        numpy.ndarray: Indices of the significant peaks.
        numpy.ndarray: Values of the significant peaks.
    """
    line_wave = line['wave'][0]
    global filtered_flux
    if filtered_flux is None:
        filtered_flux = gaussian_filter1d(spectrum.flux, sigma=SIGMA)
    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)
    peak_indices, peaks = find_peaks(spectrum, line, z)
    peak_fluxes = filtered_flux[peak_indices]
    significant_bool_arr = peak_fluxes < FLUX_CUT
    sig_peaks = peaks[significant_bool_arr]
    sig_peak_indices = peak_indices[significant_bool_arr]

    return sig_peak_indices, sig_peaks


def remove_values_within_range(A, B, n):
    """
    Removes values from A that are within a certain range of values in B.
    
    Args:
        A (numpy.ndarray): Array containing values to filter.
        B (numpy.ndarray): Array containing reference values.
        n (float): Range threshold.
        
    Returns:
        tuple: Filtered array A and corresponding indices.
    """
    # Calculate the absolute difference between each element of A and each element of B
    diff = np.abs(A[:, None] - B)

    # Check if each difference is less than or equal to n
    mask = np.any(diff <= n, axis=1)

    # Invert the mask to get the indices of the elements of A that are not within n of any element of B
    selected_indices = np.where(~mask)[0]

    # Use the selected indices to filter A
    filtered_A = A[selected_indices]

    return filtered_A, selected_indices


def remove_current_line(line, z, known_lines):
    """
    Creates a mask to remove the current line from a set of known lines.
    
    Args:
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        known_lines (numpy.ndarray): Array of known lines.
        
    Returns:
        numpy.ndarray: Mask indicating the lines to be removed.
    """
    mask = np.absolute(known_lines['wave_r'] - line['wave'][0] * u.Angstrom) > 0.01 * u.Angstrom
    mask2 = np.absolute(known_lines['z'] - z) > 0.03
    return np.array(mask + mask2, dtype=bool)


def find_significant_unknown_peaks(spectrum, line, z, known_lines):
    """
    Finds significant unknown peaks in a spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        known_lines (numpy.ndarray): Array of known lines.
        
    Returns:
        tuple: Indices and values of significant unknown peaks.
    """
    line_wave = line['wave'][0]
    sig_peak_indices, sig_peaks = find_significant_peaks(spectrum, line, z)
    known_velocities = to_velocity(known_lines['wave_o'], line_wave, z)

    new_known_mask = remove_current_line(line, z, known_lines)

    sig_uk_peaks, indices_temp = remove_values_within_range(sig_peaks, known_velocities[new_known_mask], n=50 * u.km/u.second)
    sig_uk_peak_indices = sig_peak_indices[indices_temp]

    return sig_uk_peak_indices, sig_uk_peaks


def find_significant_unknown_inrange_peaks(spectrum, line, z, known_lines, dv):
    """
    Finds significant unknown peaks within a specified velocity range.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        known_lines (numpy.ndarray): Array of known lines.
        dv (tuple): Velocity range.
        
    Returns:
        tuple: Indices and values of significant unknown peaks within the range.
    """
    sig_uk_peak_indices, sig_uk_peaks = find_significant_unknown_peaks(spectrum, line, z, known_lines)
    mask = (sig_uk_peaks > dv[0]) * (sig_uk_peaks < dv[1])

    sig_uk_inrange_indices = sig_uk_peak_indices[mask]
    sig_uk_inrange_peaks = sig_uk_peaks[mask]

    return sig_uk_inrange_indices, sig_uk_inrange_peaks


def find_inrange_troughs(spectrum, line, z, dv):
    """
    Finds troughs within a specified velocity range.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        dv (tuple): Velocity range.
        
    Returns:
        tuple: Indices and values of troughs within the range.
    """
    trough_indices, troughs = find_troughs(spectrum, line, z)
    mask = (troughs > 2*dv[0]) * (troughs < 2*dv[1])
    inrange_trough_indices = trough_indices[mask]
    inrange_troughs = troughs[mask]  

    return  inrange_trough_indices, inrange_troughs


def find_right(arr, n):
    """
    Finds the value in an array that is closest to the given value and greater than it.
    
    Args:
        arr (numpy.ndarray): Input array.
        n (float): Reference value.
        
    Returns:
        float: Closest value in the array greater than the reference value.
    """
    diff = arr - n  # Calculate the difference between the array and n
    mask = diff <= 0  # Create a mask for positive differences (values greater than n)
    masked_diff = np.ma.masked_array(diff, mask)  # Apply the mask to the differences
    closest_index = np.argmin(masked_diff)  # Find the index of the minimum positive difference
    
    return arr[closest_index]


def find_left(arr, n):
    """
    Finds the value in an array that is closest to the given value and less than it.
    
    Args:
        arr (numpy.ndarray): Input array.
        n (float): Reference value.
        
    Returns:
        float: Closest value in the array less than the reference value.
    """
    diff = arr - n  # Calculate the difference between the array and n
    mask = diff >= 0  # Create a mask for positive differences (values greater than n)
    masked_diff = np.ma.masked_array(diff, mask)  # Apply the mask to the differences
    closest_index = np.argmax(masked_diff)  # Find the index of the minimum positive difference
    
    return arr[closest_index]


def find_uk_peak_range(spectrum, line, z, known_lines, dv):
    """
    Finds the range of significant unknown peaks within a spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        known_lines (numpy.ndarray): Array of known lines.
        dv (tuple): Velocity range.
        
    Returns:
        tuple: Left peak bound, right peak bound, center peak, left peak bound index, 
               right peak bound index, center peak index.
    """
    line_wave = line['wave'][0]
    sig_uk_inrange_indices, sig_uk_inrange_peaks = find_significant_unknown_inrange_peaks(spectrum, line, z, known_lines, dv)
    inrange_trough_indices, inrange_troughs = find_inrange_troughs(spectrum, line, z, dv)

    min_peak_index = sig_uk_inrange_indices.min()
    min_peak = sig_uk_inrange_peaks.min()
    max_peak_index = sig_uk_inrange_indices.min()
    max_peak = sig_uk_inrange_peaks.max()

    left_peak_bound = find_left(inrange_troughs, min_peak)
    right_peak_bound = find_right(inrange_troughs, max_peak)
    left_peak_bound_index = find_left(inrange_trough_indices, min_peak_index)
    right_peak_bound_index = find_right(inrange_trough_indices, max_peak_index)

    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)

    center_peak_index = spectrum.flux[left_peak_bound_index:right_peak_bound_index].argmin() + left_peak_bound_index
    center_peak = velocity_axis[center_peak_index]

    return left_peak_bound, right_peak_bound, center_peak, left_peak_bound_index, right_peak_bound_index, center_peak_index
   

def mask_known_line_range(spectrum, line, z, known_lines, dv):
    """
    Masks the known line range within a spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        known_lines (numpy.ndarray): Array of known lines.
        dv (tuple): Velocity range.
        
    Returns:
        numpy.ndarray: Mask indicating the known line range.
    """
    line_wave = line['wave'][0]

    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)
    new_known_mask = remove_current_line(line, z, known_lines)
    known_velocities = to_velocity(known_lines['wave_o'], line_wave, z)
    new_known_velocities = known_velocities[new_known_mask]

    velocity_check_mask = (new_known_velocities > 3*dv[0]) * (new_known_velocities < 3*dv[1])

    new_known_velocities = new_known_velocities[velocity_check_mask]
    inrange_trough_indices, inrange_troughs = find_inrange_troughs(spectrum, line, z, dv=[3*DV[0], 3*DV[1]])

    
    mask = np.ones(shape=velocity_axis.shape)
    for known_velocity in new_known_velocities:
        left_peak_bound = find_left(inrange_troughs, known_velocity)
        right_peak_bound = find_right(inrange_troughs, known_velocity)
        mask *= (velocity_axis > right_peak_bound) + (velocity_axis < left_peak_bound)
    #mask *= (velocity_axis > dv[0]) * (velocity_axis < dv[1])

    return np.array(mask, dtype='bool')


def mask_uk_peak_range(spectrum, line, z, known_lines, dv):
    """
    Masks the range of significant unknown peaks within a spectrum.
    
    Args:
        spectrum (numpy.ndarray): Input spectrum.
        line (numpy.ndarray): Current line information.
        z (float): Redshift value.
        known_lines (numpy.ndarray): Array of known lines.
        dv (tuple): Velocity range.
        
    Returns:
        numpy.ndarray: Mask indicating the range of significant unknown peaks.
    """
    left_peak_bound, right_peak_bound, center_peak, left_peak_bound_index, right_peak_bound_index, center_peak_index = find_uk_peak_range(spectrum, line, z, known_lines, dv)

    line_wave = line['wave'][0]

    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)

    velocity_check_mask = (velocity_axis > left_peak_bound) * (velocity_axis < right_peak_bound)
    return np.array(velocity_check_mask, dtype='bool')

def default_range_free(spectrum, line, z, known_lines, dv, upper_limit_dv):
    line_wave = line['wave'][0]

    velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z)
    new_known_mask = remove_current_line(line, z, known_lines)
    known_velocities = to_velocity(known_lines['wave_o'], line_wave, z)


    sig_uk_inrange_peak_indices, sig_uk_inrange_peaks = find_significant_unknown_inrange_peaks(spectrum, line, z, known_lines, dv)
    if len(sig_uk_inrange_peaks) == 0:
        mask1 = np.zeros(shape=velocity_axis.shape[0])
    else:
        mask1 = mask_uk_peak_range(spectrum, line, z, known_lines, dv)

    mask2 = mask_known_line_range(spectrum, line, z, known_lines, dv)
    bound_mask = (velocity_axis > upper_limit_dv[0]) * (velocity_axis < upper_limit_dv[1])

    total_mask = (mask1 + np.invert(mask2))[bound_mask]
    total_mask = np.array(total_mask, dtype=bool)

    if total_mask.sum() == 0:
        return True
    else:
        return False
