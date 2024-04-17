from header import *  # Importing necessary dependencies
from line_info import *  # Importing additional dependencies from other modules
from helpers import *  # Importing helper functions from a different module
from select_absorbers import *
deltav = 60

def find_open_area(spectrum, known_lines, line, z, dv=60*u.km/u.second):
    wavelength = spectrum.wavelength
    mask = full_mask(spectrum, known_lines)
    uncertainty = np.ma.array(spectrum.uncertainty.array, mask=mask)

    line_wave = line['wave']
    velocity_axis = to_velocity(wavelength, line_wave, z)
    v_resolution = velocity_axis[1::] - velocity_axis[0:-1]
    v_resolution = np.append(v_resolution, v_resolution[-1])
    #set v_resolution to 0 where the uncertainty is masked
    v_resolution[uncertainty.mask] = 0
    # print(line_wave, z)
    # print(velocity_axis)
    # print(np.abs(velocity_axis))
    zero_index = np.argmin(np.abs(velocity_axis))
    # Find the indices of the unmasked uncertainty values
    unmasked_indices = np.where(~uncertainty.mask)[0]

    
    # Calculate the velocity difference between each unmasked point and the line wavelength
    summed_velocities_from_zeros_pos = np.cumsum(v_resolution[zero_index::])
    #print(summed_velocities_from_zeros_pos.to('km/s'))
    # print(zero_index, np.round(z, decimals=3), len(velocity_axis), summed_velocities_from_zeros_pos)
    num_points_to_dv_pos = np.where(summed_velocities_from_zeros_pos > dv)[0][0]
    summed_velocities_from_zeros_neg = np.cumsum(v_resolution[zero_index::-1])
    num_points_to_dv_neg = np.where(summed_velocities_from_zeros_neg > dv)[0][0]
    pos_velocities = velocity_axis[zero_index:zero_index+num_points_to_dv_pos]
    neg_velocities = velocity_axis[zero_index-num_points_to_dv_neg:zero_index][::-1]
    pos_vrange = [min(pos_velocities), max(pos_velocities)]
    neg_vrange = [min(neg_velocities), max(neg_velocities)]

    return pos_vrange, neg_vrange


def calc_Wr(spectrum, line, z, v=0*u.km/u.second, dv=[-deltav*u.km/u.second, deltav*u.km/u.second], mask=1.0):
    line = line['wave'][0]
    velocity_axis = to_velocity(spectrum.spectral_axis, line, z)

    v_cut = np.logical_and((velocity_axis > (dv[0])), velocity_axis < (dv[1]))

    F_cut = ((1 - spectrum.flux) * mask) [v_cut]
    F_err_cut = (spectrum.uncertainty.array * mask) [v_cut]

    rest_spectral_axis = spectrum.spectral_axis / (1 + z)
    wave_cut = rest_spectral_axis[v_cut]
    index = find_nearest_index(spectrum.spectral_axis.value, line*(1 + z))
    dW = (rest_spectral_axis[index] - rest_spectral_axis[index - 1])

    Wr = np.trapz(F_cut, wave_cut)
    #Wr = np.sum(F_cut) * dW
    WrErr = np.sqrt(np.trapz(F_err_cut**2, wave_cut) * dW)
    return Wr, WrErr

def calc_AOD(spectrum, line, z, v=0*u.km/u.second, dv=[-deltav*u.km/u.second, deltav*u.km/u.second], mask=1.0):

    
    f = line['f'][0]
    line = line['wave'][0]
    velocity_axis = to_velocity(spectrum.spectral_axis, line, z)
    v_cut = np.logical_and((velocity_axis > (dv[0])), velocity_axis < (dv[1]))

    F_cut = (1 - ((1 - spectrum.flux) * mask)) [v_cut]
    F_err_cut = (spectrum.uncertainty.array * mask) [v_cut]

    tau = np.log(1 / F_cut)
    tau_err = F_err_cut / F_cut
    
    tau[tau == np.inf] = 0
    tau[np.isnan(tau)] = 0
    tau_err[tau_err == np.inf] = 0

    rest_spectral_axis = spectrum.spectral_axis / (1 + z)
    wave_cut = spectrum.spectral_axis[v_cut]
    index = find_nearest_index(spectrum.spectral_axis.value, line*(1 + z))
    dW = (rest_spectral_axis[index] - rest_spectral_axis[index - 1])

    AOD = np.trapz(tau, wave_cut)
    AODErr = np.sqrt(np.trapz(tau_err**2, wave_cut) * dW)
    return AOD, AODErr

def AOD_to_N(AOD, AODErr, line):
    N = m_e * c**2 / np.pi / e.gauss**2 / line['f'][0] / (line['wave'][0] * u.Angstrom)**2 * AOD
    NErr = m_e * c**2 / np.pi / e.gauss**2 / line['f'][0] / (line['wave'][0] * u.Angstrom)**2 * AODErr

    return N.to(u.cm**(-2)), NErr.to(u.cm**(-2))

def Wr_to_N(Wr, WrErr, line):
    N = m_e * c**2 / np.pi / e.gauss**2 / line['f'][0] / (line['wave'][0] * u.Angstrom)**2 * Wr
    NErr = m_e * c**2 / np.pi / e.gauss**2 / line['f'][0] / (line['wave'][0] * u.Angstrom)**2 * WrErr

    return N.to(u.cm**(-2)), NErr.to(u.cm**(-2))





