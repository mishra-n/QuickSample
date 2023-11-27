from header import *
from helpers import *
from line_info import *
from select_absorbers import *
from selector import *
def plot_absorber_quick_look(fig, ax, spectrum, quick_look_table, identification_table, vrange=300, n=1):
    """
    Plot the absorber features in the given spectrum.

    Parameters:
    - fig (matplotlib.figure.Figure): The matplotlib figure for the plot.
    - ax (array of matplotlib.axes.Axes): An array of subplots for individual absorber plots.
    - spectrum (Spectrum1D): The input spectrum.
    - line_list (list): List of spectral lines to plot.
    - z (float): The redshift value.

    Returns:
    - None
    """

    line_list = np.unique(quick_look_table['line'])

    for i, line_name in enumerate(line_list):
        subtable = quick_look_table[quick_look_table['line'] == line_name]
        z = subtable['z'][0]

        # Get the wavelength of the spectral line
        line_wave = SEARCH_LINES[SEARCH_LINES['tempname']==line_name]['wave'][0]
        ion_name = SEARCH_LINES[SEARCH_LINES['tempname']==line_name]['name'][0]
        line_name = SEARCH_LINES[SEARCH_LINES['tempname']==line_name]['tempname'][0]
        line = SEARCH_LINES[SEARCH_LINES['tempname']==line_name]

        temp = 0
        if (line_wave * (1 + z) < 1100) or (line_wave * (1 + z) > 1790):
            temp += 1
            ax[i].set_visible(False)

        # Convert the full velocity axis to km/s
        full_velocity_axis = to_velocity(spectrum.spectral_axis, line_wave, z).to('km/s')

        # Averaging for better visualization
        length = spectrum.shape[0]
        remainder = length % n
        avg_flux = np.average(spectrum.flux[remainder::].reshape(-1, n), axis=1)
        avg_std = np.sqrt(np.sum((spectrum.uncertainty.array[remainder::].reshape(-1, n))**2, axis=1)/n**2)
        #avg_std = np.sqrt(avg_variance)
        avg_velocity_axis = np.average(full_velocity_axis[remainder::].reshape(-1, n), axis=1)

        # Plotting
        # if len(line_list) > 1:
        ax[i].text(0.98, 0.12, str(ion_name) + ' [' + str(np.round(line_wave, decimals=2)) + ' $\mathrm{\AA}$]', horizontalalignment='right', verticalalignment='bottom', transform=ax[i].transAxes, fontsize='x-large')
        ax[i].step(avg_velocity_axis.value, avg_flux, alpha=0.5, color='blue', where='mid')
        ax[i].fill_between(avg_velocity_axis.value,
                            (avg_flux-avg_std).value,
                            (avg_flux+avg_std).value,
                            color='blue', alpha=0.15, step='mid') 
        
        ax[i].vlines(0, 0, 1.5, color='black', linestyles='-.', alpha=0.5)
        ax[i].hlines(1, -vrange, vrange, color='slategrey', linestyles='dashed', alpha=0.9)
        # Grey out specific velocity range
        # if (gray_out is not None) and (gray_out[i] is not None):
        #     ax[i].fill_betweenx([0,1.5], [gray_out[i][0]]*2, [gray_out[i][1]]*2, color='lightgrey', alpha=0.5)

        # else:
        #     ax.text(0.98, 0.12, str(ion_name) + ' [' + str(np.round(line_wave, decimals=2)) + ' $\mathrm{\AA}$]', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, fontsize='x-large')
        #     ax.step(avg_velocity_axis.value, avg_flux, alpha=0.5, color='blue', where='mid')
        #     ax.fill_between(avg_velocity_axis.value,
        #                     (avg_flux-avg_std).value,
        #                     (avg_flux+avg_std).value,
        #                     color='blue', alpha=0.15, step='mid')
            
        #     ax.vlines(0, 0, 1.5, color='black', linestyles='-.', alpha=0.5)
        #     ax.hlines(1, -vrange, vrange, color='slategrey', linestyles='dashed', alpha=0.9)
        #     # if (gray_out is not None):
        #     #     ax[i].fill_betweenx([0,1.5], [gray_out[i][0]]*2, [gray_out[i][1]]*2, color='lightgrey', alpha=0.5)


        vmin = subtable['vmin']
        vmax = subtable['vmax']
        v = subtable['v']
        if subtable['v'] != 0 * u.km/u.s:
            ax[i].axvline(x=vmin.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
            ax[i].axvline(x=vmax.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
            ax[i].axvline(x=v.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
        if subtable['default_upper_limit'] == True:
            ax[i].axvline(x=-60, ls='--', ymin=0, ymax=2, color='orange')
            ax[i].axvline(x=60, ls='-.', ymin=0, ymax=2, color='orange')
            ax[i].axvline(x=0, ls='--', ymin=0, ymax=2, color='orange')

        inner_search_dv = [-250 * u.km/u.s, 350 * u.km/u.s]
        outer_search_dv = [-400 * u.km/u.s, 400 * u.km/u.s]

        sig_uk_inrange_peak_indices, sig_uk_inrange_peaks = find_significant_unknown_inrange_peaks(spectrum, line, z, identification_table, dv=inner_search_dv)
        inrange_trough_indices, inrange_troughs = find_inrange_troughs(spectrum, line, z, dv=outer_search_dv)
        
        
        ax[i].vlines(sig_uk_inrange_peaks.to('km/s'), ls='-.', ymin=0, ymax=2, color='black', alpha=0.2)
        ax[i].vlines(inrange_troughs.to('km/s'), ls='--', ymin=0, ymax=2, color='black', alpha=0.2)



        mask_out = mask_known_line_range(spectrum, line, z, identification_table, dv=outer_search_dv)
        new_known_mask = remove_current_line(line, z, identification_table)


        mask_in = np.invert(mask_out)

        ax[i].fill_between(full_velocity_axis.value, 0,2, where=mask_in, facecolor='red', alpha=.1)
        
        known_velocities = to_velocity(identification_table['wave_o'], line_wave, z)

        ax[i].vlines(known_velocities[new_known_mask].to('km/s'), ymin=0, ymax=2, color='red')


        ax[i].set_xlim(-vrange,vrange)
        ax[i].set_ylim(0, 1.5)
        if i == (len(line_list) - 1 - temp):
            ax[i].set_xlabel('Velocity (km/s)')

    fig.add_subplot(111, frameon=False)
    plt.grid(False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False, length=0, width=0, color='none', grid_alpha=0, which='both')

    plt.ylabel('Normalized Flux')


# Define the cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.30)

def plot_nearby_galaxies(fig, ax, central_galaxy, nearby_galaxies, dvs=[-500, 500] * u.km/u.s):

    #select nearby galaxies within the velocity range
    velocities = get_dvs(nearby_galaxies['REDSHIFT'], central_galaxy['REDSHIFT'])
    nearby_galaxies1 = nearby_galaxies[(velocities > dvs[0]) * (velocities < dvs[1])]
    nearby_galaxies2 = nearby_galaxies[(velocities > 2*dvs[0]) * (velocities < 2*dvs[1])]
    nearby_galaxies3 = nearby_galaxies[(velocities > 3*dvs[0]) * (velocities < 3*dvs[1])]
    # Define the zero point coordinates
    zero_point = SkyCoord(ra=central_galaxy['RA']*u.deg, dec=central_galaxy['DEC']*u.deg)
    
    # # Convert RA and DEC of nearby galaxies to distance coordinates
    # coords = SkyCoord(ra=nearby_galaxies['RA']*u.deg, dec=nearby_galaxies['DEC']*u.deg)
    # distances = coords.separation(zero_point)
    # print(distances)
    # Convert distances to angular diameter distance
    angular_diameter_distance = cosmo.angular_diameter_distance(z=central_galaxy['REDSHIFT']) #* distances.to('radian').value).to(u.Mpc)
    
    # Calculate dx and dy in angular diameter distance

    dx1 = (angular_diameter_distance * (nearby_galaxies1['RA']*u.deg - zero_point.ra).to('radian').value).to('Mpc') * np.cos(zero_point.dec.to('radian').value)
    dy1 = (angular_diameter_distance * (nearby_galaxies1['DEC']*u.deg - zero_point.dec).to('radian').value).to('Mpc')

    dx2 = (angular_diameter_distance * (nearby_galaxies2['RA']*u.deg - zero_point.ra).to('radian').value).to('Mpc') * np.cos(zero_point.dec.to('radian').value)
    dy2 = (angular_diameter_distance * (nearby_galaxies2['DEC']*u.deg - zero_point.dec).to('radian').value).to('Mpc')
    
    dx3 = (angular_diameter_distance * (nearby_galaxies3['RA']*u.deg - zero_point.ra).to('radian').value).to('Mpc') * np.cos(zero_point.dec.to('radian').value)
    dy3 = (angular_diameter_distance * (nearby_galaxies3['DEC']*u.deg - zero_point.dec).to('radian').value).to('Mpc')

    dx_qso = (angular_diameter_distance * (nearby_galaxies['RA_QSO']*u.deg - zero_point.ra).to('radian').value).to('Mpc') * np.cos(zero_point.dec.to('radian').value)
    dy_qso = (angular_diameter_distance * (nearby_galaxies['DEC_QSO']*u.deg - zero_point.dec).to('radian').value).to('Mpc')
    # Plot the nearby galaxies in distance space
    ax.scatter(dx1, dy1, color='blue', label='within 500 km/s', s=9)
    ax.scatter(dx2, dy2, color='blue', label='within 1000 km/s', alpha=0.5, s=9)
    ax.scatter(dx3, dy3, color='blue', label='within 1500 km/s', alpha=0.3, s=9)

    # ax.scatter(dx1, dy1, color='blue', label='within 500 km/s')
    # ax.scatter(dx2, dy2, color='blue', label='within 1000 km/s', alpha=0.5)
    # ax.scatter(dx3, dy3, color='blue', label='within 1500 km/s', alpha=0.3)

    # Draw dotted line circles at the virial radius of each galaxy
    for i in range(len(nearby_galaxies1)):
        #print(dx1[i].to('Mpc').value, dy1[i].to('Mpc').value, (nearby_galaxies1['RVIR'][i] * u.kpc).to('Mpc'))
        ax.add_patch(plt.Circle((dx1[i].to('Mpc').value, dy1[i].to('Mpc').value), (nearby_galaxies1['RVIR'][i] * u.kpc).to('Mpc').value, fill=False, linestyle='dotted', edgecolor='blue'))

    for i in range(len(nearby_galaxies2)):
        ax.add_patch(plt.Circle((dx2[i].to('Mpc').value, dy2[i].to('Mpc').value), (nearby_galaxies2['RVIR'][i] * u.kpc).to('Mpc').value, fill=False, linestyle='dotted', edgecolor='blue', alpha=0.5))

    for i in range(len(nearby_galaxies3)):
        ax.add_patch(plt.Circle((dx3[i].to('Mpc').value, dy3[i].to('Mpc').value), (nearby_galaxies3['RVIR'][i] * u.kpc).to('Mpc').value, fill=False, linestyle='dotted', edgecolor='blue', alpha=0.3))

    #for i in range(len(nearby_galaxies)):
    ax.add_patch(plt.Circle((0,0), (central_galaxy['D'] * u.kpc).to('Mpc').value, fill=False, linestyle='dotted', edgecolor='orange', alpha=0.3))
    ax.scatter(dx_qso, dy_qso, color='orange', label='QSO', s=9)
    ax.scatter(0, 0, color='red', label='Central Galaxy', s=9)
    ax.add_patch(plt.Circle((0, 0), (central_galaxy['RVIR'] * u.kpc).to('Mpc').value, fill=False, linestyle='dotted', edgecolor='red'))

    ax.set_xlabel('x [Mpc]', fontsize='large')
    ax.set_ylabel('y [Mpc]', fontsize='large')
    QSO = central_galaxy['QSO'][0]
    lstar = np.round(central_galaxy['LSTAR'], 5)[0]
    redshift = np.round(central_galaxy['REDSHIFT'], 4)[0]
    ax.set_title('QSO:' + str(QSO) + ', $L_* = $' + str(lstar) + ', $z = $' + str(redshift), fontsize='x-large')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.legend(fontsize='large')
    #plt.show()
