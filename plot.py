from header import *
from helpers import *
from line_info import *



def plot_velocity_dispersion(ax, spectrum, line_object, z, label=False, xlabels=True, ylabels=True, calculations=None):
    line = line_object['wave']
    line_name = line_object['tempname']
    n=2
    avg_spec = np.average(spectrum.spectral_axis[:-1].reshape(-1, n), axis=1)
    avg_flux = np.average(spectrum.flux[:-1].reshape(-1, n), axis=1)
    avg_variance = np.average((spectrum.uncertainty.array**2)[:-1].reshape(-1, n), axis=1)
    avg_std = np.sqrt(avg_variance)

    line_velocity = to_velocity(avg_spec, line, z).to('km/s')
    
    x_lambda = np.linspace(line*(1+z) - 3*u.Angstrom, line*(1+z) + 3*u.Angstrom, 1000)
    line_velocity_unbinned = to_velocity(x_lambda, line, z).to('km/s')
    
    ax.step(line_velocity, avg_flux, 'blue', alpha=0.8)
    #ax.plot(line_velocity_unbinned, fit(x_lambda.value), 'red')
    ax.set_xlim(-400,400)
    ax.set_ylim(0,1.4)
    if xlabels: ax.set_xlabel("Velocity ({})".format(line_velocity.unit))
    ax.fill_between(line_velocity.value,
                     avg_flux-avg_std,
                     avg_flux+avg_std,
                     color='blue', alpha=0.15, step='pre')
    #kernel = Box1DKernel(width=5)
    #smoothed = convolve(spectrum['normalized_flux'], kernel)
    ax.axvline(x=0, ls='--', ymin=0, ymax=2, color='black', alpha=0.6)
    ax.axhline(y=1, ls='--', xmin=-500, xmax=500, color='black', alpha=0.6)

    #if calculations != None:
    if isinstance(calculations, QTable):
        vmin = calculations[calculations['line'] == line_name]['vmin']
        vmax = calculations[calculations['line'] == line_name]['vmax']
        v = calculations[calculations['line'] == line_name]['v']
        if calculations[calculations['line'] == line_name]['interloper_no_measurement'] == False:
            ax.axvline(x=vmin.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
            ax.axvline(x=vmax.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
            ax.axvline(x=v.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
        if calculations[calculations['line'] == line_name]['default_upper_limit'] == True:
            ax.axvline(x=-60, ls='--', ymin=0, ymax=2, color='red')
            ax.axvline(x=60, ls='--', ymin=0, ymax=2, color='red')
            ax.axvline(x=0, ls='--', ymin=0, ymax=2, color='red')

    if xlabels==False:
        ax.set_xticklabels([])
    if ylabels==False:
        ax.set_yticklabels([])
    ax.set_title(line_name + ',obs:' + str(np.round(line.value*(1+z), decimals=1) ) + ',tru:' + str(np.round(line.value, decimals=1) ))


def plot_absorber(fig, ax, spectrum, line_list, z, vrange=300, n=1, gray_out=None, calculations=None):
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
    nrows = len(line_list)

    for i, line in enumerate(line_list):
        # Get the wavelength of the spectral line
        line_wave = SEARCH_LINES[SEARCH_LINES['tempname']==line]['wave'][0]
        ion_name = SEARCH_LINES[SEARCH_LINES['tempname']==line]['name'][0]
        line_name = SEARCH_LINES[SEARCH_LINES['tempname']==line]['tempname'][0]

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
        if len(line_list) > 1:
            ax[i].text(0.98, 0.12, str(ion_name) + ' [' + str(np.round(line_wave, decimals=2)) + ' $\mathrm{\AA}$]', horizontalalignment='right', verticalalignment='bottom', transform=ax[i].transAxes, fontsize='x-large')
            ax[i].step(avg_velocity_axis.value, avg_flux, alpha=0.5, color='blue', where='mid')
            ax[i].fill_between(avg_velocity_axis.value,
                                (avg_flux-avg_std).value,
                                (avg_flux+avg_std).value,
                                color='blue', alpha=0.15, step='mid') 
            
            ax[i].vlines(0, 0, 1.5, color='black', linestyles='-.', alpha=0.5)
            ax[i].hlines(1, -vrange, vrange, color='slategrey', linestyles='dashed', alpha=0.9)
            # Grey out specific velocity range
            if (gray_out is not None) and (gray_out[i] is not None):
                ax[i].fill_betweenx([0,1.5], [gray_out[i][0]]*2, [gray_out[i][1]]*2, color='lightgrey', alpha=0.5)

        else:
            ax.text(0.98, 0.12, str(ion_name) + ' [' + str(np.round(line_wave, decimals=2)) + ' $\mathrm{\AA}$]', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, fontsize='x-large')
            ax.step(avg_velocity_axis.value, avg_flux, alpha=0.5, color='blue', where='mid')
            ax.fill_between(avg_velocity_axis.value,
                            (avg_flux-avg_std).value,
                            (avg_flux+avg_std).value,
                            color='blue', alpha=0.15, step='mid')
            
            ax.vlines(0, 0, 1.5, color='black', linestyles='-.', alpha=0.5)
            ax.hlines(1, -vrange, vrange, color='slategrey', linestyles='dashed', alpha=0.9)
            if (gray_out is not None):
                ax[i].fill_betweenx([0,1.5], [gray_out[i][0]]*2, [gray_out[i][1]]*2, color='lightgrey', alpha=0.5)


        if isinstance(calculations, QTable):
            vmin = calculations[calculations['line'] == line_name]['vmin']
            vmax = calculations[calculations['line'] == line_name]['vmax']
            v = calculations[calculations['line'] == line_name]['v']
            if calculations[calculations['line'] == line_name]['interloper_no_measurement'] == False:
                ax.axvline(x=vmin.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
                ax.axvline(x=vmax.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
                ax.axvline(x=v.to('km/s').value, ls='--', ymin=0, ymax=2, color='green')
            if calculations[calculations['line'] == line_name]['default_upper_limit'] == True:
                ax.axvline(x=-60, ls='--', ymin=0, ymax=2, color='red')
                ax.axvline(x=60, ls='--', ymin=0, ymax=2, color='red')
                ax.axvline(x=0, ls='--', ymin=0, ymax=2, color='red')

        plt.xlim(-vrange,vrange)
        plt.ylim(0, 1.5)
        plt.grid(False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False, length=0, width=0, color='none', grid_alpha=0, which='both')

        plt.ylabel('Normalized Flux')
        plt.xlabel('Velocity [km/s]')
