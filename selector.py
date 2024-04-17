from header import *

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def get_dvs(zs, z):
    return (c*(zs-z)/(1+z)).to(u.km/u.s)

def get_dls(z, coords, coord):
    radian_to_Mpc = cosmo.angular_diameter_distance(z)
    return radian_to_Mpc*coord.separation(coords).to(u.rad).value

def selectIsolated(dqso_cut, dl_cut, v_cut, vqso_cut, dwarfs, avoid_these, testing=False, test_list=None):
    isolated = Table(names=dwarfs.colnames, dtype=dwarfs.dtype)
    start = time.time()
    
    for i, galaxy in enumerate(dwarfs):
        if test_list is not None:
            if galaxy['SHORTNAME'] in test_list:
                testing = True
            else:
                testing = False

        if testing:
            print('Testing ', galaxy['SHORTNAME'])
       # if i>10:
       #     continue
        # if i % 100 == 0:
        #     print('checked ', i, 'galaxies')
        z = galaxy['REDSHIFT']
        z_QSO = galaxy['REDSHIFT_QSO']

        #############
        
        QSO_coordinate = SkyCoord(galaxy['RA_QSO']*u.degree, galaxy['DEC_QSO']*u.degree)
        coordinate = SkyCoord(galaxy['RA']*u.degree, galaxy['DEC']*u.degree)

        d_qso = get_dls(z, QSO_coordinate, coordinate)

        if d_qso > dqso_cut:
            if testing:
                print('Outside quasar range')
            continue

        ############## 
        avoid_these_temp = avoid_these[avoid_these['NAME'] != galaxy['NAME']]
        ##############

        deltav = get_dvs(avoid_these_temp['REDSHIFT'], galaxy['REDSHIFT'])
        dv_bool = np.abs(deltav) < v_cut

        v_selected = avoid_these_temp[dv_bool]

        ##############
        coordinates = SkyCoord(v_selected['RA']*u.degree, v_selected['DEC']*u.degree)
        ##############

        dls = get_dls(z, coordinates, QSO_coordinate)

        dl_bool = np.abs(dls) < dl_cut

        if dl_bool.sum() > 0:
            if  testing:
                print('Too close to another galaxy')
            continue

        ##############
        dls = get_dls(z, coordinates, coordinate)

        dl_bool = np.abs(dls) < dl_cut

        if dl_bool.sum() > 0:
            if  testing:
                print('Another galaxy to close to quasar')
            continue


        ##############

        qso_deltav = get_dvs(galaxy['REDSHIFT_QSO'], galaxy['REDSHIFT'])
        qso_dv_bool = np.abs(qso_deltav) < vqso_cut

        if qso_dv_bool.sum() > 0:
            if testing:
                print('Too close to quasar redshift')
            continue


        isolated.add_row(galaxy)

    end = time.time()
    #print(end-start, 'seconds')
        ##############


    return isolated
