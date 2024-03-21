# take CAMB file (e.g. test_lensedCls.dat) and produce dataset with given noise for testing
# Use in cosmomc .ini file using e.g.
# cmb_dataset[MyForecast]=data/forecast_vmodes/prova_06_12.dataset

from __future__ import absolute_import
import shutil
import os
import numpy as np
from getdist import IniFile
from CMBlikes import lastTopComment
import math

def make_so_cmb_dataset(input_cl_file, output_root, output_dir=None, beta_degrees=None, so=None,
                              lmin=2, lmax=None, fields_use=None, cl_data_cols=''): #NR
    """
    Make a simulated .dataset and associated files with 'data' set at the input fiducial model.

    :param input_cl_file: input fiducial CL
    :param output_root: root name for output files, e.g. 'my_sim1'
    :param output_dir: output directory
    :param lmin: l_min (not smaller than 40 for SO LAT, not smaller than 2 for SAT)
    :param lmax: l_max (not larger than 7979 for LAT, not larger than 999 for SAT)
    :param fsky: sky fraction for SO LAT is 40%, 10% for SAT
    :param fields_use: optional list of fields to restict to (e.g. 'T E')
    :param cl_data_cols: if not specified in file header, order of columns in input CL file (e.g. 'TT TE EE BB PP')
    :param beta_degrees: departure (phase-shift) from non-ideality of a half-wave plate for circular polarization #NR

    # SO LAT (baseline sensitivity) noise power spectra: post-component separation noise for CMB T
    #  units = uK^2 
    #  Deproj-0: standard ILC
    #  ell, Dlnoise_TT, Dlnoise_EE

    # SO SAT (baseline sensitivity) noise power spectra
    #  B and V modes 
    #  units = uK^2 
    #  L, dl_B, dl_V
    """

    ini = IniFile()
    dataset = ini.params

    if not cl_data_cols:
        cl_data_cols = lastTopComment(input_cl_file)
        if not cl_data_cols:
            raise Exception('input CL file must specific names of columns (TT TE EE..)')
    else:
        dataset['cl_hat_order'] = cl_data_cols

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(__file__), '..', 'data', output_root)
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    if so == 'lat':
        lmin = 50
        lmax = 3000
        fsky = 0.4
        fields_use = 'T E'
        noise_cols = 'TT           EE'
        lmin_so = 40 #for lat the smallest ell is 40
        dl_1, dl_2 = np.loadtxt('/marconi_work/INF24_indark/nraffuzz/cosmoMC/lat_so_noise.txt',unpack=True,usecols=(1,2)) #ell TT EE

    elif so == 'sat':
        lmin = 50
        lmax = 300
        fsky = 0.1
        lmin_so = 2 #for sat the smallest ell is 2

        if beta_degrees == 0:
            fields_use = 'B'
            noise_cols = 'BB'
            beta_degrees = 10**(-9) #precaution to avoid divergencies in the 1/sin(beta)
        else:
            fields_use = 'B V'
            noise_cols = 'BB           VV'
        beta = beta_degrees * np.pi / 180  # degrees to radians
        cos = math.cos(beta/2)**4
        sin = math.sin(beta)**2
        dl_1, dl_2 = np.loadtxt('/marconi_work/INF24_indark/nraffuzz/cosmoMC/sat_so_noise.txt',unpack=True,usecols=(1,2)) #ell BB VV
        dl_1 = dl_1/cos #hwp non ideality for B modes
        dl_2 = dl_2/sin #hwp non ideality for V modes

    noise_file = output_root + '_Noise.dat'

    dataset['fields_use'] = fields_use

    with open(os.path.join(output_dir, noise_file), 'w') as f:
        f.write('#L %s\n' % noise_cols)

        for l in range(lmin, lmax + 1):
            noises = []
            noises += [dl_1[l-lmin_so], dl_2[l-lmin_so]] 
            f.write("%d " % l + " ".join("%E" % elem for elem in noises) + "\n")

    #(l-lmin_so) serve per traslare l'elemento di array in posizione 'l' (non multipolo 'l'),
    #perche' gli ell partono non da 0 ma da lmin_so, percio' lmin > lmin_so

    dataset['fullsky_exact_fsky'] = fsky
    dataset['dataset_format'] = 'CMBLike2'
    dataset['like_approx'] = 'exact'

    dataset['cl_lmin'] = lmin
    dataset['cl_lmax'] = lmax

    dataset['binned'] = False

    dataset['cl_hat_includes_noise'] = False

    shutil.copy(input_cl_file, os.path.join(output_dir, output_root + '.dat'))
    dataset['cl_hat_file'] = output_root + '.dat'
    dataset['cl_noise_file '] = noise_file

    ini.saveFile(os.path.join(output_dir, output_root + '.dataset'))
    print('Saving dataset in data/forecast_vmodes...')

#--------------------------------------------------------------------------------------------------------------# 

if __name__ == "__main__":
    import tempfile

    delens = '1' # 1 or 0.3
    input_cl_file = '/marconi_work/INF24_indark/nraffuzz/cosmoMC/data/fiducial_delensing'+delens+'.txt' #cl_out_fid.txt
    output_dir = '/marconi_work/INF24_indark/nraffuzz/cosmoMC/data/forecast_vmodes/3_delens'+delens
    beta_degrees = 10
    output_root = '3_BV_sosat_b%.0f'%(beta_degrees) #'3_B_sosat_delens0p3_b%.0f'%(beta_degrees)
    so = 'sat' # sat or lat

    make_so_cmb_dataset(input_cl_file=input_cl_file, output_root=output_root, output_dir=output_dir, 
                              beta_degrees=beta_degrees, so=so, cl_data_cols=None)
    print('Made ' + os.path.join(output_dir, output_root + '.dataset'))
    #lmin=lmin, lmax=lmax,