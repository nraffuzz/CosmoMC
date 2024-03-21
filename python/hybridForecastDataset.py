# take CAMB file (e.g. test_lensedCls.dat) and produce dataset with given noise for testing
# Use in cosmomc .ini file using e.g.
# cmb_dataset[MyForecast]=data/forecast_vmodes/prova_06_12.dataset

from __future__ import absolute_import
import shutil
import os
import numpy as np
import healpy as hp
from getdist import IniFile
from CMBlikes import lastTopComment


def make_forecast_cmb_dataset(input_cl_file, output_root, output_dir=None, beta_degrees=None, noise_muK_arcmin_T_exp1=None,
                              noise_muK_arcmin_P_exp1=None, noise_muK_arcmin_T_exp2=None, noise_muK_arcmin_P_exp2=None,
                              fwhm_arcmin_exp1=None, fwhm_arcmin_exp2=None, lmin=2, lmax=None, fsky=1, fields_use=None,
                              lens_recon_noise=None, cl_data_cols=''):
    """
    Make a simulated .dataset and associated files with 'data' set at the input fiducial model.

    :param input_cl_file: input fiducial CL
    :param output_root: root name for output files, e.g. 'my_sim1'
    :param output_dir: output directory
    :param noise_muK_arcmin_T: temperature noise in muK-arcmin
    :param noise_muK_arcmin_P: polarization noise in muK-arcmin
    :param NoiseVar: effective isotropic noise variance for the temperature (N_L=NoiseVar with no beam)
    :param ENoiseFac: factor by which polarization noise variance is higher (usually 2, for Planck about 4
                        as only half the detectors polarized)
    :param fwhm_arcmin: beam fwhm in arcminutes
    :param lmin: l_min
    :param lmax: l_max
    :param fsky: sky fraction
    :param fields_use: optional list of fields to restict to (e.g. 'T E')
    :param lens_recon_noise: optional array, starting at L=0, for the PP lensing reconstruction noise, in [L(L+1)]^2C_L^phi/2pi units
    :param cl_data_cols: if not specified in file header, order of columns in input CL file (e.g. 'TT TE EE BB PP')
    :param beta_degrees: departure (phase-shift) from non-ideality of a half-wave plate for circular polarization #NR
    :return:
    """

    use_lensing = lens_recon_noise
    use_CMB = noise_muK_arcmin_T_exp1 is not None

    ini = IniFile()
    dataset = ini.params

    if not cl_data_cols:
        cl_data_cols = lastTopComment(input_cl_file)
        if not cl_data_cols:
            raise Exception('input CL file must specific names of columns (TT TE EE..)')
    else:
        dataset['cl_hat_order'] = cl_data_cols

    if use_CMB:
        NoiseVar_exp1 = (noise_muK_arcmin_T_exp1 * np.pi / 180 / 60.) ** 2 # experiment 1
        NoiseVar_exp2 = (noise_muK_arcmin_T_exp2 * np.pi / 180 / 60.) ** 2 # experiment 2
        if noise_muK_arcmin_P_exp1 or noise_muK_arcmin_P_exp2 is not None:
            ENoiseFac_exp1 = (noise_muK_arcmin_P_exp1 / noise_muK_arcmin_T_exp1) ** 2
            ENoiseFac_exp2 = (noise_muK_arcmin_P_exp2 / noise_muK_arcmin_T_exp2) ** 2
        if not fields_use:
            fields_use = ''
            if 'TT' or 'TE' in cl_data_cols: fields_use = 'T'
            if 'EE' or 'TE' in cl_data_cols: fields_use += ' E'
            if 'BB' in cl_data_cols: fields_use += ' B'
            if 'VV' in cl_data_cols: fields_use += ' V'  #NR
            if 'PP' in cl_data_cols and use_lensing: fields_use += ' P'
    else:
        fields_use = fields_use or 'P'

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(__file__), '..', 'data', output_root)
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    dataset['fields_use'] = fields_use

    if use_CMB:
        fwhm_exp1 = fwhm_arcmin_exp1 / 60 * np.pi/180 #arcmin --> degrees --> radians for experiment 1
        fwhm_exp2 = fwhm_arcmin_exp2 / 60 * np.pi/180 #arcmin --> degrees --> radians for experiment 2
        if beta_degrees == 0:
            beta_degrees = 10**(-9)
        beta = beta_degrees * np.pi / 180  #degrees to radians
        noise_cols = 'TT           EE'  #NR
        if use_lensing: noise_cols += '          PP'
    elif use_lensing:
        noise_cols = 'PP'
    noise_file = output_root + '_Noise.dat'
    bl2_exp1 = (hp.gauss_beam(fwhm_exp1, lmax = lmax, pol = True))**2 #beam window func bl(ell,colonne)
    bl2_exp2 = (hp.gauss_beam(fwhm_exp2, lmax = lmax, pol = True))**2

    fac_exp1 = NoiseVar_exp1
    fac_exp2 = NoiseVar_exp2
    cos = np.cos(beta/2)**4 #to speed up the process I store the cos or sin values here
    #sin = np.sin(beta)**2
    
    with open(os.path.join(output_dir, noise_file), 'w') as f:
        f.write('#L %s\n' % noise_cols)

        for l in range(lmin, lmax + 1):
            Nl_hybrid = []
            Nl_hybrid += [
                l*(l+1.)/2/np.pi * (bl2_exp1[l,0] /fac_exp1 + bl2_exp2[l,0] /fac_exp2 )**(-1),
                l*(l+1.)/2/np.pi * (bl2_exp1[l,1] / (fac_exp1 * ENoiseFac_exp1) + bl2_exp2[l,1] * cos /(fac_exp2 * ENoiseFac_exp2) )**(-1)
                ]
            if use_lensing: Nl_hybrid += [lens_recon_noise[l]]
            f.write("%d " % l + " ".join("%E" % elem for elem in Nl_hybrid) + "\n")

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

if __name__ == "__main__":
    import tempfile

    # Noise var is N_l in muK^2 for white noise
    # note  NoiseVar = (muKArcmin * np.pi / 180 / 60.) ** 2
    # Pol noise var = ENoiseFac * NoiseVar
    # ENoiseFac = 2 normally, but 4 for Planck only half detectors are polarized

    #noise, exp2 is the one using HWP, COMPLETE HERE
    noise_muK_arcmin_T_exp1 = 3 # uK arcmin
    noise_muK_arcmin_P_exp1 = np.sqrt(2)*noise_muK_arcmin_T_exp1
    noise_muK_arcmin_T_exp2 = 2 # uK arcmin
    noise_muK_arcmin_P_exp2 = np.sqrt(2)*noise_muK_arcmin_T_exp2
    
    #beam
    fwhm_arcmin_exp1 = 1.5 #arcmin
    fwhm_arcmin_exp2 = 30  #arcmin

    delens = ['1','0p3'] # 1 or 0.3 #COMPLETE HERE
    root_dir = '/marconi_work/INF24_indark/nraffuzz/cosmoMC/data'
    
    lmin = 51    #COMPLETE HERE
    lmax = 3000  #COMPLETE HERE
    fsky = 0.6   #COMPLETE HERE
    prefix = 'S4LAT_LB_' #COMPLETE HERE
    which_dataset = '5'  #COMPLETE HERE
    fields_use         = [      'T E',          'T E'] #COMPLETE HERE
    use_non_ideal_hwp  = [      False,           True] #COMPLETE HERE
    name_distinguisher = ['_idealHWP', '_nonidealHWP'] #COMPLETE HERE
    for j in range(len(delens)):
        for i in range(len(use_non_ideal_hwp)):
            if use_non_ideal_hwp[i]:
                beta_degrees = 10
            else:
                beta_degrees = 0
            
            if 'V' in fields_use[i] and use_non_ideal_hwp[i]==False:
                raise Exception('Must use non-ideal half-wave plate if using V-modes')
            
            input_cl_file = root_dir+'/fiducial_delensing'+delens[j]+'.txt' #cl_out_fid.txt
            output_dir = root_dir+'/forecast_vmodes/'+which_dataset+'_delens'+delens[j]
            
            output_root = prefix+fields_use[i].replace(' ', '')+'_ell'+str(lmin)+'_'+str(lmax)\
                +'_delens'+delens[j]+'_fsky_p'+str(fsky)[2:]+name_distinguisher[i] #COMPLETE HERE
            
            make_forecast_cmb_dataset(
                input_cl_file=input_cl_file, 
                output_root=output_root, 
                output_dir=output_dir, 
                beta_degrees=beta_degrees, 
                noise_muK_arcmin_T_exp1=noise_muK_arcmin_T_exp1, 
                noise_muK_arcmin_P_exp1=noise_muK_arcmin_P_exp1, 
                noise_muK_arcmin_T_exp2=noise_muK_arcmin_T_exp2, 
                noise_muK_arcmin_P_exp2=noise_muK_arcmin_P_exp2, 
                fwhm_arcmin_exp1=fwhm_arcmin_exp1, 
                fwhm_arcmin_exp2=fwhm_arcmin_exp2, 
                lmin=lmin, 
                lmax=lmax, 
                fsky=fsky, 
                fields_use=fields_use[i]
                )
            print('Made ' + os.path.join(output_dir, output_root + '.dataset'))
