# take CAMB file (e.g. test_lensedCls.dat) and produce dataset with given noise for testing
# Use in cosmomc .ini file using e.g.
# cmb_dataset[MyForecast]=data/forecast_vmodes/prova_06_12.dataset

from __future__ import absolute_import
import shutil
import os
import numpy as np
from getdist import IniFile
from CMBlikes import lastTopComment#, DatasetLikelihood, ClsArray
import math

def make_forecast_cmb_dataset(input_cl_file, output_root, output_dir=None, beta_degrees=None, noise_muK_arcmin_T=None,
                              noise_muK_arcmin_P=None, NoiseVar=None, ENoiseFac=2, fwhm_arcmin=None,
                              lmin=2, lmax=None, fsky=1, fields_use=None,
                              lens_recon_noise=None, cl_data_cols=''): #NR
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
    use_CMB = noise_muK_arcmin_T or NoiseVar is not None

    ini = IniFile()
    dataset = ini.params

    if not cl_data_cols:
        cl_data_cols = lastTopComment(input_cl_file)
        if not cl_data_cols:
            raise Exception('input CL file must specific names of columns (TT TE EE..)')
    else:
        dataset['cl_hat_order'] = cl_data_cols

    if use_CMB:
        if NoiseVar is None:
            if noise_muK_arcmin_T is None:
                raise ValueError('Must specify noise')
            NoiseVar = (noise_muK_arcmin_T * np.pi / 180 / 60.) ** 2  
            if noise_muK_arcmin_P is not None:
                ENoiseFac = (noise_muK_arcmin_P / noise_muK_arcmin_T) ** 2
        elif noise_muK_arcmin_T is not None or noise_muK_arcmin_P is not None:
            raise ValueError('Specific either noise_muK_arcmin or NoiseVar')
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
        fwhm = fwhm_arcmin / 60 #degrees
        if beta_degrees == 0:
            beta_degrees = 10**(-9)
        beta = beta_degrees * np.pi / 180  #NR degrees to radians
        xlc = 180 * np.sqrt(8. * np.log(2.)) / np.pi 
        sigma2 = (fwhm / xlc) ** 2
        cos = math.cos(beta/2)**4
        sin = math.sin(beta)**2
        noise_cols = 'TT           EE          BB          VV'  #NR
        if use_lensing: noise_cols += '          PP'
    elif use_lensing:
        noise_cols = 'PP'
    noise_file = output_root + '_Noise.dat'
    with open(os.path.join(output_dir, noise_file), 'w') as f:
        f.write('#L %s\n' % noise_cols)

        for l in range(lmin, lmax + 1):
            fac = l*(l+1.)/2/np.pi * np.exp(l*(l+1)*sigma2) * NoiseVar #NR np.exp(l*(l+1)*sigma2) this is already beam squared
            NoiseCl = [fac, fac/cos, fac/cos, fac/sin] #NR 
            noises = []
            if use_CMB: noises += [NoiseCl[0], ENoiseFac * NoiseCl[1], ENoiseFac * NoiseCl[2], ENoiseFac * NoiseCl[3]] #NR
            if use_lensing: noises += [lens_recon_noise[l]]
            f.write("%d " % l + " ".join("%E" % elem for elem in noises) + "\n")

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
    # 2 normally, but for Planck only half detectors are polarized

    delens = ['1','0p3'] # 1 or 0.3
    root_dir = '/marconi_work/INF24_indark/nraffuzz/cosmoMC/data'
    lmin = 50    #COMPLETE HERE
    lmax = 300   #COMPLETE HERE
    fsky = 0.03  #COMPLETE HERE
    noise_muK_arcmin_T = 1 #COMPLETE HERE
    fwhm_arcmin = 30       #COMPLETE HERE
    prefix = 'S4SAT'       #COMPLETE HERE
    which_dataset = '5'    #COMPLETE HERE
    fields_use         = [        'B',          'B V'] #COMPLETE HERE # 'T E B V' 
    use_non_ideal_hwp  = [      False,           True]
    name_distinguisher = ['_idealHWP', '_nonidealHWP']
    

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
            output_root = prefix+'_'+fields_use[i].replace(' ', '')+'_ell'+str(lmin)+'_'+str(lmax)+'_delens'\
                +delens[j]+'_fsky_p'+str(fsky)[2:]+name_distinguisher[i]
            noise_muK_arcmin_P = np.sqrt(2)*noise_muK_arcmin_T # np.sqrt(2) or 2 (for Planck, as half dets are polarized) 
            
            make_forecast_cmb_dataset(
                input_cl_file=input_cl_file, 
                output_root=output_root, 
                output_dir=output_dir, 
                beta_degrees=beta_degrees, 
                noise_muK_arcmin_T=noise_muK_arcmin_T, 
                noise_muK_arcmin_P=noise_muK_arcmin_P, 
                NoiseVar=None,
                ENoiseFac=None, 
                fwhm_arcmin=fwhm_arcmin, 
                lmin=lmin, 
                lmax=lmax, 
                fsky=fsky, 
                fields_use=fields_use[i],
                lens_recon_noise=None, 
                cl_data_cols=None
                )
            print('Made ' + os.path.join(output_dir, output_root + '.dataset'))

            # The rest is just a test on files produced above
            #like = DatasetLikelihood(os.path.join(output_dir, output_root + '.dataset'))
            #cls = ClsArray(input_cl_file) # input_cl_file lensedTotClFileRoot
            #cls.cls_array[0, 0] *= 1.004
            #cls.cls_array[1, 1] *= 0.991
            #import time

            #start = time.time()
            #chi2 = like.chi_squared(cls)
            #end = time.time() - start
            #print('Test chi2 = %s' % chi2)
            #print('Time: %s' % end)
            #if not np.allclose(49.055, chi2, rtol=1e-5): raise Exception('likelihood test failed')