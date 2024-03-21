# take CAMB file (e.g. test_lensedCls.dat) and produce dataset with given noise for testing
# Use in cosmomc .ini file using e.g.
# cmb_dataset[MyForecast]=data/forecast_vmodes/prova_06_12.dataset

import shutil
import os
import numpy as np
from getdist import IniFile
from CMBlikes import lastTopComment

def make_forecast_cmb_dataset(input_cl_file, output_root, output_dir=None, lmin=2, lmax=None, fsky=1, fields_use=None, cl_data_cols=''): #NR
    """
    Make a simulated .dataset and associated files with 'data' set at the input fiducial model.

    :param input_cl_file: input fiducial CL
    :param output_root: root name for output files, e.g. 'my_sim1'
    :param output_dir: output directory
    :param lmin: l_min
    :param lmax: l_max
    :param fsky: sky fraction
    :param fields_use: optional list of fields to restict to (e.g. 'T E')
    :param cl_data_cols: if not specified in file header, order of columns in input CL file (e.g. 'TT TE EE BB PP')
    :return:
    """

    ini = IniFile()
    dataset = ini.params

    if not cl_data_cols:
        cl_data_cols = lastTopComment(input_cl_file)
        if not cl_data_cols:
            raise Exception('input CL file must specific names of columns (TT TE EE..)')
    else:
        dataset['cl_hat_order'] = cl_data_cols

    if not fields_use:
        fields_use = ''
        if 'TT' or 'TE' in cl_data_cols: fields_use = 'T'
        if 'EE' or 'TE' in cl_data_cols: fields_use += ' E'

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(__file__), '..', 'data', output_root)
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    dataset['fields_use'] = fields_use

    noise_cols = 'TT           EE'  #NR
    noise_file = output_root + '_Noise.dat'

    noise = np.loadtxt('/marconi_work/INF22_indark/nraffuzz/cosmoMC/data/noise_litebird_only_b30.txt', unpack=True, usecols=(0,1,2))
    noise[1] = noise[1] *noise[0]*(noise[0] + 1)/2/np.pi
    noise[2] = noise[2] *noise[0]*(noise[0] + 1)/2/np.pi
    np.savetxt(os.path.join(output_dir, noise_file), np.column_stack((noise[0],noise[1],noise[2])), header='L ' + noise_cols)

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

    # Pol noise var = ENoiseFac * NoiseVar
    # 2 normally, but for Planck only half detectors are polarized
    input_cl_file = '/marconi_work/INF22_indark/nraffuzz/cosmoMC/data/cl_out_fid.txt' #cl_out_fid.txt
    output_dir = '/marconi_work/INF22_indark/nraffuzz/cosmoMC/data/forecast_vmodes'
    output_root = 'noise_TE1350_serena_sigmatau'
    lmin = 2
    lmax = 1350
    fsky = 0.6
    fields_use = 'T E'

    make_forecast_cmb_dataset(input_cl_file=input_cl_file, output_root=output_root, output_dir=output_dir, 
                              lmin=lmin, lmax=lmax, fsky=fsky, fields_use=fields_use, cl_data_cols=None)
    print('Made ' + os.path.join(output_dir, output_root + '.dataset'))
