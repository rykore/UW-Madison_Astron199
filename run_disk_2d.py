import os,sys
import numpy as np


common_dir = '../scripts/'
sys.path.append(common_dir)
import nc
import r3d_disk as rdisk
import copy
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

ref_dir = os.path.abspath('./')

# set number of cores to run RADMC3D
ncore  = 2

# load basic disk para and do modification
para_file = 'Sz71_model_params.txt'


sed_obs = np.load('Sz71_obs_SED.npz')

#make modifications from default parameters in MAPS_template
# we use spherical coordinates, nz=0 means no azimuthal difference
r3d_base = {'scattering_mode_max':0,'nphot':'4000000',\
            'nz':0,'sigma_dust_file':'["",""]',\
            'writeimage_unformatted':'1'}
            

disk_dic = rdisk.rac2d_disk_para_2_r3d_dic(\
            fname=para_file, r3d_base=r3d_base)



# setup backup dir 
# the relative path of backupdir to running_dir
backupdir = 'backup/'


# setup disk object
dd = rdisk.setup_disk(disk_dic=disk_dic,sed_obs=sed_obs,\
                      model_name='model_MAPS_template')



argv = sys.argv[1:]



if argv[0] =='run_all':     
    
    dd.setup_dust(use_rac2d_dust=False)
    dd.run_mcthermal(ncore=ncore)
    dd.plot_sed(sed_filename='sed.txt')
    
if argv[0]=='setup_dust': dd.setup_dust(use_rac2d_dust=False)

if argv[0]=='mcthermal': dd.run_mcthermal(ncore=ncore)

if argv[0]=='do_sed': dd.run_sed(ncore=ncore)

if argv[0]=='save_td': dd.save_tdust_2npz()

if argv[0]=='plot_sed':
    if len(argv)>=2:
        sed_filename = argv[1]
    else:
        sed_filename='sed.txt'
        
    dd.plot_sed(sed_filename=sed_filename)

if argv[0]=='plot_dust': dd.plot_dust_analysis()

if argv[0]=='cont_image': dd.do_cont_image(wave_um=1300, sizeau=400,npix=1000)

if argv[0]=='plot_cont': 
    if source=='as209':
        beam = [0.038,0.036,68]
        obs_radial_file = '/mnt/HDDRAID/kezhang/work/mwc480/common_dir/radial_profiles/AS209.profile.txt'

    dd.plot_cont_image(beam=beam,obs_radial_file=obs_radial_file)
        #rdisk.compare_radial_profile(obs_radial_file=obs_radial_file,dpc = disk_dic['dpc'],plot_show=1)

if argv[0] == 'quick_td':
    dd.setup_dust()
    dd.run_mcthermal(ncore=ncore)
    dd.save_tdust_2npz()
    

if argv[0]=='backup': 
    #backupdir=argv[1]
    dd.backup_model(backupdir=backupdir)


if argv[0]=='chi2':
    sed_data = np.loadtxt('GWLup.SED.txt',usecols=(0,1,2,3))

    wavelength_obs = sed_data[:,0]
    fnu_obs        = sed_data[:,1] #Jy
    delta_fnu_obs  = fnu_obs*sed_data[:,3]

    try:
        sed_file = argv[1]
    except:
        print('please put sed file name')

    sed_model = np.loadtxt(sed_file)

    nu_model = nc.cc*1e4/sed_model[:,0] 

    fnu_model = sed_model[:,2]/nu_model*1e23 #Jy

    
    func_interp = interp1d(np.log(sed_model[:,0]),np.log(fnu_model))

    fnu_model_interp = np.exp(func_interp(np.log(wavelength_obs)))

    chi_2 = (fnu_model_interp-fnu_obs)**2/delta_fnu_obs**2

    chi_2_tot = np.sum(chi_2)

    print('chi square value = ', chi_2_tot)

    plt.loglog(wavelength_obs,fnu_obs,'o',label='Data')
    plt.loglog(sed_model[:,0],fnu_model,'r',linewidth=2,label='model')
    plt.text(1,1,sed_file)
    plt.xlabel('wavelength (um)')
    plt.ylabel('fnu(Jy)')

    plt.legend()

    plt.ylim(1e-3,10)

    plt.show()




    

    



    
