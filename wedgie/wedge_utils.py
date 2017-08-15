"""
This module contains functions to create and plot wedges/pitchforks from HERA data files.

Authors:
Paul Chichura <pchich_at_sas.upenn.edu>
Austin Fox Fortino <fortino_at_sas.upenn.edu>
Amy Igarashi <igarashiamy_at_gmail.com>
Saul Aryeh Kohn <saulkohn_at_sas.upenn.edu>
"""

import capo, aipy, decimal
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
import scipy.constants as sc
import gen_utils as gu
import cosmo_utils as cu
import matplotlib.image as mpimg

# For Interactive Development
from IPython import embed

# Repeated Operations for Wedge Creation:
def cleanfft(d, f, pair, pol, clean):
    """Performs the fft with an antenna pair for a certain polarization with a blackman-harris window, then performs the CLEAN."""

    # fft
    w = aipy.dsp.gen_window(d[pair][pol].shape[-1], window='blackman-harris')
    _dw = np.fft.ifft(d[pair][pol] * w)

    # CLEAN
    _ker= np.fft.ifft(f[pair][pol] * w)
    gain = aipy.img.beam_gain(_ker)
    for time in range(_dw.shape[0]):
        _dw[time, :], info = aipy.deconv.clean(_dw[time, :], _ker[time, :], tol=clean)
        _dw[time, :] += info['res'] / gain

    ftd_2D_data = np.ma.array(_dw)

    return ftd_2D_data

def name_npz(args, files, pol, ex_ants, tag):
    """Formats the name of the npz file depending on which type of wedge is being made."""
    args = vars(args)

    file_start = files[0].split('/')[-1].split('.')
    file_end = files[-1].split('/')[-1].split('.')

    zen_day = file_start[:2]
    time_range = ['{}_{}'.format(file_start[2], file_end[2])]
    HH = [file_start[4]]
    ext = [file_start[5]]
    pol = [pol]
    tag = [tag]

    if ex_ants:
        ex_ants = ['_'.join(args['ex_ants'].split(','))]
    else:
        ex_ants = []

    npz_name = zen_day + time_range + pol + HH + ex_ants + ext + tag + ["npz"]
    npz_name = '.'.join(npz_name)
    print npz_name

    return npz_name

# Calfile Specific Operations:
def calculate_slope(antennae, pair):
    """Returns the slope between two antennae given the dictionary of antennae and their positions and the antenna pair in question"""
    decimal.getcontext().prec = 6
    dy = decimal.Decimal(antennae[pair[1]]['top_y'] - antennae[pair[0]]['top_y'])
    dx = decimal.Decimal(antennae[pair[1]]['top_x'] - antennae[pair[0]]['top_x'])
    if dx != 0:
        slope = float(dy / dx)
    else:
        slope = np.inf

    return slope

def calculate_baseline(antennae, pair):
    """Returns the baseline length between two antennae given the dictionary of antennae and their positions and the antenna pair in question"""
    decimal.getcontext().prec = 6
    dx = decimal.Decimal(antennae[pair[0]]['top_x'] - antennae[pair[1]]['top_x'])
    dy = decimal.Decimal(antennae[pair[0]]['top_y'] - antennae[pair[1]]['top_y'])
    baseline = (dx**2 + dy**2).sqrt()
    
    return float(baseline)

def get_baselines(calfile, ex_ants):
    """
    Returns a dictionary of baseline lengths and the corresponding pairs. The data is based 
    on a calfile. ex_ants is a list of integers that specify antennae to be exlcuded from 
    calculation.
    
    Requires cal file to be in PYTHONPATH.
    """
    try:
        print 'Reading calfile: %s.' %calfile
        exec("import {cfile} as cal".format(cfile=calfile))
        antennae = cal.prms['antpos_ideal']
    except ImportError:
        raise Exception("Unable to import {cfile}.".format(cfile=calfile))
    
    # Remove all placeholder antennae from consideration
    # Remove all antennae from ex_ants from consideration
    ants = []
    for ant in antennae.keys():
        if (not antennae[ant]['top_z'] < 0) and (ant not in ex_ants):
            ants.append(ant)

    # Form pairs of antennae
    # Store unique baselines and slopes for later use
    pairs, baselines, slopes =  {}, [], []
    for ant_i in ants:
        for ant_j in ants:
            if (ant_i >= ant_j): continue
            pair = (ant_i, ant_j)

            baseline = calculate_baseline(antennae, pair)
            baselines.append(baseline)
            
            slope = calculate_slope(antennae, pair)
            slopes.append(slope)
            
            pairs[pair] = (baseline, slope)

    # Remove duplicates baseline and slope values
    baselines = set(baselines)
    slopes = set(slopes)

    # Initalize antdict with baselines as keys and empty lists as values
    antdict = {baseline: [] for baseline in baselines}

    # Add pairs to the list of their respective baseline
    for pair in pairs:
        baseline = pairs[pair][0]
        antdict[baseline].append(pair)

    # Initialize slopedict with baselines for keys and the dictionary of slopes for each value
    slopedict = {}
    for baseline in baselines:
        slopedict[baseline] = {slope: [] for slope in slopes}

    # Add pairs to their respective slope within their respective baseline
    for pair in pairs:
        baseline = pairs[pair][0]
        slope = pairs[pair][1]
        slopedict[baseline][slope].append(pair)

    for baseline in slopedict.copy():
        for slope in slopedict[baseline].copy():
            if slopedict[baseline][slope] == []:
                del slopedict[baseline][slope]

    return (antdict, slopedict, pairs, sorted(list(baselines)), sorted(list(slopes)))

# Forming Stokes Paramters:
def stokesI(files):
    """Calculate I (VI = Vxx + Vyy)"""
    txx, dxx, fxx = capo.miriad.read_files(files['xx'], antstr='cross', polstr='xx')
    tyy, dyy, fyy = capo.miriad.read_files(files['yy'], antstr='cross', polstr='yy')

    tI = txx
    dI = {}
    fI = {}
    for key in dxx.keys():
        dI[key] = {'I': dxx[key]['xx'] + dyy[key]['yy'] }
        fI[key] = {'I': fxx[key]['xx'] + fyy[key]['yy'] }

    del txx, dxx, fxx
    del tyy, dyy, fyy

    return tI, dI, fI

def stokesQ(files):
    """Calculate Q (VQ = Vxx - Vyy)"""
    txx, dxx, fxx = capo.miriad.read_files(files['xx'], antstr='cross', polstr='xx')
    tyy, dyy, fyy = capo.miriad.read_files(files['yy'], antstr='cross', polstr='yy')

    tQ = tyy
    dQ = {}
    fQ = {}
    for key in dxx.keys():
        dQ[key] = {'Q': dxx[key]['xx'] - dyy[key]['yy'] }
        fQ[key] = {'Q': fxx[key]['xx'] + fyy[key]['yy'] }

    del txx, dxx, fxx
    del tyy, dyy, fyy

    return tQ, dQ, fQ

def stokesU(files):
    """Calculate U (VU = Vxy + Vyx)"""
    tyx, dyx, fyx = capo.miriad.read_files(files['yx'], antstr='cross', polstr='yx')
    txy, dxy, fxy = capo.miriad.read_files(files['xy'], antstr='cross', polstr='xy')

    tU = tyx
    dU = {}
    fU = {}
    for key in dxy.keys():
        dU[key] = {'U': dxy[key]['xy'] + dyx[key]['yx'] }
        fU[key] = {'U': fxy[key]['xy'] + fyx[key]['yx'] }

    del tyx, dyx, fyx
    del txy, dxy, fxy

    return tU, dU, fU

def stokesV(files):
    """Calculate V (VV = -i*Vxy + i*Vyx)"""
    tyx, dyx, fyx = capo.miriad.read_files(files['yx'], antstr='cross', polstr='yx')
    txy, dxy, fxy = capo.miriad.read_files(files['xy'], antstr='cross', polstr='xy')

    tV = txy
    dV = {}
    fV = {}
    for key in dxy.keys():
        dV[key] = {'V': -1j*dxy[key]['xy'] + 1j*dyx[key]['yx'] }
        fV[key] = {'V': fxy[key]['xy'] + fyx[key]['yx'] }

    del tyx, dyx, fyx
    del txy, dxy, fxy

    return tV, dV, fV

# Data Analysis Functions:
def wedge_timeavg(args, files, pol, calfile, history, freq_range, ex_ants):
    npz_name = name_npz(args, files[files.keys()[0]], pol, ex_ants, 'timeavg')

    if pol in ['I', 'Q', 'U', 'V']:
        exec('t, d, f = stokes{}(files)'.format(pol))
    elif pol in ['xx', 'xy', 'yx', 'yy']:
        t, d, f = capo.miriad.read_files(files, antstr='cross', polstr=pol)
    else:
        raise Exception('Polarization Not Recognized')

    t['freqs'] = t['freqs'][freq_range[0]: freq_range[1]]

    for key in d.keys():
        d[key][pol] = d[key][pol][:, freq_range[0]: freq_range[1]]
        f[key][pol] = f[key][pol][:, freq_range[0]: freq_range[1]]

    #stores vis^2 for each baselength averaged over time
    wedgeslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = get_baselines(calfile, ex_ants)[0]
    baselengths = antdict.keys()
    baselengths.sort()

    #get number of times and number of channels
    ntimes, nchan = len(t['times']), len(t['freqs'])

    dt = np.diff(t['times'])[0] #dJD

    lst_range = []

    #get vis^2 for each baselength
    for baselength in baselengths:
        vissq_per_bl = np.zeros((ntimes // 2 ,nchan), dtype=np.complex128)
        
        #go through each individual antenna pair
        for antpair in antdict[baselength]:
            vissq_per_antpair = np.zeros((ntimes // 2 ,nchan), dtype=np.complex128)

            # Create Metadata    
            uv = aipy.miriad.UV(files[files.keys()[0]][0])
            aa = aipy.cal.get_aa(calfile, uv['sdf'], uv['sfreq'], uv['nchan']) 
            del(uv)

            # FFT and CLEAN
            ftd_2D_data = cleanfft(d, f, antpair, pol, clean=1e-3)

            #multiply at times (1*2, 3*4, etc...) 
            for i in range(ntimes):
                 
                #set up phasing    
                aa.set_active_pol(pol) 
                if i!=0:
                    old_zenith = zenith
                time = t['times'][i]    
                aa.set_jultime(time)
                lst = aa.sidereal_time()
                zenith = aipy.phs.RadioFixedBody(lst, aa.lat)
                zenith.compute(aa)
                if i==0 and baselength==baselengths[0] and antpair==antdict[baselength][0]:
                    lst_range.append(str(lst))
                if i==ntimes-1 and baselength==baselengths[0] and antpair==antdict[baselength][0]:
                    lst_range.append(str(lst))
                if i==0:
                    continue

                #phase and multiply, store in vissq_per_antpair
                if i % 2:
                    _v1 = ftd_2D_data[i-1,:]
                    phase_correction = np.conj(aa.gen_phs(zenith,antpair[0],antpair[1]))*aa.gen_phs(old_zenith,antpair[0],antpair[1])
                    _v2 = ftd_2D_data[i,:]*phase_correction[freq_range[0]:freq_range[1]]
                    vissq_per_antpair[i // 2,:] = np.conj(_v1)*_v2

            #store time average for this baseline length
            vissq_per_bl += vissq_per_antpair

        #compute average for baseline length, average over time, and store in wedgeslices
        vissq_per_bl /= len(antdict[baselength])
        wedgeslices.append(np.log10(np.fft.fftshift(np.abs(np.mean(vissq_per_bl, axis=0)))))
        print 'Wedgeslice for baseline {} complete.'.format(baselength)

    #get delays
    channel_width = (t['freqs'][1] - t['freqs'][0])*10**3 # Channel width in units of GHz
    num_bins = len(t['freqs'])
    delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width / num_bins))
    
    np.savez(npz_name, wdgslc=wedgeslices, dlys=delays, pol=pol, bls=baselengths, lst=lst_range, hist=history)
    return npz_name

def wedge_blavg(args, files, pol, calfile, history, freq_range, ex_ants):
    npz_name = name_npz(args, files[files.keys()[0]], pol, ex_ants, 'blavg')

    if pol in ['I', 'Q', 'U', 'V']:
        exec('t, d, f = stokes{}(files)'.format(pol))
    elif pol in ['xx', 'xy', 'yx', 'yy']:
        t, d, f = capo.miriad.read_files(files, antstr='cross', polstr=pol)
    else:
        raise Exception('Polarization Not Recognized')

    t['freqs'] = t['freqs'][freq_range[0]:freq_range[1]]
    ntimes,nchan = len(t['times']),len(t['freqs'])

    for key in d.keys():
        d[key][pol] = d[key][pol][:,freq_range[0]:freq_range[1]]
        f[key][pol] = f[key][pol][:,freq_range[0]:freq_range[1]]

    #create variable to store the wedge subplots in
    wedgeslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = get_baselines(calfile, ex_ants)[0]
    baselengths = antdict.keys()
    baselengths.sort()
    
    lst_range = []
    
    #for each baselength in the dictionary
    for length in baselengths:
        
        totald = np.zeros_like(d[antdict[length][0]][pol], dtype=np.complex128) #variable to store cumulative fft data
        vissq_per_bl = np.zeros((ntimes // 2, nchan), dtype=np.complex128)

        #cycle through every ant pair of the given baselength
        for antpair in antdict[length]:

            #create/get metadata    
            uv = aipy.miriad.UV(files[files.keys()[0]][0])
            aa = aipy.cal.get_aa(calfile, uv['sdf'], uv['sfreq'], uv['nchan']) 
            del(uv)


            #CLEAN and fft the data
            totald += cleanfft(d, f, antpair, pol, clean=1e-3)

            #holds our data
            vissq_per_antpair = np.zeros((ntimes // 2 ,nchan), dtype=np.complex128)

            #multiply at times (1*2, 3*4, etc...) 
            for i in range(ntimes):
                 
                #set up phasing    
                aa.set_active_pol(pol) 
                if i!=0:
                    old_zenith = zenith
                time = t['times'][i]    
                aa.set_jultime(time)
                lst = aa.sidereal_time()
                zenith = aipy.phs.RadioFixedBody(lst, aa.lat)
                zenith.compute(aa)
                if i==0 and length==baselengths[0] and antpair==antdict[length][0]:
                    lst_range.append(str(lst))
                if i==ntimes-1 and length==baselengths[0] and antpair==antdict[length][0]:
                    lst_range.append(str(lst))
                if i==0:
                    continue

                #phase and multiply, store in vissq_per_antpair
                if i % 2:
                    _v1 = totald[i-1,:]
                    phase_correction = np.conj(aa.gen_phs(zenith,antpair[0],antpair[1]))*aa.gen_phs(old_zenith,antpair[0],antpair[1])
                    _v2 = totald[i,:]*phase_correction[freq_range[0]:freq_range[1]]
                    vissq_per_antpair[i // 2,:] = np.conj(_v1)*_v2

            #store time average for this baseline length
            vissq_per_bl += vissq_per_antpair
        
        
        #get average of all values for this baselength, store in wedgeslices
        vissq_per_bl /= len(antdict[length])
        wedgeslices.append(np.log10(np.fft.fftshift(np.abs(vissq_per_bl), axes=1)))
        print 'baseline {} complete.'.format(length)

    channel_width = (t['freqs'][1] - t['freqs'][0])*10**3 # Channel width in units of GHz
    num_bins = len(t['freqs'])
    #get data to recalculate axes  
    delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width / num_bins))

    #save filedata as npz
    #NB: filename of form like "zen.2457746.16693.xx.HH.uvcOR"
    np.savez(npz_name, times=t['times'], wdgslc=wedgeslices, dlys=delays, pol=pol, bls=baselengths, lst=lst_range, hist=history)
    return npz_name

def wedge_flavors(args, files, pol, calfile, history, freq_range, ex_ants):
    npz_name = name_npz(args, files[files.keys()[0]], pol, ex_ants, 'flavors')

    if pol in ['I', 'Q', 'U', 'V']:
        exec('t, d, f = stokes{}(files)'.format(pol))
    elif pol in ['xx', 'xy', 'yx', 'yy']:
        t, d, f = capo.miriad.read_files(files, antstr='cross', polstr=pol)
    else:
        raise Exception('Polarization Not Recognized')

    t['freqs'] = t['freqs'][freq_range[0]: freq_range[1]]
    for key in d.keys():
        d[key][pol] = d[key][pol][:, freq_range[0]:freq_range[1]]
        f[key][pol] = f[key][pol][:, freq_range[0]:freq_range[1]]

    baseline_info = get_baselines(calfile, ex_ants)
    antdict, slopedict, pairs = baseline_info[0], baseline_info[1], baseline_info[2]
    baselines, slopes = baseline_info[3], baseline_info[4]

    ntimes, nchan = len(t['times']), len(t['freqs'])

    dt = np.diff(t['times'])[0]

    lst_range = []

    wedgeslices = []
    for baseline in sorted(slopedict.keys()):
        vis_sq_baseline = np.zeros((ntimes // 2 ,nchan), dtype=np.complex128)
        for slope in sorted(slopedict[baseline].keys()):
            vis_sq_slope = np.zeros((ntimes // 2 ,nchan), dtype=np.complex128)
            for pair in slopedict[baseline][slope]:
                vis_sq_pair = np.zeros((ntimes // 2, nchan), dtype=np.complex128)

                uv = aipy.miriad.UV(files[files.keys()[0]][0])
                aa = aipy.cal.get_aa(calfile.split('.')[0], uv['sdf'], uv['sfreq'], uv['nchan'])
                del(uv)
                
                ftd_2D_data = cleanfft(d, f, pair, pol, clean=1e-3)

                for i in range(ntimes):                     
                    aa.set_active_pol(pol) 
                    if i != 0:
                        old_zenith = zenith
                    time = t['times'][i]    
                    aa.set_jultime(time)
                    lst = aa.sidereal_time()
                    zenith = aipy.phs.RadioFixedBody(lst, aa.lat)
                    zenith.compute(aa)
                    if i==0 and baseline==sorted(slopedict.keys())[0] and slope==sorted(slopedict[baseline])[0] and pair==slopedict[baseline][slope][0]:
                        lst_range.append(str(lst))
                    if i==ntimes-1 and baseline==sorted(slopedict.keys())[0] and slope==sorted(slopedict[baseline])[0] and pair==sorted(slopedict[baseline][slope])[0]:
                        lst_range.append(str(lst))
                    if i == 0:
                        continue

                    if i % 2:
                        _v1 = ftd_2D_data[i-1,:]
                        phase_correction = np.conj(aa.gen_phs(zenith, pair[0], pair[1])) * aa.gen_phs(old_zenith, pair[0], pair[1])
                        _v2 = ftd_2D_data[i, :] * phase_correction[freq_range[0]: freq_range[1]]

                        vis_sq_pair[i // 2, :] = np.conj(_v1) * _v2

                vis_sq_slope += vis_sq_pair

            vis_sq_slope /= len(slopedict[baseline])

            wedgeslices.append(np.log10(np.fft.fftshift(np.abs(np.mean(vis_sq_slope, axis=0)))))
            print 'Wedgeslice for baseline {} and slope {} complete.'.format(baseline, slope)
   
    channel_width = (t['freqs'][1] - t['freqs'][0])*10**3 # Channel width in units of GHz
    num_bins = len(t['freqs'])
    delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width / num_bins))
    
    np.savez(npz_name, wdgslc=wedgeslices, dlys=delays, pol=pol, bls=baselines, slps=slopes, prs=pairs, slpdct=slopedict, lst=lst_range, hist=history)
    return npz_name

def wedge_bltype(args, files, pol, calfile, history, freq_range, ex_ants):
    bl_num = args.bl_num
    npz_name = name_npz(args, files[files.keys()[0]], pol, ex_ants, 'bl{}'.format(bl_num))


    if pol in ['I', 'Q', 'U', 'V']:
        exec('t, d, f = stokes{}(files)'.format(pol))
    elif pol in ['xx', 'xy', 'yx', 'yy']:
        t, d, f = capo.miriad.read_files(files, antstr='cross', polstr=pol)
    else:
        raise Exception('Polarization Not Recognized')

    bl_num -= 1

    t['freqs'] = t['freqs'][freq_range[0]:freq_range[1]]
    ntimes,nchan = len(t['times']),len(t['freqs'])

    for key in d.keys():
        d[key][pol] = d[key][pol][:,freq_range[0]:freq_range[1]]
        f[key][pol] = f[key][pol][:,freq_range[0]:freq_range[1]]

    #create variable to store the wedge subplots in
    antpairslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = get_baselines(calfile, ex_ants)[0]
    baselengths = antdict.keys()
    baselengths.sort()

    totald = np.zeros_like(d[antdict[baselengths[bl_num]][0]][pol])

    length = baselengths[bl_num]

    #access antenna tuples in baselengths dictionary, antpair is antenna tuple
    for antpair in antdict[baselengths[bl_num]]:
        
        #an array to store visibilities^2 per antpair
        vissq_per_antpair = np.zeros((ntimes // 2 ,nchan), dtype=np.complex128)
        
        #create/get metadata    
        uv = aipy.miriad.UV(files[files.keys()[0]][0])
        aa = aipy.cal.get_aa(calfile.split('.')[0], uv['sdf'], uv['sfreq'], uv['nchan']) 
        del(uv)

        #CLEAN and fft the data
        totald = cleanfft(d, f, antpair, pol, clean=1e-3)


        #multiply at times (1*2, 3*4, etc...) 
        for i in range(ntimes):
             
            #set up phasing    
            aa.set_active_pol(pol) 
            if i!=0:
                old_zenith = zenith
            time = t['times'][i]    
            aa.set_jultime(time)
            lst = aa.sidereal_time()
            zenith = aipy.phs.RadioFixedBody(lst, aa.lat)
            zenith.compute(aa)
            if i==0:
                continue

            #phase and multiply, store in vissq_per_antpair
            if i % 2:
                _v1 = totald[i-1,:]
                phase_correction = np.conj(aa.gen_phs(zenith,antpair[0],antpair[1]))*aa.gen_phs(old_zenith,antpair[0],antpair[1])
                _v2 = totald[i,:]*phase_correction[freq_range[0]:freq_range[1]]
                vissq_per_antpair[i // 2,:] = np.conj(_v1)*_v2
    
        antpairslices.append(np.log10(np.fft.fftshift(np.abs(vissq_per_antpair), axes=1)))
        print "antpair {antpair} done!!!"
    
    antpairs = antdict[baselengths[bl_num]]

    channel_width = (t['freqs'][1] - t['freqs'][0])*10**3
    num_bins = len(t['freqs'])
    #get data to recalculate axes  
    delays = np.fft.fftshift(np.fft.fftfreq(num_bins, channel_width / num_bins))

    np.savez(npz_name, antpairslc=antpairslices, dlys=delays, pol=pol, antprs=antpairs, length=length, hist=history)
    return npz_name

# Create Wedge from Pitchfork
def wedge_delayavg(npz_name, multi=False):

    plot_data = np.load(npz_name)
    delays, wedgevalues, baselines = plot_data['dlys'], plot_data['wdgslc'], plot_data['bls']
    d_start = plot_data['dlys'][0]
    d_end = plot_data['dlys'][-1]
    split = (len(wedgevalues[0,:])/2)

    wedgevalues2 = np.zeros((len(wedgevalues),len(delays)))

    for baselength in range(len(wedgevalues)):        
        for i in range(split):
            avg = ((wedgevalues[baselength,(split-1+i)]+wedgevalues[baselength,(split+i)])/2)
            wedgevalues2[baselength][split-i] = avg         
    delayavg_wedgevalues = wedgevalues2.T.T.T
    npz_delayavg = (npz_name[:-11] + 'delayavg.npz')
    np.savez(npz_delayavg, wdgslc=delayavg_wedgevalues, dlys=delays, bls=baselines)
                                #saving to longer arrary fml??? idk
    print "got here!!!"
    return npz_delayavg

# Plotting Functions:
def plot_timeavg(npz_name, multi=False):
    plot_data = np.load(npz_name)

    d_start = plot_data['dlys'][0]
    d_end = plot_data['dlys'][-1]
    plot = plt.imshow(plot_data['wdgslc'], aspect='auto',interpolation='nearest',extent=[d_start,d_end,len(plot_data['wdgslc']),0], vmin=-3.0, vmax=1.0)
    plt.xlabel("Delay (ns)", size='medium')
    plt.ylabel("Baseline length (m)", size='medium')
    cbar = plt.colorbar()
    cbar.set_label("log10((mK)^2)")
    plt.xlim((-450,450))
    npz_name = npz_name.split('/')[-1]
    plt.suptitle('JD '+npz_name.split('.')[1]+' LST '+plot_data['lst'][0]+' to '+plot_data['lst'][-1], size='large')
    plt.title(npz_name.split('.')[3], size='medium')
    plt.yticks(np.arange(len(plot_data['bls'])), [round(n,1) for n in plot_data['bls']])

    #calculate light travel time for each baselength
    light_times = []
    for length in plot_data['bls']:
        light_times.append(length/sc.c*10**9)

    #plot lines on plot using the light travel time
    for i in range(len(light_times)):
       x1, y1 = [light_times[i], light_times[i]], [i, i+1] 
       x2, y2 = [-light_times[i], -light_times[i]], [i, i+1]
       plt.plot(x1, y1, x2, y2, color = 'white')
    
    if multi:
        return
    else:
        plt.savefig(npz_name[:-3] + 'png')
        plt.show()

def plot_multi_timeavg(npz_names):
    #set up multiple plots
    nplots = len(npz_names)
    plt.figure(figsize=(4*nplots-3,3))
    G = gridspec.GridSpec(3, 4*nplots-4)

    #plot each plot in its own gridspec area   
    for i in range(len(npz_names)):
        axes = plt.subplot(G[:, (i*3):(i*3)+3])
        plot_timeavg(npz_names[i], multi=True)

    plt.tight_layout()
    plt.savefig(npz_names[0][:-3] + "multi.png")
    # plt.show()

def plot_blavg(npz_name): 
    plot_data = np.load(npz_name)

    d_start = plot_data['dlys'][0]
    d_end = plot_data['dlys'][-1]
    # t_start = plot_data['wdgslc'][0].shape[0]
    # t_end = 0

    #format time scale in minutes
    t_start = (plot_data['times'][-1] - plot_data['times'][0]) * 24*60
    t_end = 0

    #create subplot to plot data
    f,axarr = plt.subplots(len(plot_data['wdgslc']),1,sharex=True,sharey=True)
    
    #add axes labels
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', 
                    right='off')
    plt.xlabel("Delay (ns)")
    plt.ylabel("Time (min)")
    plt.suptitle('JD '+npz_name.split('.')[1]+' LST '+str(plot_data['lst'][0])+' to '+str(plot_data['lst'][-1]), size='large')
    plt.title(npz_name.split('.')[3], size='medium')

    #calculate light travel time for each baselength
    light_times = []
    for length in plot_data['bls']:
        light_times.append(length/sc.c*10**9)
 
    #plot individual wedge slices
    for i in range(len(plot_data['wdgslc'])):
        #plot the graph
        im = axarr[i].imshow(plot_data['wdgslc'][i], aspect='auto',interpolation='nearest', vmin=-9, vmax= 1, extent=[d_start,d_end,t_start,t_end])
        #plot light delay time lines
        light_time = (plot_data['bls'][i])/sc.c*10**9
        x1, y1 = [light_time, light_time], [t_start, t_end] 
        x2, y2 = [-light_time, -light_time], [t_end, t_start]
        axarr[i].plot(x1, y1, x2, y2, color = 'white') 

    cax,kw = mpl.colorbar.make_axes([ax for ax in axarr.flat])
    cbar = plt.colorbar(im, cax=cax, **kw)
    cbar.set_label('log10(mK)')
    
    #scale x axis to the significant information
    axarr[0].set_xlim(-450,450)

    plt.savefig(npz_name[:-3] + 'png')
    f.set_size_inches(5, 11, forward=True)
    plt.show()

def plot_flavors(npz_name, multi=False):
    npz_name = "".join(npz_name)
    data = np.load(npz_name)

    npz_name = npz_name.split('/')[-1]

    axis_delay_start = data['dlys'][0]
    axis_delay_end = data['dlys'][-1]
    plot = plt.imshow(
        data['wdgslc'],
        aspect='auto',
        # interpolation='nearest',
        extent=[axis_delay_start, axis_delay_end, len(data['wdgslc']), 0],
        vmin=-3.0,
        vmax=1.0)
    
    plt.xlabel("Delay (ns)")

    plt.xlim((-450, 450))

    ticks = []
    slopedict = data['slpdct'].tolist()
    for baseline in sorted(slopedict.keys()):
        for slope in sorted(slopedict[baseline].keys()):
            ticks.append("{:.3}: {:8.3}".format(baseline, slope))

    light_times = []
    for baseline in sorted(slopedict.keys()):
        for slope in sorted(slopedict[baseline].keys()):
            light_times.append(baseline/sc.c*10**9)

    for i in range(len(light_times)):
        x1, y1 = [light_times[i], light_times[i]], [i, i+1] 
        x2, y2 = [-light_times[i], -light_times[i]], [i, i+1]
        plt.plot(x1, y1, 'w', x2, y2, 'w')

    plt.axvline(x=0, color='k', linestyle='dashed', linewidth=0.5)


    if multi:
        plt.title(".".join(npz_name.split('.')[2:4]))
        return ticks
    else:
        plt.title(npz_name.split('.')[3])
        plt.ylabel("Baseline length (short to long)")
        color_bar = plt.colorbar().set_label("log10((mK)^2)")
        plt.yticks(np.arange(len(ticks)), ticks)
        plt.tight_layout()
        plt.savefig(".".join(npz_name.split('.')[:-1]) + '.png')
        plt.show()

def plot_multi_flavors(npz_names):
    nplots = len(npz_names)
    plt.figure(figsize=(15, 5))
    G = gridspec.GridSpec(3, (4 * nplots) - 4)

    for i in range(nplots):
        axes = plt.subplot(G[:, (i*3):(i*3)+3])
        ticks = plot_flavors(npz_names[i], multi=True)
        if i != 0:
            plt.yticks([])
            continue

        plt.yticks(np.arange(len(ticks)), ticks)
        plt.ylabel("Baseline length (short to long)")
    color_bar = plt.colorbar().set_label("log10((mK)^2)")

    # plt.tight_layout()
    npz_names = [file.split('/')[-1] for file in npz_names]
    print npz_names
    plt.savefig(".".join(npz_names[0].split('.')[0:3] + npz_names[0].split('.')[4:8]) + '.multi.png')
    # plt.show()

def plot_bltype(npz_name):

    plot_data = np.load(npz_name)

    d_start = plot_data['dlys'][0]
    d_end = plot_data['dlys'][-1]
    t_start = plot_data['antpairslc'][0].shape[0]

    #create subplot to plot data
    f,axarr = plt.subplots(len(plot_data['antpairslc']),1,sharex=True,sharey=True)
    
    #add axes labels
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', 
                    right='off')
    plt.xlabel("Delay (ns)")
    plt.ylabel("Time", labelpad=15)

    plt.suptitle(npz_name.split('.')[1]+'.'+npz_name.split('.')[2]+'.'+npz_name.split('.')[3]+'.baseline'+npz_name.split('.')[7])
    lengthstr = str(plot_data['length'])    
    plt.title('baseline length:'+lengthstr)

    #plot individual wedge slices
    for i in range(len(plot_data['antpairslc'])):
        #plot the graph
        im = axarr[i].imshow(plot_data['antpairslc'][i], aspect='auto',interpolation='nearest', vmin=-9, vmax= 1, extent=[d_start,d_end,t_start,0])
        #plot light delay time lines
        light_time = (plot_data['length'])/sc.c*10**9
        x1, y1 = [light_time, light_time], [0, np.shape(plot_data['antpairslc'][i])[0]] 
        x2, y2 = [-light_time, -light_time], [0, np.shape(plot_data['antpairslc'][i])[0]]
        axarr[i].plot(x1, y1, x2, y2, color = 'white')
        axarr[i].set_ylabel(plot_data['antprs'][i], fontsize=6) 

    cax,kw = mpl.colorbar.make_axes([ax for ax in axarr.flat])
    plt.colorbar(im, cax=cax, **kw)
    
    #scale x axis to the significant information
    axarr[0].set_xlim(-450,450)
    
    f.set_size_inches(6, 9, forward=True)
    plt.savefig(npz_name[:-3] + 'png')
    plt.show()

def plot_delayavg(npz_delayavg):
    
    plot_data = np.load(npz_delayavg)
    delays, wedgevalues, baselines = plot_data['dlys'], plot_data['wdgslc'], plot_data['bls']
    d_start = plot_data['dlys'][0]
    d_end = plot_data['dlys'][-1]
    plot = plt.imshow(wedgevalues, aspect='auto', interpolation='nearest',extent=[0,len(npz_delayavg),d_start,d_end], vmin=-3.0, vmax=1.0)      
    #plot = plt.imshow(npz_delayavg, aspect='auto', interpolation='nearest',extent=[0,len(wedgevalues),d_start,d_end], vmin=-3.0, vmax=1.0)
   
    plt.xlabel("Baseline length (short to long)")
    plt.ylabel("Delay (ns)")
    cbar = plt.colorbar()
    cbar.set_label("log10((mK)^2)")
    plt.xlim((0,len(baselines)))
    plt.ylim(0,450)
    plt.title(npz_delayavg.split('.')[1]+'.'+npz_delayavg.split('.')[2]+'.'+npz_delayavg.split('.')[3])

    #calculate light travel time for each baselength
    light_times = []
    for length in plot_data['bls']:
        light_times.append(length/sc.c*10**9)

    #plot lines on plot using the light travel time
    for i in range(len(light_times)):
       y1, x1 = [light_times[i], light_times[i]], [i, i+1] 
       y2, x2 = [-light_times[i], -light_times[i]], [i, i+1]
       plt.plot(x1, y1, x2, y2, color = 'white')
    
    print "got here1"   
    plt.savefig(npz_delayavg[:-12]+'delayavg.png')
    plt.show()

# 1D Plotting for timeavg wedges
def plot_1D(npz_name, baselines=[]):

    """
    Plots all baselines overlapped on a 1D plot.
    If baselines is a specified argument (start indexing with baseline length #1),
    then only plots the provided baselines.
    """

    plot_data = np.load(npz_name)

    if len(baselines):
        baselines = [i-1 for i in baselines]
    else:
        baselines = range(len(plot_data['wdgslc']))
    
    plt.figure(figsize=(12,6))
    G = gridspec.GridSpec(2, 9)
    
    axes = plt.subplot(G[:,0:4])
    for i in baselines:
        plt.plot(plot_data['dlys'], plot_data['wdgslc'][i], label='bl len '+str(plot_data['bls'][i]))
    plt.xlabel('Delay (ns)')
    plt.ylabel('log10((mK)^2)')
    plt.legend(loc='upper left')
    plt.ylim((-3.5,2.0)) 
   
    axes = plt.subplot(G[:,5:9])
    for i in baselines:
        plt.plot(plot_data['dlys'], plot_data['wdgslc'][i])
    if len(baselines)==1:
        light_time = plot_data['bls'][baselines[0]]/sc.c*10**9
        plt.axvline(light_time, color='#d3d3d3', linestyle='--')
        plt.axvline(-1*light_time, color='#d3d3d3', linestyle='--')
        plt.axvline(0, color='#d3d3d3', linestyle='--')
    plt.xlim((-450,450))
    plt.ylim((-3.5,2.0))
    
    plt.xlabel('Delay (ns)')
    plt.ylabel('log10((mK)^2)')
    npz_name = npz_name.split('/')[-1]
    plt.suptitle(npz_name.split('.')[1]+'.'+npz_name.split('.')[2]+'.'+npz_name.split('.')[3])
        
    plt.show()
    
def plot_multi_1D(npz_names, baselines=[]):
    """
    Plots four 1D plots next to each other.
    If baselines is a specified argument (start indexing with baseline lengt #1),
    then only plots the the provided baselines.
    """

    plot_data = np.load(npz_names[0])
    #set up baselines
    if len(baselines):
        baselines = [i-1 for i in baselines]
    else:
        baselines = range(len(plot_data['wdgslc']))

    #set up the plotting space
    plt.figure(figsize=(18,4))
    G = gridspec.GridSpec(1,4)
    
    #plot each 1D plot
    polorder = ''
    for n in range(len(npz_names)):
        
        #load data, format plotting section
        plot_data = np.load(npz_names[n])
        axes = plt.subplot(G[:,n:n+1])
        
        #plot the data
        for i in baselines:
            plt.plot(plot_data['dlys'], plot_data['wdgslc'][i], label='bl len '+str(plot_data['bls'][i]))
        if len(baselines)==1:
            light_time = plot_data['bls'][baselines[0]]/sc.c*10**9
            plt.axvline(light_time, color='#d3d3d3', linestyle='--')
            plt.axvline(-1*light_time, color='#d3d3d3', linestyle='--')
        plt.axvline(0, color='#d3d3d3', linestyle='--')
        plt.xlim((-450,450))
        plt.ylim((-3.0,2.0))
        
        if n==0:
            plt.legend(loc='upper left')
        plt.xlabel('Delay (ns)')
        plt.ylabel('log10((mK)^2)')
        pol = npz_names[n].split('/')[-1].split('.')[3]
        plt.title(pol)
        polorder += pol
    
    npz_name = npz_names[0].split('/')[-1]
    plt.suptitle(npz_name.split('.')[1]+'.'+npz_name.split('.')[2])
    
    if len(baselines) == len(plot_data['wdgslc']):
        blstr = 'allbls'
    else:
        blstr = 'bl'
        for bl in baselines: blstr += str(bl+1)
    
    plt.tight_layout()
    plt.savefig(npz_name.split(polorder[0])[0]+polorder+npz_name.split(polorder[0])[-1][:-3] + "multi1D." + blstr + ".png")
    plt.show()

# Inside/Outside Averaging
def in_out_avg(npz_name):
    data = np.load(npz_name)
    history = data['hist'].tolist()
    num_files = len(history['filenames'])

    light_times = []
    for length in data['bls']:
        light_times.append(length / (sc.c * (10**9)))

    total_in, total_out = 0, 0
    total_in_count, total_out_count = 0, 0
    for i in range(len(data['bls'])):
        for index, delay in enumerate(data['dlys']):
            if abs(delay) >= light_times[i]:
                total_out += data['wdgslc'][i][index]
                total_out_count += 1
            else:
                total_in += data['wdgslc'][i][index]
                total_in_count += 1

    avg_in = total_in / total_in_count
    avg_out = total_out / total_out_count

    return (avg_in, avg_out, num_files)

def plot_avgs(npz_names, rng):
    name = ".".join(npz_names[-1].split('/')[-1].split('.')[:7])

    total_files = []
    avgs_in = []
    avgs_out = []
    for npz_name in npz_names:
        avgs_in_out = in_out_avg(npz_name)
        total_files.append(avgs_in_out[2])
        avgs_in.append(avgs_in_out[0])
        avgs_out.append(avgs_in_out[1])

    plt.figure(figsize=(15, 10))

    plot_avgs_out = plt.scatter(total_files, avgs_out)
    plot_avgs_in = plt.scatter(total_files, avgs_in)
    plt.legend((plot_avgs_out, plot_avgs_in), ('Averages Outside Wedge', 'Averages Inside Wedge'))
  
    plt.xlim(0, len(total_files))
    plt.xticks(np.arange(0, len(total_files)+4, 2))
    plt.xlabel('Number of Files Used in Averaging')

    plt.ylim(-3.5, 1.5)
    plt.yticks(np.arange(-3, 2))
    plt.ylabel('log10((mK)^2)')

    plt.title(name)
    plt.suptitle('How Well does the Pitchfork Average Down?')

    # Setting range to calculate a fit over.
    # rng = (7, 26)
    x = np.arange(rng[0], rng[1])

    # Outside RMSE and fit calculation.
    m, b = np.polyfit(x, avgs_out[rng[0]:rng[1]], 1)
    fit = np.poly1d([m, b])
    outside = avgs_out[rng[0]:rng[1]]

    RMSE = gu.RMSE(fit(x), outside)
    plt.plot(x, m*x + b, '--', color='k')
    plt.text(8, 0, 'Outside Pitchfork Fit RMSE: {:.3}\nOutside Pitchfork Fit Slope: {:.3}'.format(RMSE, m))


    # Inside RMSE and fit calculation.
    m, b = np.polyfit(x, avgs_in[rng[0]:rng[1]], 1)
    fit = np.poly1d([m, b])
    inside = avgs_in[rng[0]:rng[1]]

    RMSE = gu.RMSE(fit(x), inside)
    plt.plot(x, m*x + b, '--', color='k')
    plt.text(8, -3, 'Inside Pitchfork Fit RMSE: {:.3}\nInside Pitchfork Fit Slope: {:.3}'.format(RMSE, m))

    plt.savefig(name + ".avg_val.png")
    plt.show()