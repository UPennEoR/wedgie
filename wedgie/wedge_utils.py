"""
Module for wedge-creation methods
"""
import capo, aipy, os, pprint
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.constants as sc
import gen_utils as gu
import cosmo_utils as cu

def in_out_avg(npz_name):
    plot_data = np.load(npz_name)

    light_times = []
    for length in plot_data['bls']:
        light_times.append(length / (sc.c * (10**9)))

    total_in, total_out = 0, 0
    total_in_count, total_out_count = 0, 0
    for i in range(len(plot_data['bls'])):
        for index, delay in enumerate(plot_data['dlys']):
            if abs(delay) >= light_times[i]:
                total_out += plot_data['wdgslc'][i][index]
                total_out_count += 1
            else:
                total_in += plot_data['wdgslc'][i][index]
                total_in_count += 1

    avg_in = total_in / total_in_count
    avg_out = total_out / total_out_count
    num_files = len(plot_data['hist'].item()['filenames'][0])

    return (avg_in, avg_out, num_files)

def plot_avgs(npz_names):
    total_files = []
    avgs_in = []
    avgs_out = []
    for npz_name in npz_names:
        avgs_in_out = in_out_avg(npz_name)
        total_files.append(avgs_in_out[2])
        avgs_in.append(avgs_in_out[0])
        avgs_out.append(avgs_in_out[1])

    pprint.pprint(total_files)
    pprint.pprint(avgs_in)
    pprint.pprint(avgs_out)
    plot_avgs_out = plt.scatter(total_files, avgs_out)
    plot_avgs_in = plt.scatter(total_files, avgs_in)

    plt.legend((plot_avgs_out, plot_avgs_in), ('Averages Outside Wedge', 'Averages Inside Wedge'))

    plt.xlim((0, 20))
    plt.ylim(-3.5, 1.5)
    # plt.savefig('fig1.png')
    plt.show()


# Calfile-specific manipulations

def calculate_baseline(antennae, pair):
    """
    XXX This problem has been "solved" with numpy instead.

    The decimal module is necessary for keeping the number of decimal places small.
    Due to small imprecision, if more than 8 or 9 decimal places are used, 
    many baselines will be calculated that are within ~1 nanometer to ~1 picometer of each other.
    Because HERA's position precision is down to the centimeter, there is no 
    need to worry about smaller imprecision.
    """

    dx = antennae[pair[0]]['top_x'] - antennae[pair[1]]['top_x']
    dy = antennae[pair[0]]['top_y'] - antennae[pair[1]]['top_y']
    baseline = np.around([np.sqrt(dx**2. + dy**2.)],3)[0] #XXX this may need tuning
    
    return baseline

def get_baselines(calfile, ex_ants=[]):
    """
    Returns a dictionary of redundant baselines based on a calfile, and
    excluding bad antennae, specified by ex_ants, a list of integers.
    
    Requires cal file to be in PYTHONPATH.
    """
    try:
        print 'reading, %s'%calfile
        exec("import {cfile} as cal".format(cfile=calfile))
        antennae = cal.prms['antpos_ideal']
    except ImportError:
        raise Exception("Unable to import {cfile}.".format(cfile=calfile))
    
    """
    determines the baseline and places them in the dictionary.
    excludes antennae with z-position < 0 or if in ex_ants list
    """
    
    baselines = {}
    for antenna_i in antennae:
        if antennae[antenna_i]['top_z'] < 0.:
            continue
        if antenna_i in ex_ants:
            continue
        
        for antenna_j in antennae:
            if antennae[antenna_j]['top_z'] < 0.:
                continue
            if antenna_j in ex_ants:
                continue

            if antenna_i == antenna_j:
                continue
            elif antenna_i < antenna_j:
                pair = (antenna_i, antenna_j)
            elif antenna_i > antenna_j:
                pair = (antenna_j, antenna_i)

            baseline = calculate_baseline(antennae, pair)

            if (baseline not in baselines):
                baselines[baseline] = [pair]
            elif (pair in baselines[baseline]):
                continue
            else:
                baselines[baseline].append(pair)
    # pprint.pprint(baselines)
    return baselines

# wedge/pitchfork calculation methods
def wedge_blavg(filenames, pol, calfile, history, ex_ants=[]):
    """
    Plots wedges per baseline length, averaged over baselines.
    Remember to not include the ".py" in the name of the calfile
    """
    
    #get data from file
    t,d,f = capo.miriad.read_files(filenames,antstr='cross', polstr=pol)
    
    #create variable to store the wedge subplots in
    wedgeslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = get_baselines(calfile, ex_ants)
    baselengths = antdict.keys()
    baselengths.sort()
    
    #for each baselength in the dictionary
    for length in baselengths:
        
        totald = np.zeros_like(d[antdict[length][0]][pol]) #variable to store cumulative fft data
        
        #cycle through every ant pair of the given baselength
        for antpair in antdict[length]:

            #CLEAN and fft the data
            clean=1e-3
            w = aipy.dsp.gen_window(d[antpair][pol].shape[-1], window='blackman-harris')
            _dw = np.fft.ifft(d[antpair][pol]*w)
            _ker= np.fft.ifft(f[antpair][pol]*w)
            gain = aipy.img.beam_gain(_ker)
            for time in range(_dw.shape[0]):
                _dw[time,:],info = aipy.deconv.clean(_dw[time,:], _ker[time,:], tol=clean)
                _dw[time,:] += info['res']/gain
          
            totald += np.ma.array(_dw)
        
        #get average of all values for this baselength, store in wedgeslices
        totald /= len(antdict[length])
        wedgeslices.append(np.log10(np.fft.fftshift(np.abs(totald), axes=1)))

    #get data to recalculate axes  
    delays = np.fft.fftshift(np.fft.fftfreq(1024, .1/1024)) #XXX hard coded #1024 bins, channel width of 0.1 GHz/1024

    #save filedata as npz
    #NB: filename of form like "zen.2457746.16693.xx.HH.uvcOR"
    fn1, fn2 = filenames[0].split('.'), filenames[-1].split('.')
    npz_name = fn1[0]+'.'+fn1[1]+'.'+fn1[2]+'_'+fn2[2]+'.'+fn1[3]+'.'+fn1[4]+'.'+fn1[5]+'.blavg.npz'
    np.savez(npz_name, wdgslc=wedgeslices, dlys=delays, pol=pol, bls=baselengths, hist=history)
    return npz_name

def wedge_timeavg(filenames, pol, calfile, history, ex_ants=[], stokes=[]):
    """
    Plots wedges per baseline length, averaged over baselines and time
    if stokes is specified, then it should be of form [t, d, f]
    """
    if len(stokes):
        #label given data if provided
        t,d,f = stokes[0], stokes[1], stokes[2]
    else: 
        #get data from file
        t,d,f = capo.miriad.read_files(filenames,antstr='cross',polstr=pol) 
    
    #stores vis^2 for each baselength averaged over time
    wedgeslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = get_baselines(calfile, ex_ants)
    baselengths = antdict.keys()
    baselengths.sort()

    #get number of times and number of channels
    ntimes,nchan = len(t['times']),len(t['freqs'])
    dt = np.diff(t['times'])[0] #dJD

    #get vis^2 for each baselength
    for baselength in baselengths:

        vissq_per_bl = np.zeros((ntimes // 2 ,nchan))
        
        #go through each individual antenna pair
        for antpair in antdict[baselength]:

            #create/get metadata    
            uv = aipy.miriad.UV(filenames[0])
            aa = aipy.cal.get_aa(calfile.split('.')[0], uv['sdf'], uv['sfreq'], uv['nchan']) 
            del(uv)

            #CLEAN and fft the data
            clean=1e-3
            w = aipy.dsp.gen_window(d[antpair][pol].shape[-1], window='blackman-harris')
            _dw = np.fft.ifft(d[antpair][pol]*w)
            _ker= np.fft.ifft(f[antpair][pol]*w)
            gain = aipy.img.beam_gain(_ker)
            for time in range(_dw.shape[0]):
                _dw[time,:],info = aipy.deconv.clean(_dw[time,:], _ker[time,:], tol=clean)
                _dw[time,:] += info['res']/gain
            ftd_2D_data = np.ma.array(_dw)

            #holds our data
            vissq_per_antpair = np.zeros((ntimes // 2 ,nchan))

            #multiply at times (1*2, 3*4, etc...) 
            for i in range(ntimes):
                 
                #set up phasing    
                aa.set_active_pol(pol) 
                if i!=0:
                    old_zenith = zenith
                else:    
                    time = t['times'][i]    
                aa.set_jultime(time)
                lst = aa.sidereal_time()
                zenith = aipy.phs.RadioFixedBody(lst, aa.lat)
                zenith.compute(aa)
                if i==0:
                    continue

                #phase and multiply, store in vissq_per_antpair
                if i % 2:
                    _v1 = ftd_2D_data[i-1,:]
                    phase_correction = np.conj(aa.gen_phs(zenith,antpair[0],antpair[1]))*aa.gen_phs(old_zenith,antpair[0],antpair[1])
                    _v2 = ftd_2D_data[i,:]*phase_correction
                    vissq_per_antpair[i // 2,:] = np.conj(_v1)*_v2

            #store time average for this baseline length
            vissq_per_bl += vissq_per_antpair

        #compute average for baseline length, average over time, and store in wedgeslices
        vissq_per_bl /= len(antdict[baselength])
        wedgeslices.append(np.log10(np.fft.fftshift(np.mean(np.abs(vissq_per_bl), axis=0))))
        print 'finished a wedgeslice!'

    #get delays
    delays = np.fft.fftshift(np.fft.fftfreq(1024, .1/1024)) #XXX hardcoded #1024 bins, channel width of 0.1 GHz/1024 
    
    #save filedata as npz
    #NB: filename of form like "zen.2457746.16693.xx.HH.uvcOR"
    fn1, fn2 = (filenames[0].split('/')[-1]).split('.'), (filenames[-1].split('/')[-1]).split('.')
    npz_name = fn1[0]+'.'+fn1[1]+'.'+fn1[2]+'_'+fn2[2]+'.'+pol+'.'+fn1[4]+'.'+fn1[5]+'.timeavg.npz'
    np.savez(npz_name, wdgslc=wedgeslices, dlys=delays, pol=pol, bls=baselengths, hist=history)
    return npz_name

def wedge_stokes(filenames, calfile, history, ex_ants=[]):
    """
    calls wedge_timeavg for each stokes parameter
    assumes filenames is a list of lists separated by pol:
        [[xx files],[xy files],[yx files],[yy files]]
    """
    
    txx,dxx,fxx = capo.miriad.read_files(filenames[0],antstr='cross',polstr='xx')
    # txy,dxy,fxy = capo.miriad.read_files(filenames[1],antstr='cross',polstr='xy')
    # tyx,dyx,fyx = capo.miriad.read_files(filenames[2],antstr='cross',polstr='yx')
    tyy,dyy,fyy = capo.miriad.read_files(filenames[3],antstr='cross',polstr='yy')
    
    #calculate I (VI = Vxx + Vyy)
    tI = txx
    dI = {}
    fI = {}
    for key in dxx.keys():
        dI[key] = {'I': dxx[key]['xx'] + dyy[key]['yy'] }
        fI[key] = {'I': fxx[key]['xx'] + fyy[key]['yy'] }
    nameI = wedge_timeavg(filenames[0], 'I', calfile, history, ex_ants, stokes=[tI, dI, fI])
    print 'Stokes I completed'

    """
    #calculate Q (VQ = Vxx - Vyy)
    tQ = txx
    dQ = {}
    fQ = {}
    for key in dxx.keys():
        dQ[key] = {'Q': dxx[key]['xx'] - dyy[key]['yy'] }
        fQ[key] = {'Q': fxx[key]['xx'] + fyy[key]['yy'] }
    nameQ = wedge_timeavg(filenames[0], 'Q', calfile, history, ex_ants, stokes=[tQ, dQ, fQ])
    print 'Stokes Q completed'
    
    #calculate U (VU = Vxy + Vyx)
    tU = txy
    dU = {}
    fU = {}
    for key in dxy.keys():
        dU[key] = {'U': dxy[key]['xy'] + dyx[key]['yx'] }
        fU[key] = {'U': fxy[key]['xy'] + fyx[key]['yx'] }
    nameU = wedge_timeavg(filenames[2], 'U', calfile, history, ex_ants, stokes=[tU, dU, fU])
    print 'Stokes U completed'
    
    #calculate V (VV = -i*Vxy + i*Vyx)
    tV = txy
    dV = {}
    fV = {}
    for key in dxy.keys():
        dV[key] = {'V': -1j*dxy[key]['xy'] + 1j*dyx[key]['yx'] }
        fV[key] = {'V': fxy[key]['xy'] + fyx[key]['yx'] }
    nameV = wedge_timeavg(filenames[2], 'V', calfile, history, ex_ants, stokes=[tV, dV, fV])
    print 'Stokes V completed'
    """
    return [nameI]
    # , nameQ, nameU, nameV]

# Plotting Routines
def plot_blavg(npz_name, path='./'): 
    plot_data = np.load(npz_name)

    d_start = plot_data['dlys'][0]
    d_end = plot_data['dlys'][-1]
    t_start = plot_data['wdgslc'][0].shape[0]

    #create subplot to plot data
    f,axarr = plt.subplots(len(plot_data['wdgslc']),1,sharex=True,sharey=True)
    
    #add axes labels
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', 
                    right='off')
    plt.xlabel("Delay (ns)")
    plt.ylabel("Time")

    #calculate light travel time for each baselength
 
    #plot individual wedge slices
    for i in range(len(plot_data['wdgslc'])):

        #plot the graph
        axarr[i].imshow(plot_data['wdgslc'][i], aspect='auto',interpolation='nearest',vmin=-6,vmax=-1, extent=[d_start,d_end,t_start,0])

        #plot light delay time lines
        light_time = (plot_data['bls'][i])/sc.c*10**9
        x1, y1 = [light_time, light_time], [0, np.shape(plot_data['wdgslc'][i])[0]] 
        x2, y2 = [-light_time, -light_time], [0, np.shape(plot_data['wdgslc'][i])[0]]
        axarr[i].plot(x1, y1, x2, y2, color = 'white') 
     
    #scale x axis to the significant information
    axarr[0].set_xlim(-450,450)

    plt.savefig(path + npz_name[:-3] + 'png')
    plt.show()

def plot_timeavg(npz_name, path='./', multi=False):
    plot_data = np.load(npz_name)
    
    d_start = plot_data['dlys'][0]
    d_end = plot_data['dlys'][-1]
    plot = plt.imshow(plot_data['wdgslc'], aspect='auto',interpolation='nearest',extent=[d_start,d_end,len(plot_data['wdgslc']),0], vmin=-3.0, vmax=1.0)
    plt.xlabel("Delay (ns)")
    plt.ylabel("Baseline length (short to long)")
    cbar = plt.colorbar()
    cbar.set_label("log10((mK)^2)")
    plt.xlim((-450,450))
    plt.title(npz_name.split('.')[1]+'.'+npz_name.split('.')[2]+'.'+npz_name.split('.')[3])

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
        plt.savefig(path + npz_name[:-3] + 'png')
        plt.show()

def plot_multi_timeavg(npz_names, path='./'):
    #set up multiple plots
    nplots = len(npz_names)
    plt.figure(figsize=(4*nplots-3,3))
    G = gridspec.GridSpec(3, 4*nplots-4)

    #plot each plot in its own gridspec area   
    for i in range(len(npz_names)):
        axes = plt.subplot(G[:, (i*3):(i*3)+3])
        plot_timeavg(npz_names[i], path, multi=True)

    plt.tight_layout()
    plt.savefig(path + npz_names[0][:-3] + "multi.png")
    plt.show()

def wedge_delayavg(npz_name):

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

def plot_multi_delayavg(npz_names):
    #set up multiple plots
    nplots = len(npz_names)
    plt.figure(figsize=(4*nplots-3,3))
    G = gridspec.GridSpec(3, 4*nplots-4)

    #plot each plot in its own gridspec area   
    for i in range(len(npz_names)):
        axes = plt.subplot(G[:, (i*3):(i*3)+3])
        wedge_delayavg(npz_names[i], multi=True)

    plt.tight_layout()
    plt.show()  

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
    
    axes = plt.subplot(G[:,5:9])
    for i in baselines:
        plt.plot(plot_data['dlys'], plot_data['wdgslc'][i])
    if len(baselines)==1:
	light_time = plot_data['bls'][baselines[0]]/sc.c*10**9
	plt.axvline(light_time, color='#d3d3d3', linestyle='--')
	plt.axvline(-1*light_time, color='#d3d3d3', linestyle='--')
	plt.axvline(0, color='#d3d3d3', linestyle='--')
    plt.xlim((-450,450))
    
    plt.xlabel('Delay (ns)')
    plt.ylabel('log10((mK)^2)')
    plt.suptitle(npz_name.split('.')[1]+'.'+npz_name.split('.')[2]+'.'+npz_name.split('.')[3])
        
    plt.show()
