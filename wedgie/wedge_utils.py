"""
Module for wedge-creation methods
"""
import baseline_lengths #XXX this should exists here, not in another library
import capo
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import aipy

def plot_wedge_blavg(filenames , pol):
    """
    Plots wedges per baseline length, averaged over baselines
    """
    
    #get data from file
    t,d,f = capo.miriad.read_files(filenames,antstr='cross', polstr=pol)
    
    #create variable to store the wedge subplots in
    wedgeslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = baseline_lengths.get_baselines()
    baselengths = antdict.keys()
    baselengths.sort()
    
    #for each baselength in the dictionary
    for length in baselengths:
        
        totald = np.zeros_like(d[antdict[length][0]][pol]) #variable to store cumulative fft data
        
        #cycle through every ant pair of the given baselength
        for antpair in antdict[length]:
            
            #fft the data wrt freq and sum it with the previous iterations
            ftd_2D_data = np.fft.ifft(d[antpair][pol],axis=1) 
            totald += ftd_2D_data 
        
        #get average of all values for this baselength, store in wedgeslices
        totald /= len(antdict[length])
        wedgeslices.append(totald)

	#get data to recalculate axes  
	delays = np.fft.fftshift(np.fft.fftfreq(1024, .1/1024)) #XXX hard coded #1024 bins, channel width of 0.1 GHz/1024
	d_start = delays[0]
	d_end = delays[-1]
	t_start = wedgeslices[0].shape[0]

    #create subplot to plot data
    f,axarr = plt.subplots(len(wedgeslices),1,sharex=True,sharey=True)
    
    #add axes labels
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', 
                    right='off')
    plt.xlabel("Delay (ns)")
    plt.ylabel("Time")
    
    #plot individual wedge slices
    for i in range(len(wedgeslices)):
        axarr[i].imshow(np.log10(np.fft.fftshift(np.abs(wedgeslices[i]),
             axes=1)),aspect='auto',interpolation='nearest',vmin=-6,vmax=-1, extent=[d_start,d_end,t_start,0])
				    
    #scale x axis to the significant information
    axarr[0].set_xlim(-450,450)
    
    plt.show()

def plot_wedge_timeavg(filenames, calfile, pol):

    """
    Plots wedges per baseline length, averaged over baselines and time
    """
    #get data from file
    t,d,f = capo.miriad.read_files(filenames,antstr='cross',polstr=pol) 
    
    #stores vis^2 for each baselength averaged over time
    wedgeslices = []
    
    #get dictionary of antennae pairs
    #keys are baseline lengths, values are list of tuples (antenna numbers)
    antdict = baseline_lengths.get_baselines()
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
            aa = aipy.cal.get_aa(calfile, uv['sdf'], uv['sfreq'], uv['nchan']) 
            del(uv)

            #fourier transform the data
            ftd_2D_data = np.fft.ifft(d[antpair][pol],axis=1) 

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
        
    
	#plot wedge
    delays = np.fft.fftshift(np.fft.fftfreq(1024, .1/1024)) #XXX hardcoded #1024 bins, channel width of 0.1 GHz/1024 
    d_start = delays[0]
    d_end = delays[-1]
    plot = plt.imshow(wedgeslices, aspect='auto',interpolation='nearest',extent=[d_start,d_end,len(wedgeslices),0])
    plt.xlabel("Delay (ns)")
    plt.ylabel("Baseline length (shortest to longest)")
    cbar = plt.colorbar()
    cbar.set_label("log10((mK)^2)")
    plt.xlim((-450,450))

	#calculate light travel time for each baselength
    light_times = []
    for length in baselengths:
        light_times.append(length/sc.c*10**9)

	#plot lines on plot using the light travel time
    for i in range(len(light_times)):
       x1, y1 = [light_times[i], light_times[i]], [i, i+1] 
       x2, y2 = [-light_times[i], -light_times[i]], [i, i+1]
       plt.plot(x1, y1, x2, y2, color = 'white')

    plt.show()
