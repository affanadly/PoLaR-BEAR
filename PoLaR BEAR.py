#PoLaR BEAR by Affan
progress = True
if progress: print('PoLaR BEAR - by Affan\n')

#initialization
if progress: print('Initializing libraries and parameters,')
import numpy as np
import struct
from bitarray import bitarray
import matplotlib.pyplot as plt
from tabulate import tabulate

input_file = 'bj1.fil'
output_directory = ''
ifchannel = 0 #if channel

start = 0 # % for start time
end = 0.2 # % for end time

ds = 20 #downsampling coefficient
#frequency channels to zap [[fstart1,fend1],[fstart2,fend2],...]
zaplist = [[1435,1440],[1495,1505]]

dmstart, dmend, dminc = 1, 2000, 1 #dedispersion dms

snrloss = 0.15 #maximum snr loss
minw = 0.2e-3 #minimum width
nbox = 12 #number of boxes

maxpeak = 100 #maximum number of peaks
prob = 1e-7 #probability for threshold

waterfall = -1 #waterfall candidate plot to output (-1 for all)
dpi = 300 #resolution of plots and outputs

#header dict initialization
header_string = {
    'source_name': '',
    'rawdatafile': ''
}
header_int = {
    'machine_id': '',
    'telescope_id': '',
    'data_type': '',
    'nchans': '',
    'nifs': '',
    'nbits': '',
    'barycentric': '',
    'pulsarcentric': '',
    'nsamples': '',
    'nbeams': '',
    'ibeam': '',
    'start_sample': '',
    'end_sample': '',
    'cand_location': '',
    'MJD_start': '',
    'MJD_hour': '',
    'MJD_minute': ''
}
header_float = {
    'fch1': '',
    'foff': '',
    'tstart': '',
    'tsamp': '',
    'az_start': '',
    'za_start': '',
    'src_raj': '',
    'src_dej': '',
    'MJD_second': ''
}

#functions for reading header
def read_string(length,byte):
    string = struct.unpack('{:d}s'.format(length),byte)
    string = string[0].decode('utf-8')
    return string
def read_int(byte):
    integer = struct.unpack('<i',byte)
    return integer[0]
def read_float(byte):
    floating = struct.unpack('<d',byte)
    return floating[0]
    
#initialize file reading
if progress: print('Reading filterbank,')
f = open(input_file,'rb')
i = 4
data_type = 0 
#0-length_int 1-header_string 2-header_string_length 3-string 4-int 5-float
header_done = False
data = b''

#header/file reading
while(byte := f.read(i)):
    if header_done:
        data += byte
    elif data_type == 0:
        leng = read_int(byte)
        data_type = 1
        i = leng
    elif data_type == 1:
        temp_header = read_string(leng,byte)
        if temp_header in header_string:
            data_type = 2
            i = 4
        elif temp_header in header_int:
            data_type = 4
            i = 4
        elif temp_header in header_float:
            data_type = 5
            i = 8
        elif temp_header == 'HEADER_END':
            header_done = True
            i = None
        else:
            i = 4
            data_type = 0
    elif data_type == 2:
        leng = read_int(byte)
        if leng != 0:
            data_type = 3
            i = leng
        else:
            data_type = 0
            i = 4
    elif data_type == 3:
        string = read_string(leng,byte)
        header_string[temp_header] = string
        i = 4
        data_type = 0
    elif data_type == 4:
        integer = read_int(byte)
        header_int[temp_header] = integer
        i = 4
        data_type = 0
    elif data_type == 5:
        floating = read_float(byte)
        header_float[temp_header] = floating
        i = 4
        data_type = 0
f.close()

#header dictionary
header = dict(header_string, **header_int, **header_float)
header_string = None
header_int = None
header_float = None

#conversion to bitarray
data_bit = bitarray('')
data_bit.frombytes(data)
data = None

scale = header['nbits']*header['nchans']*header['nifs']
n = np.floor(len(data_bit)/scale)

#slicing bitarray
starti = int(np.floor(start*n)*scale)
endi = int(np.ceil(end*n)*scale)

data_bit = data_bit[starti:endi]

#exporting data from bitarray
if progress: print('Converting data to numpy array,')
if header['nbits'] == 32:
    data = np.frombuffer(data_bit,dtype=np.float32)
elif header['nbits'] == 16:
    data = np.frombuffer(data_bit,dtype=np.uint16)
elif header['nbits'] == 8:
    data = np.frombuffer(data_bit,dtype=np.uint8)
elif header['nbits'] == 4:
    data = np.frombuffer(data_bit,dtype=np.uint8)
    data = (np.unpackbits(data,bitorder='big')).reshape((-1,4))
    data[:,0] = data[:,0]*8
    data[:,1] = data[:,1]*4
    data[:,2] = data[:,2]*2
    data = np.flip((data.sum(axis=1)).reshape((-1,2)),axis=1).flatten()
elif header['nbits'] == 2:
    data = np.frombuffer(data_bit,dtype=np.uint8)
    data = (np.unpackbits(data,bitorder='big')).reshape((-1,2))
    data[:,0] = data[:,0]*2
    data = np.flip((data.sum(axis=1)).reshape((-1,4)),axis=1).flatten()
elif header['nbits'] == 1:
    data = np.unpackbits(np.frombuffer(data_bit,dtype=np.uint8),bitorder='little')
data_bit = None

data = np.reshape(data,(-1,header['nifs'],header['nchans']))
data = data[:,ifchannel,:]
data = np.copy(data).astype(np.float32)

#analysis functions
def downsample(data,ds):
    temp = np.copy(data)
    dscal = int(np.floor(np.shape(temp)[0]/ds))*ds
    return np.add.reduceat(temp[0:dscal,:],
                           np.arange(0,dscal,ds))
                           
def zap(data,zaplist,f):
    temp = np.copy(data)
    for zapf in zaplist:
        for zapch in np.where(np.logical_and(f>=zapf[0],f<=zapf[1])):
            temp[:,zapch] = 0 
    return temp

def zdmf(data):
    temp = np.copy(data)
    sdm0 = np.sum(temp,axis=1)
    sum_sdm0 = np.sum(sdm0)
    N, nchan = np.shape(temp)
    deno = np.inner(sdm0,sdm0) - sum_sdm0**2/N
    for k in range(nchan):
        sdm0si = np.inner(sdm0,temp[:,k])
        sum_si = np.sum(temp[:,k])
        temp[:,k] -= ((sdm0si - sum_si*sum_sdm0/N)/deno)*sdm0
    return np.nan_to_num(temp)

def rembaseline(data):
    temp = np.copy(data)
    basesig = np.mean(temp,axis=0)
    basesig = np.nan_to_num(basesig)
    temp -= basesig
    return temp

def dmdelay(dm,fx,fstart):
    return 4148.741601*(1.0/(fstart**2)
                          -(1.0/(fx**2)))*dm

def dedisperse(data,dm,f,fstart,tsamp,sep=False):
    temp = np.copy(data)
    delay = np.round(dmdelay(dm,f,fstart)/tsamp)
    for i in range(np.shape(temp)[1]):
        temp[:,i] = np.roll(temp[:,i],int(delay[i]))
    return (temp if sep else np.sum(temp,axis=1))
    
def lsr(timeseries,box):
    stddev = np.std(timeseries)
    mask = np.ones((1,box))
    mask = mask[0,:]
    temp = np.convolve(timeseries,mask,'same')
    return temp**2/(box*(stddev**2))

def threshold(dof,prob):
    x = 1
    for i in range(40):
        x = 2*np.log(dof)-0.4516-np.log(x)-2*np.log(prob)
    return x

def peaksearch(cube, thres, maxpeak):
    ddm, dw, dt = 10, 4, 10
    ndm, nw, nt = np.shape(cube)
    mask = np.zeros((ndm,nw,nt),dtype=bool)
    rms = np.std(cube)
    peaks = []
    for i in range(maxpeak):
        peak = np.unravel_index(np.argmax(cube[~mask]),(ndm,nw,nt))
        if cube[peak] < thres: break
        print('Peaks found: {:d}'.format(i+1),end='\r')
        peaks.append(peak+(cube[peak],))
        mask[peak] = True
        neib = [peak]
        for j,k,l in neib:
            m = np.s_[int(0   if j-ddm   <0   else j-ddm):
                      int(ndm if j+ddm+2 >ndm else j+ddm+2),
                      int(0   if k-dw    <0   else k-dw):
                      int(nw  if k+dw+2  >nw  else k+dw+2),
                      int(0   if l-dt    <0   else l-dt):
                      int(nt  if l+dt+2  >nt  else l+dt+2)]
            neib_cube = cube[m]
            cond = (neib_cube > thres) & (neib_cube < cube[peak]
                                          +0.3*rms) & (mask[m] == False)
            np.place(mask[m],cond,True)
            neib.extend(np.add(np.transpose(np.where(cond)),
                        [m[0].start,m[1].start,m[2].start]).tolist())
    return peaks

def freqdist(data,dm,width,time,freq_a,tsamp):
    temp = np.transpose(np.copy(data))
    mean = np.mean(temp,axis=1)
    ileft = np.round((time-width/2+dmdelay(dm,freq_a,
                                        freq1))/tsamp).astype(int)
    iright = np.round((time+width/2+dmdelay(dm,freq_a,
                                        freq1))/tsamp).astype(int)
    zero = np.zeros(np.shape(temp))
    r = np.arange(np.shape(temp)[1])
    imask = (ileft[:,None] <= r) & (iright[:,None] >= r)
    zero[:len(imask)][imask] = temp[:len(imask)][imask]
    return np.sum((zero.T-mean).T,axis=1)**2/np.round(width/tsamp)

def waterfall_plot(n): #######
    fig, axs = plt.subplots(5,6,figsize=(16,10),dpi=dpi,gridspec_kw={
                'height_ratios':[2,4,1,4,1],'width_ratios':[1,1,6,1,1,4],
                'hspace':0,'wspace':0}) 
    fig.subplots_adjust(top=0.92)
    fig.suptitle('Candidate {:d}, DM = {:.2f} cm⁻³ pc'.format(
                    n+1,dm_a[peaks[n][0]])
                 +'\nSource: '+str(header['source_name'])
                 +', Time (MJD): {:.15f}'.format(header['tstart']))
    
    deddata = dedisperse(data,dm_a[peaks[n][0]],
                          freq_a,freq1,tsamp,True)
    dw = 20 
    l = [int(np.round((peaks[n][2]*tsamp-(2*dw+1)*w_a[peaks[n][1]]*
                tsamp/2)/tsamp)),int(np.round((peaks[n][2]*
                tsamp+(2*dw+1)*w_a[peaks[n][1]]*tsamp/2)/tsamp))+1]
    if l[0]<0: l[0] = 0
    if l[1]>len(t_a): l[1] = len(t_a)

    for blank in [(0,0),(0,1),(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),
                  (2,0),(2,1),(2,2),(2,3),(2,4),(2,5),(3,0),(3,3),
                  (4,0),(4,1),(4,3),(4,4)]:
        axs[blank].set_visible(False)
    
    #pulse profile
    axs[0,2].get_shared_x_axes().join(axs[0,2],axs[1,2],axs[3,2],axs[4,2])
    axs[0,2].plot(t_a[l[0]:l[1]],np.sum(deddata[l[0]:l[1],:],
                                 axis=1),'-k',linewidth=0.5)
    axs[0,2].scatter(t_a[peaks[n][2]],0,
                     s=800,facecolors='none',edgecolors='r')
    axs[0,2].set(ylabel='Flux (aU)')
    axs[0,2].set_xticklabels([])
    axs[0,2].tick_params(axis='x',direction='inout')
    axs[0,2].tick_params(axis='y',direction='in')
    
    dist = freqdist(data,dm_a[peaks[n][0]],w_a[peaks[n][1]]*tsamp,
                peaks[n][2]*tsamp,freq_a,tsamp)
    
    #frequency distribution
    axs[1,0].set(ylabel='Frequency (MHz)')
    axs[1,0].get_shared_y_axes().join(axs[1,0], axs[1,1],axs[1,2])
    axs[1,0].plot(np.cumsum(dist),freq_a,'-k',lw=0.5)
    axs[1,0].invert_xaxis()
    axs[1,0].set(xlabel='F (aU)')
    axs[1,0].tick_params(axis='y',direction='inout')
    
    axs[1,1].plot(dist,freq_a,'-k',lw=0.5)
    axs[1,1].invert_xaxis()
    axs[1,1].set(xlabel='dF (aU)')
    axs[1,1].set_yticklabels([])
    axs[1,1].tick_params(axis='y',direction='inout')
    
    #waterfall plot
    axs[1,2].pcolormesh(t_a[l[0]:l[1]],freq_a,
                         deddata.T[:,l[0]:l[1]],shading='nearest')
    axs[1,2].set(xlabel='t (s)')
    axs[1,2].set_yticklabels([])
    axs[1,2].tick_params(axis='x',direction='in')
    axs[1,2].tick_params(axis='y',direction='inout')
    
    dwdm = 50
    dml = [peaks[n][0]-dwdm, peaks[n][0]+dwdm+1]
    if dml[0]<0: dml[0] = 0
    if dml[1]>len(dm_a): dml[1] = len(dm_a)
        
    #DM vs t (S) plot
    mesh = axs[3,2].pcolormesh(t_a[l[0]:l[1]],dm_a[dml[0]:dml[1]],
                np.max(cube[dml[0]:dml[1],:,l[0]:l[1]],axis=1),
                               shading='nearest')
    axs[3,2].get_shared_y_axes().join(axs[3,2],axs[3,1])
    axs[3,2].set_xticklabels([])
    axs[3,2].set_yticklabels([])
    axs[3,2].tick_params(axis='x',direction='inout')
    axs[3,2].tick_params(axis='y',direction='inout')
    axs[3,2].scatter(t_a[peaks[n][2]],dm_a[peaks[n][0]],
                     s=500,facecolors='none',edgecolors='r')
    
    #DM vs S plot
    axs[3,1].plot(cube[dml[0]:dml[1],peaks[n][1],peaks[n][2]],
                               dm_a[dml[0]:dml[1]],'-k',lw=0.5)
    axs[3,1].axvline(x=thres,color='red',lw=0.5)
    axs[3,1].tick_params(axis='x',direction='in')
    axs[3,1].tick_params(axis='y',direction='in')
    axs[3,1].set_xlabel('S')
    axs[3,1].set_ylabel('DM (cm$^{-3}$ pc)')
    
    #S vs t plot
    axs[4,2].plot(t_a[l[0]:l[1]],cube[peaks[n][0],peaks[n][1],
                                        l[0]:l[1]],'-k',lw=0.5)
    axs[4,2].axhline(y=thres,color='red',lw=0.5)
    axs[4,2].tick_params(axis='x',direction='in')
    axs[4,2].tick_params(axis='y',direction='in')
    axs[4,2].set_ylabel('S')
    axs[4,2].set(xlabel='t (s)')
    
    #W vs S plot
    axs[3,4].plot(cube[peaks[n][0],:,peaks[n][2]],
                            w_a*tsamp,'-k',lw=0.5)
    axs[3,4].axvline(x=thres,color='red',lw=0.5)
    axs[3,4].get_shared_y_axes().join(axs[3,4], axs[3,5])
    axs[3,4].set_yscale('log')
    axs[3,4].tick_params(axis='x',direction='in')
    axs[3,4].tick_params(axis='y',direction='in')
    axs[3,4].set_xlabel('S')
    axs[3,4].set_ylabel('W (s)')

    #W vs DM (S) plot
    mesh = axs[3,5].pcolormesh(dm_a[dml[0]:dml[1]],w_a*tsamp,
        np.max(cube[dml[0]:dml[1],:,l[0]:l[1]],axis=2).T,shading='nearest')
    axs[3,5].scatter(dm_a[peaks[n][0]],w_a[peaks[n][1]]*tsamp,
                     s=500,facecolors='none',edgecolors='r')
    axs[3,5].set_xticklabels([])
    axs[3,5].set_yticklabels([])
    axs[3,5].get_shared_x_axes().join(axs[3,5], axs[4,5])
    axs[3,5].tick_params(axis='x',direction='inout')
    axs[3,5].tick_params(axis='y',direction='inout')

    #S vs DM plot
    axs[4,5].plot(dm_a[dml[0]:dml[1]],cube[dml[0]:dml[1]
                        ,peaks[n][1],peaks[n][2]],'-k',lw=0.5)
    axs[4,5].axhline(y=thres,color='red',lw=0.5)
    axs[4,5].tick_params(axis='x',direction='in')
    axs[4,5].tick_params(axis='y',direction='in')
    axs[4,5].set_xlabel('DM (cm$^{-3}$ pc)')
    axs[4,5].set_ylabel('S')
    
    fig = plt.savefig(output_directory+output
                +'_C{:d}.jpeg'.format(n+1),dpi=dpi)
    #fig = plt.show()
    plt.close('all')
    return
    
tsamp = header['tsamp']
freq1 = header['fch1']
freqoff = header['foff']
nchans = header['nchans']

#downsampling
if progress: print('Downsampling,')
data = downsample(data,ds)
tsamp *= ds
tlen = np.shape(data)[0]

#frequency array
freq_a = freq1 + np.arange(nchans)*freqoff

#zapping frequency channels with RFI
if progress: print('Removing RFI and baseline,')
data = zap(data,zaplist,freq_a)

#zero dm matched filter
data = zdmf(data)

#remove baseline
data = rembaseline(data)

#dedispersion to timeseries
if progress: print('Performing dedispersion,')
dm_a = np.arange(dmstart,dmend+dminc,dminc)

timeseries = np.zeros((len(dm_a),tlen))
for dm in range(len(dm_a)):
    print('DM = {:.2f}'.format(dm_a[dm]),end='\r')
    timeseries[dm,:] = (dedisperse(data,dm_a[dm],
                                   freq_a,freq1,tsamp))
                                   
#initialize boxcar filter
if progress: print('Applying matched filter (boxcar),')
step = (snrloss+1)**4
w_a = np.arange(1,nbox+1)

if minw/tsamp > 1:
    w_a[0] = minw/tsamp
else:
    minw = tsamp
    
w_a[1:] = np.maximum(w_a[1:],
                        np.round(minw*step**(w_a[1:]-1)/tsamp))

w_a = w_a[w_a<=tlen]
nbox = len(w_a)

#apply boxcar filter as matched filter
cube = np.zeros((len(dm_a),nbox,tlen))
for dm in range(len(dm_a)):
    print('DM = {:.2f}'.format(dm_a[dm]),end='\r')
    for w in range(nbox):
        cube[dm,w,:] = lsr(timeseries[dm],int(w_a[w]))
        
#peak search
if progress: print('Searching for peaks,')
dof = np.sum(tlen*len(dm_a)/w_a)
thres = threshold(dof,prob)
peaks = peaksearch(cube, thres, maxpeak)

#main output plot
if progress: print('Plotting main analysis output,')
fig,axs = plt.subplots(7,3,figsize=(8.5,11),dpi=dpi, 
                       gridspec_kw={'height_ratios':[4,4,2,2,1,4,2],
                                    'width_ratios':[1,1,4],
                                    'hspace':0,'wspace':0})

fig.suptitle('Source: ' + str(header['source_name']) 
                + ', Time (MJD): {:.15f}'.format(header['tstart']))
fig.subplots_adjust(top=0.95)

for blank in [(1,0),(1,1),(2,0),(2,1),(3,0),(3,1),(4,0),(4,1),
                                        (4,2),(5,0),(6,0),(6,1)]:
    axs[blank].set_visible(False)

#arrays
t_a = np.arange(tlen)*tsamp

#frequency vs flux
axs[0,0].invert_xaxis()
axs[0,0].get_shared_y_axes().join(axs[0,0], axs[0,1], axs[0,2])
axs[0,0].tick_params(axis='y',direction='in')
axs[0,0].tick_params(axis='x',direction='in')
axs[0,0].set_xlabel('F (aU)')
axs[0,0].set_ylabel('f (MHz)')

#frequency vs dflux
axs[0,1].invert_xaxis()
axs[0,1].set_yticklabels([])
axs[0,1].tick_params(axis='x',direction='in')
axs[0,1].tick_params(axis='y',direction='inout')
axs[0,1].set_xlabel('dF (aU)')

#waterfall plot
axs[0,2].pcolormesh(t_a,freq_a,data.T,shading='nearest')
axs[0,2].set_xticklabels([])
axs[0,2].set_yticklabels([])
axs[0,2].get_shared_x_axes().join(axs[0,2], axs[1,2], axs[2,2], axs[3,2])
axs[0,2].tick_params(axis='x',direction='inout')
axs[0,2].tick_params(axis='y',direction='inout')
axs[0,2].set_xlim([0,t_a[-1]])

#DM vs t (S) plot
mesh = axs[1,2].pcolormesh(t_a,dm_a,np.max(cube,axis=1),
                           shading='nearest')
axs[1,2].set_xticklabels([])
axs[1,2].tick_params(axis='x',direction='inout')
axs[1,2].tick_params(axis='y',direction='in')
axs[1,2].set_ylabel('DM (cm$^{-3}$ pc)')

#S vs t plot
axs[2,2].plot(t_a,np.max(cube,axis=(0,1)),'-k',lw=0.5)
axs[2,2].axhline(y=thres,color='red',lw=0.5)
axs[2,2].set_xticklabels([])
axs[2,2].tick_params(axis='x',direction='inout')
axs[2,2].tick_params(axis='y',direction='in')
axs[2,2].set_ylabel('S')

#flux vs t plot
axs[3,2].tick_params(axis='x',direction='in')
axs[3,2].tick_params(axis='y',direction='in')
axs[3,2].set_xlabel('t (s)')
axs[3,2].set_ylabel('Flux (aU)')

#W vs S plot
axs[5,1].plot(np.max(cube,axis=(0,2)),
              w_a*tsamp,'-k',lw=0.5)
axs[5,1].axvline(x=thres,color='red',lw=0.5)
axs[5,1].get_shared_y_axes().join(axs[5,1], axs[5,2])
axs[5,1].set_yscale('log')
axs[5,1].tick_params(axis='x',direction='in')
axs[5,1].tick_params(axis='y',direction='in')
axs[5,1].set_xlabel('S')
axs[5,1].set_ylabel('W (s)')

#W vs DM (S) plot
mesh = axs[5,2].pcolormesh(dm_a,w_a*tsamp,
                    np.max(cube,axis=2).T,shading='nearest')
axs[5,2].set_xticklabels([])
axs[5,2].set_yticklabels([])
axs[5,2].get_shared_x_axes().join(axs[5,2], axs[6,2])
axs[5,2].tick_params(axis='x',direction='inout')
axs[5,2].tick_params(axis='y',direction='inout')

#S vs DM plot
axs[6,2].plot(dm_a,np.max(cube,axis=(1,2)),'-k',lw=0.5)
axs[6,2].axhline(y=thres,color='red',lw=0.5)
axs[6,2].tick_params(axis='x',direction='in')
axs[6,2].tick_params(axis='y',direction='in')
axs[6,2].set_xlabel('DM (cm$^{-3}$ pc)')
axs[6,2].set_ylabel('S')

if len(peaks)>0:
    #best timeseries
    axs[3,2].plot(t_a,timeseries[peaks[0][0]],'-k',lw=0.5)
    output = 'FRB_{:.15f}'.format(header['tstart'])
    #frequency profile
    fvsds = freqdist(data,dm_a[peaks[0][0]],w_a[peaks[0][1]]*tsamp,
                t_a[peaks[0][2]],freq_a,tsamp)
    axs[0,0].plot(np.cumsum(fvsds),freq_a,'-k',lw=0.5)
    axs[0,1].plot(fvsds,freq_a,'-k',lw=0.5)
    for peak in peaks:
        #candidate highlighting
        axs[1,2].scatter(peak[2]*tsamp,
               dm_a[peak[0]],s=500,facecolors='none',edgecolors='r')
        axs[5,2].scatter(dm_a[peak[0]],w_a[peak[1]]*tsamp,
                     s=500,facecolors='none',edgecolors='r')
        #FRB line
        frbdelay = dmdelay(dm_a[peak[0]],freq_a,
                           freq1)/tsamp
        line = (peak[2]-frbdelay)*tsamp
        axs[0,2].plot(line,freq_a,'k',lw=0.5,alpha=0.5)
else:
    output = 'non_{:.15f}'.format(header['tstart'])
    #frequency profile
    bandpass = np.sum(data,axis=0)
    axs[0,0].plot(np.cumsum(bandpass),freq_a,'-k',lw=0.5)
    axs[0,1].plot(bandpass,freq_a,'-k',lw=0.5)

fig.subplots_adjust(right=0.8)
mesh_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
fig.colorbar(mesh, cax=mesh_ax)
fig = plt.savefig(output_directory+output+'_Plot.jpeg',dpi=dpi)
#fig = plt.show()
plt.close('all')

#candidate plots
if len(peaks)>0:
    if progress: print('Plotting candidate(s) pulse profile,')
    fig, axs = plt.subplots(len(peaks),figsize=(8.5,11),dpi=dpi,
                            squeeze=False,sharex=True,
                            gridspec_kw={'hspace':0})
    fig.subplots_adjust(top=0.94,bottom=0.08)
    fig.suptitle('Candidates\nSource: ' + str(header['source_name']) 
                         + ', Time (MJD): {:.15f}'.format(header['tstart']))

    for i in range(len(peaks)):
        axs[i,0].plot(t_a,timeseries[peaks[i][0]],'-k',lw=0.5)
        axs[i,0].tick_params(axis='x',direction='in')
        axs[i,0].tick_params(axis='y',direction='in')
    fig.text(0.5, 0.04, 't (s)', ha='center')
    fig.text(0.04, 0.5, 'Flux (aU)', va='center', rotation='vertical')

    fig = plt.savefig(output_directory+output+'_Candidates.jpeg',dpi=dpi)
    #fig = plt.show()
    plt.close('all')
    
#dedispersed waterfall plotting
if len(peaks)>0:
    if progress: print('Plotting candidate(s) waterfall plot(s),')
    if waterfall < 0:
        for i in range(len(peaks)):
            waterfall_plot(i)
    else:
        waterfall_plot(waterfall)
        
#candidate properties list output
if progress: print('Writing analysis summary,')
cand = []
for i in range(len(peaks)):
    cand.append([i+1,dm_a[peaks[i][0]],w_a[peaks[i][1]]*tsamp,
                  t_a[peaks[i][2]],peaks[i][3],np.sqrt(peaks[i][3])])
c = open(output_directory+output+'_List.txt','w')
c.write('Source: '+str(header['source_name'])+'\nTime (MJD): ' 
         +'{:.15f}\n'.format(header['tstart'])
         +'\nDownsampling coefficient: '+'{:d}'.format(ds)
         +'\nStart: {:.8f}\nEnd: {:.8f}'.format(start,end)
         +'\nThreshold: {:.8f}\n\n'.format(thres))
c.write('Zapped frequencies:\n')
for i in zaplist:
    c.write('{:.2f} - {:.2f} MHz\n'.format(i[0],i[1]))
c.write('\n')
if len(peaks)>0:
    c.write(tabulate(cand,['Candidate','DM (cm-3 pc)','Width (s)',
            'Time (s)','S','SNR'],colalign=("center","center","center",
            "center","center","center")))
c.close()

data = None
cube = None
peaks = None
dm_a = None
w_a = None
t_a = None
          
if progress: print('Done!')