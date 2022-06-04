import numpy as np
from astropy.time import Time
import struct
from tabulate import tabulate

#fake parameter ranges
widthmin = 0.5 #minimum width of frb (ms)
widthmax = 5 #maximum width of frb (ms)
snrmin = 8 #minimum snr
snrmax = 50 #maximum snr
dmmin = 200 #minimum dispersion measure
dmmax = 1900 #maximum dispersion measure

#observation parameters
tobs = 0.5 #observation period (min)
tsamp = 64 #sampling time (us)
nifs = 1 #number of polarization channels
nchan = 512 #number of frequency channels
fstart = 1200.0 #frequency of first channel (MHz)
finc = 0.1 #bandwidth (MHz)

outnum = 20 #output file number

def send_string(string):
    length = len(string)
    fout.write(length.to_bytes(4,'little'))
    fout.write(bytearray(string,encoding="utf8"))

def send_int(name,integer):
    send_string(name)
    fout.write(integer.to_bytes(4,'little'))
    
def send_float(name,floating):
    send_string(name)
    fout.write(struct.pack('<d',floating))
    
def send_header(source):
    send_string('HEADER_START')
    send_string('source_name')
    send_string(source)
    send_int('machine_id',10)
    send_int('telescope_id',4)
    send_int('data_type',1)
    send_float('fch1',fstart)
    send_float('foff',-finc)
    send_int('nchans',nchan)
    send_int('nbits',32)
    send_float('tstart',Time.now().mjd)
    send_float('tsamp',tsamp*1e-6)
    send_int('nifs',nifs)
    send_string('HEADER_END')
    
def generate(width,dm,snr):
    data = np.zeros((len(t),nifs,nchan))
    delay = [4148.741601*((1.0/fstart/fstart)
                          -(1.0/n/n))*dm for n in f]
    maximum = tobs*60+delay[-1]
    rising = np.random.uniform(0,maximum)
    center = rising + width*0.5e-3
    trailing = rising + width*1e-3
    snr /= np.sqrt(width*nchan/(tsamp*1e-3)) #scaling over width and channels
    for i in range(len(t)): #looping over sampling time
        for j in range(nifs): #looping over desif
            for k in range(nchan): #looping over frequencies
                if t[i]+delay[k]>=rising and t[i]+delay[k]<=trailing:
                    data[i,j,k] = snr
    data += np.random.normal(0,1,(len(t),nifs,nchan))
    return data, center
    
def write(data):
    data = data.ravel()
    fout.write(struct.pack('<%sf' % len(data), *data))
    
#initialization for frequency, time
f = np.arange(fstart,fstart-nchan*finc,-finc)
t = np.arange(0,tobs*60,tsamp*1e-6)
prop = []

#looping over filterbanks
for i in range(outnum):
    print('FRB Fake: {:d}'.format(i),end='\r')
    #randomizer
    width = np.random.uniform(widthmin,widthmax)
    snr = np.random.uniform(snrmin,snrmax)
    dm = np.random.uniform(dmmin,dmmax)
    #initialize file
    output = 'fake_output' + str(i) + '.fil'
    source = 'PyFake DM = {:.5f}'.format(dm)
    fout = open(output,'wb')
    #writing header
    send_header(source)
    #writing fake pulsar
    data, center = generate(width,dm,snr)
    prop.append([i,center,width,snr,dm])
    write(data)
    fout.close()
    
#text file of fake properties
propout = open('fakeproperties.txt','w')
print(tabulate(prop,['Fake','Pulse location (s)','Width (ms)','SNR','DM (cm-3 pc)'],
         colalign=("center","center","center","center","center")))
propout.write(tabulate(prop,['Fake','Pulse location (s)','Width (ms)','SNR','DM (cm-3 pc)'],
         colalign=("center","center","center","center","center")))
propout.close()