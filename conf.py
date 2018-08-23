import os,sys
import time, numpy
import analysis.event
import analysis.stack
import analysis.pixel_detector
import analysis.beamline
import analysis.hitfinding
import analysis.sizing
import plotting.image
import plotting.line
import plotting.correlation
import analysis.patterson
import utils.reader
import utils.cxiwriter
import utils.array
import ipc.mpi
import numpy.random

this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(this_dir)

from backend import add_record

# Skip events with no pulse
do_skipdark = False
# Collect and write out stacks of frames
do_stacks = False
# Common mode correction along fastest changing dimension
do_cmc = True
# Send all events to the frontend
do_showall = True
# Show hybrid image
do_showhybrid = True
# Sizing of hits
do_sizing = True
# ToF analysis
do_tof = False
# MCP detector
do_MCP = False
# Save to file
do_save   = False
# Pattern analysis (for multiples)
do_patterson = False
# Validation (for amoc6914)
do_validation = False
# Use Gijs' analysis
i_like_flamingos = False
#Do hitfinding only on center speckle and perform GMD corrected hitfinding
mask_center_speckle = True
# do hitfinder with statistical hitfinder
do_stat_hitfinder = False
# GMD hitfinder
use_gmd_hitfinding = True

# Origin of data

# Switch between different datasets
data_mode = "amol3416"
#data_mode = "amol3116"
#data_mode = "amoc6914"
#data_mode = "simulation" # condor


## IMPORT FILES FOR STATISTICAL HITFINDER ##
#import h5py
#poiss_mask = h5py.File("/reg/data/ana14/amo/amol3416/scratch/hummingbird/stathitfinder/r0153/poisson_mask.h5",'r')['data'][:].astype(bool)
#sum_over_bkg_frames = h5py.File("/reg/data/ana14/amo/amol3416/scratch/hummingbird/photon_count_data/r0153_photon_data_done.h5",'r')['sum_over_frames'][-1]
#sample_params = h5py.File("/reg/data/ana14/amo/amol3416/scratch/hummingbird/stathitfinder/r0153/lambda_matrix.h5",'r')['params'][:] #fitting params of background
#thr_params = h5py.File("/reg/data/ana14/amo/amol3416/scratch/hummingbird/stathitfinder/r0150/thr_params.h5",'r')['data'][:]
#bag_bkg = h5py.File("/reg/data/ana14/amo/amol3416/scratch/hummingbird/stathitfinder/r0150/bag_bkg.h5",'r')['data'][:]
#fit_bkg = h5py.File("/reg/data/ana14/amo/amol3416/scratch/hummingbird/stathitfinder/r0153/lambda_matrix.h5",'r')['fit_bkg'][:]


# ------------------------------
# PSANA
# ------------------------------
if data_mode == "amol3416":
    import conf_amol3416 as conf
elif data_mode == "amol3116":
    import conf_amol3116 as conf
elif data_mode == "amoc6914":
    import conf_amoc6914 as conf
elif data_mode == "simulation":
    import conf_condor as conf
else:
    raise Exception("%s is not a valid data_mode." % data_mode)
state = conf.state

# Read events from idx file
idxfile = conf.idxfile
if idxfile is not None:
    import h5py
    idxdata = h5py.File(idxfile, "r")
    fiducials = idxdata['/fiducial'][()]
    times = idxdata['/timestamp2'][()]
    len = fiducials.shape[0]
    rankrange = numpy.arange(ipc.mpi.slave_rank(), len, ipc.mpi.nr_workers())
    state['fiducials'] = fiducials[rankrange]
    state['times'] = times[rankrange]

# Detector
front_type = conf.pnccd_type
front_key  = conf.pnccd_key

# MCP motors
back_type = conf.mcp_type
back_key  = conf.mcp_key

# Injector motors
injector_x_key = conf.injector_x_key
injector_z_key = conf.injector_z_key

# Nozle pressures
nozzle_pressure_1_key = conf.nozzle_pressure_1_key
nozzle_pressure_2_key = conf.nozzle_pressure_2_key

# Injector pressures
injector_pressure_1_key = conf.injector_pressure_1_key
injector_pressure_2_key = conf.injector_pressure_2_key
injector_pressure_3_key = conf.injector_pressure_3_key

# Injector voltage
injector_voltage_key = conf.injector_voltage_key

# Electrospray
electrospray_key = conf.electrospray_key

# Tof detector
acqiris_type = conf.acqiris_type
acqiris_key  = conf.acqiris_key

# KB settings (for amol3416)
if data_mode == "amol3416":
    kb_04_key = conf.kb_04_key
    kb_05_key = conf.kb_05_key
    kb_06_key = conf.kb_06_key
    kb_07_key = conf.kb_07_key
    kb_08_key = conf.kb_08_key
    kb_09_key = conf.kb_09_key
    kb_10_key = conf.kb_10_key
    # Horizontal focusing upstream bend
    kb_11_key = conf.kb_11_key
    # Horizontal focusing downstream bend
    kb_12_key = conf.kb_12_key
    # Vertical focusing upstream bend
    kb_13_key = conf.kb_13_key
    # Vertical focusing downstream bend
    kb_14_key = conf.kb_14_key

# Gain modes (gain mode gets overwritten if available in datastream)
#gain_mode = 6 # 1/1 (highest gain)
#gain_mode = 5 # 1/4
#gain_mode = 4 # 1/16
gain_mode = 3 # 1/64
gain_mode_key = conf.gain_mode_key

# Make a list of important parameters
param_dict = {'InjPosx': injector_x_key, 
              'InjPosz': injector_z_key,
              'NozzlePress1': nozzle_pressure_1_key,
              'NozzlePress2': nozzle_pressure_2_key,
              'InjPress1': injector_pressure_1_key,
              'InjPress2': injector_pressure_2_key,
              'InjVolt':injector_voltage_key,
              'Electrospray':electrospray_key,
              'GainMode': gain_mode_key}

# Backgrounds
# -----------
Nbg   = 20
rbg   = 20
obg   = 20
#Nbg   = 10
#rbg   = 2000
#obg   = 2000
bg_front = analysis.stack.Stack(name="bg_front",maxLen=Nbg,outPeriod=obg,reducePeriod=rbg)
bg_dir = "/reg/neh/home/hantke/amol3116/stacks/"

# Masks
# -----
M_front    = utils.reader.MaskReader(this_dir + "/mask/mask_front.h5","/data/data")
mask_front = M_front.boolean_mask
(ny_front,nx_front) = mask_front.shape
M_front_ext    = utils.reader.MaskReader(this_dir + "/mask/mask_front_ext4.h5","/data/data")
mask_front_ext = M_front_ext.boolean_mask

dx=conf.dx_front
dy=conf.dy_front

# Mask out a circle in the center
mask_front_circle = numpy.ones((ny_front+abs(dy), nx_front+abs(dx)), dtype=numpy.bool)
radius = 100#50#80#130
cx = 15
cy = 10
xx,yy = numpy.meshgrid(numpy.arange(nx_front+abs(dx)), numpy.arange(ny_front+abs(dy)))
rr = (xx-nx_front/2-cx)**2 + (yy-ny_front/2-cy)**2 >= (radius**2)
mask_front_circle &= rr

mask_front_edge = numpy.ones((ny_front+abs(dy), nx_front+abs(dx)), dtype=numpy.bool)
radius = 400
xx2,yy2 = numpy.meshgrid(numpy.arange(nx_front+abs(dx)), numpy.arange(ny_front+abs(dy)))
rr2 = (xx2-nx_front/2-cx)**2 + (yy2-ny_front/2-cy)**2 >= (radius**2)
mask_front_edge &= rr2
#mask_front_edge &= mask_front

# Mask sizing

# For TBSV (masks out low res)
radius_fit = 80
rr3 = (xx-nx_front/2-cx)**2 + (yy-ny_front/2-cy)**2 >= (radius_fit**2)
mask_front_fit = numpy.ones((ny_front+abs(dy), nx_front+abs(dx)), dtype=numpy.bool)
mask_front_fit &= rr3

# For carboxysomes (masks out high res)
#mask_front_fit = numpy.zeros((ny_front+abs(dy), nx_front+abs(dx)), dtype=numpy.bool)
#radius = 130
#cx = 40
#cy = 30
#xx,yy = numpy.meshgrid(numpy.arange(nx_front+abs(dx)), numpy.arange(ny_front+abs(dy)))
#rr = (xx-nx_front/2-cx)**2 + (yy-ny_front/2-cy)**2 <= (radius**2)
#mask_front_fit |= rr

# Mask out some pixel in the top
top_mask = numpy.ones_like(mask_front)
top_mask[:30] = False
#mask_front &= top_mask

# Mask a bad panel
bad_panel_mask = numpy.ones_like(mask_front)
bad_panel_mask[420:512,512+70:512+70+60] = False
mask_front &= bad_panel_mask

# Mask a bad panel2
bad_mask = numpy.ones_like(mask_front)
bad_mask[:512,512:512+256] = False
#mask_front &= bad_mask

# Gain calibration
# ----------------
tmp = numpy.loadtxt("/reg/d/psdm/amo/amol3416/calib/PNCCD::CalibV1/Camp.0:pnCCD.0/pixel_gain/0-end.data")
#tmp = numpy.loadtxt("/reg/d/psdm/AMO/amo06516/scratch/yoon82/run15_pnccdFront_flatIso_clean.data")
gain_cal = numpy.zeros(shape=(1024, 1024))
gain_cal[:,:512] = tmp[:1024,:]
gain_cal[:,512:] = tmp[1024:,:]
gain_cal = gain_cal[::-1,:]

if do_save and ipc.mpi.is_worker():
    w_dir = '/reg/d/psdm/amo/amol3416/scratch/hummingbird/runs'
    W = utils.cxiwriter.CXIWriter(w_dir + "/r%04d.h5" %conf.run_nr, chunksize=1)

# ---------------------------------------------------------
# PLOTTING PARAMETER
# ---------------------------------------------------------

# Injector position limits
x_min = 1.0025
x_max = 2.0025
x_bins = 200
z_min = -8.
z_max = -6.
z_bins = 100

# Hitrate mean map 
hitrateMeanMapParams = {
    'xmin': x_min,
    'xmax': x_max,
    'ymin': z_min,
    'ymax': z_max,
    'xbins': x_bins,
    'ybins': z_bins,
    'xlabel': 'Injector Position in x (%s)' % injector_x_key,
    'ylabel': 'Injector Position in z (%s)' % injector_z_key
}

event_number = 0
hitcounter = 0

# ---------------------------------------------------------
# ANALYSIS
# ---------------------------------------------------------

# Hit finding
# -----------
# Simple hitfinding (Count Nr. of lit pixels)
GMD_corrected_threshold = 375#800#400#220 #250
GMD_corrected_threshold_strong = 600
hitscoreThreshold = 3400
hitscoreThresholdStrong = 3000
if data_mode == "simulation":
    gain_mode = 0
    hitscoreThreshold = 1000

# Photon energy [eV]
photonEnergy = 800.

# More advanced hitfinding
polynomial = [200]
    
# Sizing
# ------
binning = 4

centerParams = {
    'x0'       : (512 + 10 - (nx_front-1)/2.)/binning,
    'y0'       : (512 + 20 - (ny_front-1)/2.)/binning,
    'maxshift' : int(numpy.ceil(10./binning)),
    'threshold': 1,
    #'threshold': 20*binning**2,
    'blur'     : 4,
}
pixelsize_native = 75E-6 
modelParams = {
    #'wavelength': 2.480, # in nm, 500 eV
    'wavelength': 1.550, # in nm, 800 eV
    'pixelsize': pixelsize_native/1E-6*binning,
    #'distance': 380., # in mm
    'distance': 130., #230., #121., # in mm
    'material': 'sucrose',
}
sizingParams = {
    'd0':200., # in nm
    'i0':1., # in mJ/um2
    'brute_evals':10,
}

# Physical constants
h = 6.62606957e-34 #[Js]
c = 299792458 #[m/s]
hc = h*c  #[Jm]
eV_to_J = 1.602e-19 #[J/eV]

#res = modelParams["distance"] * 1E-3* modelParams["wavelength"] * 1E-9 / ( pixelsize_native * nx_front )
#expected_diameter = 150

# Thresholds for good sizing fits
fit_error_threshold = 2.6E-3#4.0e-3
photon_error_threshold = 3000
diameter_min = 50  #[nm]
diameter_max = 150 #[nm]

# Patterson
# ---------
patterson_threshold = 5.
patterson_floor_cut = 50.
patterson_mask_smooth = 5.
patterson_diameter = 60. #4 * expected_diameter * 1E-9 / res
multiple_score_threshold = 200.

# ---------------------------------------------------------
# INTERFACE
# ---------------------------------------------------------

# Image
vmin_front = 0
vmax_front = 1000

# Radial averages
radial_tracelen = 300

if i_like_flamingos:
    # Only import flamingo if selected
    import redflamingo
    redflamingo.read_mask()
    redflamingo.read_coordinates()

#GMD corrected hitscore parameters applied to center speckle:
center_x = 553
center_y = 545
center_speckle_radius = 200 #160 #230

pfit = None
buffer_size = 20
gmd_buffer = numpy.zeros(buffer_size)        
hitscore_buffer = numpy.ones(buffer_size)


# ---------------------------------------------------------
# E V E N T   C A L L
# ---------------------------------------------------------
lastsent = numpy.random.uniform(-100,0)
lastsentnc = numpy.random.uniform(-100,0)
lastsentmultiple = numpy.random.uniform(-100,0)

def onEvent(evt):
    #os.system("echo $HOSTNAME")
    global event_number
    global hitcounter
    global lastsent
    global lastsentnc
    global lastsentmultiple
    event_number += 1

    global pfit
    global gmd_buffer
    global hitscore_buffer
    global buffer_size
    
    # MPI
    main_worker = ipc.mpi.is_main_worker()
    rank = ipc.mpi.rank
    
    # ------------------- #
    # INITIAL DIAGNOSTICS #
    # ------------------- #

    # Time measurement
    analysis.event.printProcessingRate()

    # Plot Fiducials
    #plotting.line.plotTimestamp(evt["eventID"]["Timestamp"])
    #print evt['camera'].keys()
    # --------------------------------------------------------
    
    # Check if all paramters are in the data stream
    param_in_stream = [(p in evt["parameters"].keys()) for p in param_dict.values()]
    if all(param_in_stream):
        msg = "All parmeters are stored. \nKeeping fingers crossed!"
    else:
        msg = "Hmm, something is missing. \nWe see %d/%d parameters!" %(sum(param_in_stream), len(param_in_stream))
    for k,p in param_dict.iteritems():
        plotting.line.plotHistory(add_record(evt["stream"], "stream", p, (p in evt["parameters"].keys())),
                                  msg=msg, alert=not (p in evt["parameters"].keys()), name="%s/%s" %(k,p), group="Data stream")
    
    # --------------------------------------------------------
    
    # Average pulse energies, in case there is no pulse energy it returns None
    analysis.beamline.averagePulseEnergy(evt, evt["pulseEnergies"]) 

    # Skip frames that do not have a pulse energy
    try:
        evt['analysis']['averagePulseEnergy']
        pulse = True
    except (TypeError,KeyError):
        pulse = False
    if do_skipdark and not pulse:
        print "No pulse energy. Skipping event."
        return

    # Average pulse energies, in case there is no pulse energy it returns None
    analysis.beamline.averagePhotonEnergy(evt, evt["photonEnergies"]) 
    add_record (evt['analysis'], 'analysis', 'PhotonEnergy', photonEnergy)

    # Skip frames that do not have a photon energy
    try:
        evt['analysis']['averagePhotonEnergy']
    except (TypeError,KeyError):
        print "No photon energy. Skipping event."
        return
    
    # Skip frames that do not have the pnCCDs
    try:
        evt[front_type][front_key]
    except (TypeError,KeyError):
        print "No front pnCCD. Skipping event."
        return

    # Plotting of average pulse energy
    if pulse:
        plotting.line.plotHistory(evt["analysis"]["averagePulseEnergy"], group='Diagnostics')

    # ------------------------------------------------------------------

    # Gain of PNCCD
    try:
        gain_mode = evt['parameters'][gain_mode_key].data
    except:
        print "Could not read gain mode from data stream, using default: %d" %gain_mode
    #analysis.pixel_detector.pnccdGain(evt, evt['analysis']['averagePhotonEnergy'],gain_mode)
    analysis.pixel_detector.pnccdGain(evt, evt['analysis']['PhotonEnergy'],gain_mode)
    aduThreshold = evt['analysis']['gain'].data/2.

    if (event_number % 100 == 0):
        print "In current gain mode (%d), the gain is %.2f"  %(gain_mode, evt['analysis']['gain'].data)
                                      
    # Common mode correction for PNCCD
    if data_mode != "simulation" and do_cmc:
        analysis.pixel_detector.commonModePNCCD(evt, front_type, front_key, transpose=False, signal_threshold=aduThreshold)
        front_type_s = "analysis"
        front_key_s  = "corrected - " + front_key
    else:
        front_type_s = front_type
        front_key_s  = front_key
    front_no_geometry = evt[front_type_s][front_key_s]
    front = evt[front_type_s][front_key_s]
    
    integrated_front=numpy.sum(front.data)
    integrated_front_rec=add_record(evt['analysis'], 'analysis', 'integrated_front', integrated_front)
    plotting.line.plotHistory(evt['analysis']['integrated_front'], name= 'integrated front' ,group='Diagnostics')

    # Fix gain difference for top right quad
    #front.data[:512,512:] *= 2
    #front.data /= gain_cal
    #front.data = gain_cal

    # Move right detector half (assemble)
    if data_mode in ["amol3416","amol3116","amoc6914"]:
        front = analysis.pixel_detector.moveHalf(evt, front, horizontal=conf.dx_front, vertical=conf.dy_front, outkey='data_half-moved')
        front_type_s = "analysis"
        front_key_s  = "data_half-moved"
        mask_front_s = analysis.pixel_detector.moveHalf(evt, add_record(evt["analysis"], "analysis", "mask", mask_front), horizontal=conf.dx_front, vertical=conf.dy_front, outkey='mask_half-moved' ).data
        mask_front_ext_s = analysis.pixel_detector.moveHalf(evt, add_record(evt["analysis"], "analysis", "mask_ext", mask_front_ext), horizontal=conf.dx_front, vertical=conf.dy_front, outkey='mask_half-moved_ext' ).data
        d = 35
        mask_front_ext_s[512-d:512+d,:] = 0
    else:
        mask_front_s = evt[conf.mask_type][conf.mask_key].data == 0
        mask_front_ext_s = evt[conf.mask_type][conf.mask_key].data == 0

    # Construct assembled mask for fitting
    mask_front_fit_s = mask_front_s & mask_front_fit

    #Multiply masks with circle masks
    mask_front_s &= mask_front_circle
    mask_front_ext_s &= mask_front_circle
    
    mask_front_edge_s=mask_front_s.copy()
    mask_front_edge_s &= mask_front_edge
    # Plot PNCCD image
    if do_showall and (event_number % 100 == 0):# and rank > 1:
        plotting.image.plotImage(evt[front_type_s][front_key_s], group='Diagnostics', 
                                 msg='', name="pnCCD front", mask=mask_front_s)
        plotting.line.plotHistogram(evt[front_type_s][front_key_s], group='Diagnostics', 
                                    mask=mask_front_s, hmin=0, hmax=200, bins=200,
                                    vline=aduThreshold, 
                                    label='PNCCD [ADU]', name='PNCCD: pixel values [ADU]')
    # Plot MCP image
    
    if do_MCP and do_showall and (event_number % 10 == 0):# and rank > 1:
        plotting.image.plotImage(evt[back_type][back_key], group='Diagnostics', 
                                 msg='', name="MCP Back") #, mask=mask_front_s)
        plotting.line.plotHistogram(evt[front_type_s][front_key_s], group='Diagnostics', 
                                    mask=mask_front_s, hmin=0, hmax=200, bins=200,
                                    vline=aduThreshold, 
                                    label='PNCCD [ADU]', name='PNCCD: pixel values [ADU]')
   
    # ------------------------------------------------------------------

    # Initialize hit flags
    hit = False
    multiple_hit = False
    correctsized_hit = False

    # ------------------------------------------------------------------

    # Finding hits
    
    mask_front_2 = analysis.hitfinding.generate_radial_mask(mask_front_s,center_x,center_y,center_speckle_radius)
    mask_front_2[:,512:] = 0

    if mask_center_speckle:
        analysis.hitfinding.countLitPixels(evt, front, aduThreshold=aduThreshold, 
                                           hitscoreThreshold=hitscoreThreshold, mask=mask_front_2)
        analysis.hitfinding.countLitPixels(evt, front, aduThreshold=aduThreshold,
                                           hitscoreThreshold=hitscoreThreshold, mask=mask_front_edge_s, outkey= "litpixel_edge: ")
    else:
        analysis.hitfinding.countLitPixels(evt, front, aduThreshold=aduThreshold, 
                                       hitscoreThreshold=hitscoreThreshold, mask=mask_front_s)

    hitscore = evt["analysis"]["litpixel: hitscore"].data
    hit = evt["analysis"]["litpixel: isHit"].data
    strong_hit = hitscore > hitscoreThresholdStrong
    hitscore_edge= evt["analysis"]["litpixel_edge: hitscore"].data
    if hitscore_edge>3:
        hitscore_edge_corrected=hitscore*1./(hitscore_edge*1.)
        add_record(evt["analysis"], "analysis", "hitscore_edge_corrected", hitscore_edge_corrected) 
        plotting.line.plotHistory(evt["analysis"]["hitscore_edge_corrected"], runningHistogram=True,
                                  hmin=0, hmax=100000, bins=100, window=100, history=1000, group='Alignment')
    if hitscore_edge>5:
        plotting.line.plotHistory(evt["analysis"]["litpixel_edge: hitscore"], runningHistogram=True,
                              hmin=0, hmax=100000, bins=100, window=100, history=1000, group='Alignment')

    
    if mask_center_speckle and len(evt['pulseEnergies'].values())>0:
        gmd_evt = evt['pulseEnergies'].values()[-1].data
        gmd_buffer[event_number % buffer_size] = gmd_evt
        hitscore_buffer[event_number % buffer_size] = hitscore
        if event_number > buffer_size-1:
            if (event_number % buffer_size) == 0:
                pfit = numpy.poly1d(numpy.polyfit(gmd_buffer[gmd_buffer > 0.5],hitscore_buffer[gmd_buffer > 0.5],1))
            if pfit:
                if gmd_evt<0.5:
                    normalised_hitscore=0
                else:
                    normalised_hitscore = hitscore - pfit(gmd_evt)
                add_record(evt["analysis"], "analysis", "norm_hitscore", normalised_hitscore)
                plotting.line.plotHistory(evt["analysis"]["norm_hitscore"], runningHistogram=True,
                                          hmin=0, hmax=100000, bins=100, window=100, history=1000,

                                          hline=GMD_corrected_threshold, group='Finding hits')
                
                # Use GMD corrected hitscore to calculate if there is a hit
                hit = normalised_hitscore > GMD_corrected_threshold
                strong_hit= normalised_hitscore>GMD_corrected_threshold_strong
    evt["analysis"]["litpixel: isHit"].data = int(hit)

    #if use_gmd_hitfinding:
    #    add_record(evt["analysis"], "analysis", "norm hitscore: isHit", hit)
    #    plotting.line.plotHistory(evt["analysis"]["norm hitscore: isHit"], name="isHit", group='Hit rates')  
    #else:
    plotting.line.plotHistory(evt["analysis"]["litpixel: isHit"], name="isHit", group='Hit rates')  

    if do_stat_hitfinder:
        pulse_energy = evt["analysis"]["averagePulseEnergy"].data
        analysis.hitfinding.photon_count_frame(evt,front_type_s,front_key_s,2.*aduThreshold)
        analysis.hitfinding.lambda_values(evt,pulse_energy,sum_over_bkg_frames,fit_bkg,sample_params)
        analysis.hitfinding.baglivo_score(evt,poiss_mask)
        analysis.hitfinding.stat_hitfinder(evt,pulse_energy,thr_params,bag_bkg)
        hitscore = evt["analysis"]["baglivo_score"].data
        hit = evt["analysis"]["bagscore: isHit"].data

        plotting.correlation.plotScatter(evt["analysis"]["averagePulseEnergy"], evt["analysis"]["baglivo_score"], 
                                         name='Baglivo score vs GMD', group='Finding hits', history=10000)
        plotting.correlation.plotScatter(evt["analysis"]["averagePulseEnergy"], evt["analysis"]["bagscore: threshold"], 
                                         name='Threshold vs GMD', group='Finding hits', history=10000)
     
        if hit:
            plotting.image.plotImage(evt[front_type_s][front_key_s], group='Finding hits', send_rate=1,
                                     msg='', name="pnCCD front (Baglivo hit)")
       


    if hit: hitcounter += 1

    if hit:
        plotting.image.plotImage(evt[front_type_s][front_key_s], group='Finding hits', send_rate=1,
                                 msg='', name="pnCCD front (hit)", mask=mask_front_s)
        #plotting.image.plotImage(evt[back_type][back_key], group='Finding hits', 
        #                         msg='', name="MCP Back") #, mask=mask_front_s)
        

    if strong_hit:
        plotting.image.plotImage(evt[front_type_s][front_key_s], group='Finding hits', send_rate=1, 
                                 msg='', name="pnCCD front strong (hit)")#, mask=mask_front_s)

    if hit:#(event_number % 100 == 0):#hit:
        photons = add_record(evt["analysis"], "analysis", "photons", numpy.round(evt[front_type_s][front_key_s].data/evt['analysis']['gain'].data))
        plotting.image.plotImage(photons, group='Finding hits', send_rate=1,
                                 msg='', name="pnCCD front photons (hit)", mask=mask_front_s)
        import spimage
        photons_ds = add_record(evt["analysis"], "analysis", "photons downsampled", spimage.downsample(photons.data, 8, mode="integrate", mask=mask_front_s)[0])
        plotting.image.plotImage(photons_ds, group='Finding hits', send_rate=1,
                                 msg='', name="pnCCD front downsampled photons (hit)")

        
    # Plot the hit events
    #plotting.line.plotHistory(evt["analysis"]["litpixel: isHit"], runningHistogram=True, 
    #                          hmin=0, hmax=1, bins=100, window=100, history=1000, group='Finding hits')

    # Plot the hitscore
    plotting.line.plotHistory(evt["analysis"]["litpixel: hitscore"], runningHistogram=True, 
                              hmin=0, hmax=100000, bins=100, window=100, history=1000, 
                              hline=hitscoreThreshold, group='Finding hits')

    # ------------------------------------------------------------------

    # Finding hits (more advanced)

    # Define a photon score
    #analysis.pixel_detector.totalNrPhotons(evt, front, aduPhoton=gain, aduThreshold=0.5, outkey="nrPhotons")
    #photon_score = evt['analysis']['nrPhotons']

    # ...here one can implement a more fancy photon score

    # Finding the hits based on score and gmd
    #analysis.hitfinding.countPhotonsAgainstEnergyPolynomial(evt, photon_score, evt['analysis']['averagePulseEnergy'],
    #                                                        energyPolynomial = polynomial, outkey='photons: '):
    
    #hitscore = evt["analysis"]["photons: hitscore"].data
    #hit = evt["analysis"]["photons: isHit"].data
    #if hit: hitcounter += 1


    # Show hits
#    if hit:
#        clustermsg = ""
#        if cluster_source_on:
#            clustermsg = "WITH CLUSTER"
#        if lastsent - event_number < numpy.random.uniform(-250,0) and cluster_source_on:
#            lastsent = event_number
#            plotting.image.plotImage(evt[front_type_s][front_key_s], group='Finding hits' ,
#                                 msg=clustermsg, name="pnCCD front (hit)") #, mask=mask_front_s) 
#        if lastsentnc - event_number < numpy.random.uniform(-250,0) and not cluster_source_on:
#            lastsentnc = event_number
#            plotting.image.plotImage(evt[front_type_s][front_key_s], group='Finding hits' ,
#                                     msg=clustermsg, name="pnCCD front (hit, no cluster)") #, mask=mask_front_s) 
#        if hitscore > 550000 and not True:
#            plotting.image.plotImage(evt[back_type][back_key], group='Finding hits',
#                                 msg=clustermsg, name="MCP Back") #, mask=mask_front_s)
#            plotting.image.plotImage(evt[front_type_s][front_key_s], group='Finding hits' ,
#                                     msg=clustermsg, name="pnCCD front (super hits)") #, mask=mask_front_s) 


    # -----------------------------------------------------------------------

    # Find multiple hits
    if do_patterson:
        analysis.patterson.patterson(evt, front_type_s, front_key_s, mask_front_ext_s, 
                                     floor_cut=patterson_floor_cut, mask_smooth=patterson_mask_smooth,
                                     threshold=patterson_threshold ,diameter_pix=patterson_diameter,
                                     crop=512, full_output=True) 
        multiple_hit = evt["analysis"]["multiple score"].data > multiple_score_threshold
        
        # Plot the multiple hitscore
        if do_showall:
            plotting.line.plotHistory(evt["analysis"]["multiple score"], history=1000, name='Multiple score (patterson)',
                                      hline=multiple_score_threshold, group="Finding multiples")
            if multiple_hit:
                plotting.image.plotImage(evt['analysis']['patterson'], group='Finding multiples',
                                     msg='', name='Patterson')
                plotting.image.plotImage(front, group='Finding multiples',
                                     msg='', name='Front PNCCD (hits)')

    # -----------------------------------------------------------------------

    # Find right-sized singles
    if strong_hit and (not multiple_hit) and do_sizing:
    #if hit and (not multiple_hit) and do_sizing:

        # Comment out for fitting mask (not for sucrose or rubisco)
        #mask_front_fit_s = mask_front_s

        #print mask_front_fit_s.shape, mask_front_s.shape
        #print mask_front_fit_s.dtype, mask_front_s.dtype
        #print mask_front_fit_s.sum()
        #mask_front_s = mask_front_fit_s

        # Feed in wavelength (based on photon energy from datastream)
        photon_energy = evt['analysis']['averagePhotonEnergy'].data
        #modelParams['wavelength'] = hc / photon_energy / eV_to_J * 1e9 # in [nm]

        # Feed in gain
        modelParams['adu_per_photon'] = evt['analysis']['gain'].data
        
        # Crop to 1024 x 1024
        analysis.pixel_detector.cropAndCenter(evt, evt[front_type_s][front_key_s], w=1024, h=1024, outkey='data-cropped')
        front_key_s = "data-cropped"
        analysis.pixel_detector.cropAndCenter(evt, add_record(evt["analysis"], "analysis", "mask", mask_front_fit_s), w=1024, h=1024, outkey='mask-cropped')
        mask_front_fit_s = evt['analysis']['mask-cropped'].data

        # Binning
        analysis.pixel_detector.bin(evt, front_type_s, front_key_s, binning, mask_front_fit_s)
        mask_front_b = evt["analysis"]["binned mask - " + front_key_s].data
        front_type_b = "analysis"
        front_key_b = "binned image - " + front_key_s
        
        # CENTER DETERMINATION
        analysis.sizing.findCenter(evt, front_type_b, front_key_b, mask=mask_front_b, **centerParams)
        # RADIAL AVERAGE
        analysis.pixel_detector.radial(evt, front_type_b, front_key_b, mask=mask_front_b, cx=evt["analysis"]["cx"].data, cy=evt["analysis"]["cy"].data)          
        # FIT SPHERE MODEL
        analysis.sizing.fitSphereRadial(evt, "analysis", "radial distance - " + front_key_b, "radial average - " + front_key_b, **dict(modelParams, **sizingParams))
        # DIFFRACTION PATTERN FROM FIT
        analysis.sizing.sphereModel(evt, "analysis", "offCenterX", "offCenterY", "diameter", "intensity", (ny_front/binning,nx_front/binning), poisson=False, **modelParams)
        # RADIAL AVERAGE FIT
        analysis.pixel_detector.radial(evt, "analysis", "fit", mask=mask_front_b, cx=evt["analysis"]["cx"].data, cy=evt["analysis"]["cy"].data)
        # ERRORS
        analysis.sizing.photon_error(evt, front_type_b, front_key_b, "analysis", "fit", 144.)
        analysis.sizing.absolute_error(evt, front_type_b, front_key_b, "analysis", "fit", "absolute error")

        if data_mode == "simulation":
            # Compare with simulated values
            cx_sim = evt["particleParameters"]["cx"].data
            cx_rec = evt["analysis"]["cx"].data * binning
            cx_err = cx_rec-cx_sim
            if abs(cx_err) > binning:
                print "CENTER X RECOVERY FAILED: d_cx = %.1f (binning = %i)" % (round(cx_err,1), binning)
            cy_sim = evt["particleParameters"]["cy"].data
            cy_rec = evt["analysis"]["cy"].data * binning
            cy_err = cy_rec-cy_sim
            if abs(cy_err) > binning:
                print "CENTER Y RECOVERY FAILED: d_cy = %.1f (binning = %i)" % (round(cy_err,1), binning)
            d_sim = evt["particleParameters"]["particleDiameter"].data
            d_rec = evt["analysis"]["diameter"].data * 1E-9
            d_err = d_rec-d_sim
            if abs(d_err/d_sim) > 0.01:
                print "SIZE RECONVERY FAILED: d_d = %.1f nm (d = %.1f nm)" % (round(d_err/1E-9,1),
                                                                              round(d_sim/1E-9,1))
            i_sim = evt["particleParameters"]["particleIntensity"].data 
            i_rec = evt["analysis"]["intensity"].data * (1E-3/1E-12)
            i_err = i_rec-i_sim
            if abs(i_err/i_sim) > 0.1: # Center errors seem to propagate badly into the intensity determination,
                                       # this is why I allow myself to use here such a generous error threshold
                print "INTENSITY RECONVERY FAILED: d_i = %.3f mJ/um2 (i = %.3f mJ/um2)" % (round(i_err/(1E-3/1E-12),3),
                                                                                           round(i_sim/(1E-3/1E-12),3))           # Image of fit
        msg = "diameter: %.2f nm \nIntensity: %.4f mJ/um2\nError: %.2e; %.2e; %.2e" %(evt["analysis"]["diameter"].data, evt["analysis"]["intensity"].data, evt["analysis"]["photon error"].data, evt["analysis"]["absolute error"].data.sum(), evt["analysis"]["fit error"].data)
             
        # Selection of good fits
        small_fit_error    = evt['analysis']['fit error'].data    < fit_error_threshold
        small_photon_error = evt['analysis']['photon error'].data < photon_error_threshold
        correctsized_hit = small_fit_error and small_photon_error

        # Select only events in a certain diameter window
        diameter = evt['analysis']['diameter'].data
        #correctsized_hit &= ((diameter > diameter_min) & (diameter < diameter_max))

        # Plot errors
        plotting.line.plotHistory(evt["analysis"]["fit error"], history=1000, hline=fit_error_threshold, group='Sizing')
        plotting.line.plotHistory(evt["analysis"]["photon error"], history=1000, hline=photon_error_threshold, group='Sizing')
        #plotting.line.plotHistory(evt["analysis"]["absolute error"], history=1000, group='Sizing')

        if do_showhybrid:
            # HYBRID PATTERN
            hybrid = evt["analysis"]["fit"].data.copy()
            hybrid[:,512/binning:] = evt[front_type_b][front_key_b].data[:,512/binning:]
            add_record(evt["analysis"], "analysis", "Hybrid pattern", hybrid)

            if correctsized_hit:
                plotting.image.plotImage(evt["analysis"]["Hybrid pattern"], 
                                         name="Hybrid pattern (good)", vmin=evt['analysis']['gain'].data, 
                                         msg=msg, group='Sizing')

                plotting.image.plotImage(evt["analysis"]["Hybrid pattern"], 
                                         name="Hybrid pattern (good) / log", vmin=evt['analysis']['gain'].data, vmax=1e4, log=True, 
                                         msg=msg, group='Sizing')


            else:
                plotting.image.plotImage(evt["analysis"]["Hybrid pattern"], mask=mask_front_b, 
                                         name="Hybrid pattern (bad)", vmin=evt['analysis']['gain'].data, vmax=1e4, log=True, 
                                         msg=msg, group='Sizing')
            
        if do_showall:

            error = evt["analysis"]["photon error"].data
            #plotting.image.plotImage(evt["analysis"]["Hybrid pattern"], mask=mask_front_b, name="Hybrid pattern", 
            #                         vmin=20, vmax=1e4, log=True, msg=msg, group='Sizing')
        
            #plotting.image.plotImage(evt["analysis"]["fit"], log=True, mask=mask_front_b, name="pnCCD: Fit result (radial sphere fit)", msg=msg)
            # Plot measurement radial average
            plotting.line.plotTrace(evt["analysis"]["radial average - " + front_key_b],
                                    evt["analysis"]["radial distance - " + front_key_b],
                                    tracelen=radial_tracelen, group='Sizing')
            # Plot fit radial average
            plotting.line.plotTrace(evt["analysis"]["radial average - fit"],
                                    evt["analysis"]["radial distance - fit"],
                                    tracelen=radial_tracelen, group='Sizing')
            #plotting.line.plotHistory(evt["analysis"]["diameter"], history=1000, name_extension="single hits")

        if correctsized_hit:

            # Plot Correct sized hits
            plotting.image.plotImage(evt[front_type_s][front_key_s], group='Sizing', 
                                 msg=msg, name="pnCCD front (correct hit)")#, mask=mask_front_s)

            # Plot Intensity
            plotting.line.plotHistory(evt["analysis"]["intensity"], history=10000, 
                                      name ='Intensity (from sizing)', group='Results')
            # Plot size (in nm)
            plotting.line.plotHistory(evt["analysis"]["diameter"], history=10000, 
                                      name = 'Size in nm (from sizing)', group='Results')
            # Normalizing intensity to pulse energy (assuming 1mJ on average)
            intensity_normalized = (evt['analysis']['intensity'].data / evt['analysis']['averagePulseEnergy'].data) * 1.0
            add_record(evt['analysis'], 'analysis', 'intensity_normalized', intensity_normalized)

            # Plot Intensity (normalized)
            plotting.line.plotHistory(evt['analysis']['intensity'], history=10000, 
                                      name = 'Intensity normalized (from sizing)', group='Results')

            # Center position
            plotting.correlation.plotMeanMap(evt["analysis"]["cx"], evt["analysis"]["cy"],
                                             intensity_normalized, group='Results',
                                             name='Wavefront (center vs. intensity)', 
                                             xmin=-10, xmax=10, xbins=21,
                                             ymin=-10, ymax=10, ybins=21,
                                             xlabel='Center position in x', ylabel='Center position in y')
            if evt["analysis"]["diameter"].data<33. and evt["analysis"]["diameter"].data>30.:
                corrected_front=front.data.copy()
                corrected_front[corrected_front<aduThreshold]=0
                add_record(evt['analysis'], 'analysis', 'corrected front', corrected_front.astype(numpy.float16))
                plotting.image.plotImage(evt['analysis']['corrected front'], group='Results', name='TBSV hit')


    if hit:
        corrected_front=front.data.copy()
        corrected_front[corrected_front<aduThreshold]=0
        add_record(evt['analysis'], 'analysis', 'corrected front', corrected_front.astype(numpy.float16))
        plotting.image.plotImage(evt['analysis']['corrected front'], mask=mask_front_s, group='Results', name='Rubisco hit')
        

    # --------------------------------------------------------------------------------
    
    # Look at TOF detector
    if do_tof:
        tof_trace = evt[acqiris_type][acqiris_key]
        plotting.line.plotTrace(tof_trace, group="ToF detector")
        
        # The first 1000 points seem to arrive before the x-rays so we'll use
        # the first 900 to normalize the background to 0
        corrected_tof = tof_trace.data-numpy.mean(tof_trace.data[:900])
        add_record(evt['analysis'], 'analysis', 'Corrected ToF', corrected_tof)
        plotting.line.plotTrace(evt["analysis"]['Corrected ToF'], group="ToF detector")

        # ToF intensity score.
        ToF_score = numpy.sum(corrected_tof[:5000] * numpy.linspace(1,0,5000))
        add_record(evt['analysis'], 'analysis', 'ToF score', ToF_score)
        plotting.line.plotHistory(evt["analysis"]['ToF score'], 
                                  name='ToF score', group="ToF detector")

        # ToF score vs KB parameter
        plotting.correlation.plotScatter(evt["analysis"]["ToF score"], evt["parameters"][kb_13_key], 
                                         name='ToF score vs. KB 13', group='KB Optimization', history=100)

        plotting.correlation.plotScatter(evt["analysis"]["ToF score"], evt["parameters"][kb_11_key], 
                                         name='ToF score vs. KB 11', group='KB Optimization', history=100)


        
        
        
    # --------------------------------------------------------------------------------
    # Compute and plot hitrate
    analysis.hitfinding.hitrate(evt, hit, history=25, unit='percent', outkey='hitrate')
    analysis.hitfinding.hitrate(evt, correctsized_hit, history=25, unit='percent', outkey='hitrate_correctsize')
    if main_worker:
         plotting.line.plotHistory(evt["analysis"]["hitrate"], 
                                  name='Hit rate - overall [%]', group="Hit rates")
         plotting.line.plotHistory(evt["analysis"]["hitrate_correctsize"], 
                                  name='Hit rate - correct size [%]', group="Hit rates")

    # -----------------------------------------------------------------------------

    # Plot injector positions
    if data_mode != "simulation":
       plotting.line.plotHistory(evt["parameters"][injector_x_key], group='Injector')
       plotting.line.plotHistory(evt["parameters"][injector_z_key], group='Injector')
    
       # Plot Injector position vs. hitscore
       plotting.correlation.plotScatter(evt["analysis"]["litpixel: hitscore"], evt["parameters"][injector_x_key], 
                                    name='Injector x position vs. hitscore', group="Injector", history=100000)
    
    
       plotting.correlation.plotScatter(evt["analysis"]["litpixel: isHit"], evt["parameters"][injector_x_key], 
                                        name='Injector x position vs. ishit', group="Injector", history=100000)
    
       if hit:
           plotting.line.plotHistory(evt["parameters"][injector_x_key], 
                                   name='Injector position of hits', group='Injector')
           plotting.line.plotHistory(evt["parameters"][injector_z_key], 
                                   name='Z Injector position of hits', group='Injector')
    
       # Plot Injector position vs. hitscore
       if pulse:
           plotting.correlation.plotMeanMap(evt["parameters"][injector_x_key], 
                                            evt["parameters"][injector_z_key], hitscore, group='Injector',
                                            name='hitscoreMeanMap', **hitrateMeanMapParams)
    
       # Plot Injector position vs. hitrate
       if pulse and main_worker:
           plotting.correlation.plotScatter(evt["analysis"]["hitrate"], evt["parameters"][injector_x_key], 
                                            name='Injector x position vs. Hit rate [%]', history=int(1e6),
                                            group="Injector")

       # Plot Injector current
           plotting.correlation.plotScatter(evt["analysis"]["hitrate"], evt["parameters"]["IOC:AMO:HV1:VHS2:CH0:CurrentMeasure"],
                                            name='Injector current vs. hitrate', group="Injector", history=100000)

           plotting.line.plotHistory(evt["parameters"]["IOC:AMO:HV1:VHS2:CH0:CurrentMeasure"],
                                            name='Injector current', group="Injector", history=100000)

        
    # Plot MeanMap of hitrate(y,z)
    #if pulse and main_worker:
    #    plotting.correlation.plotMeanMap(evt["parameters"][injector_x_key], evt["parameters"][injector_z_key], 
    #                                     evt["analysis"]["hitrate"].data, group="Injector", 
    #                                     name='hitrateMeanMap', **hitrateMeanMapParams)

    # -------------------------------------------------------------------------------------

    # Plot Injector/Nozzle pressures
    #plotting.line.plotHistory(evt["parameters"][nozzle_pressure_1_key], group='Injector')
    #plotting.line.plotHistory(evt["parameters"][nozzle_pressure_2_key], group='Injector')

    #plotting.line.plotHistory(evt["parameters"][injector_pressure_1_key], group='Injector')
    #plotting.line.plotHistory(evt["parameters"][injector_pressure_2_key], group='Injector')
    #plotting.line.plotHistory(evt["parameters"][injector_pressure_3_key], group='Injector')

    # -------------------------------------------------------------------------------------

    # Plot KB settings
    #plotting.line.plotHistory(evt["parameters"][kb_04_key], group='KB settings')
    #plotting.line.plotHistory(evt["parameters"][kb_05_key], group='KB settings')
    #plotting.line.plotHistory(evt["parameters"][kb_06_key], group='KB settings')
    #plotting.line.plotHistory(evt["parameters"][kb_07_key], group='KB settings')
    #plotting.line.plotHistory(evt["parameters"][kb_08_key], group='KB settings')
    #plotting.line.plotHistory(evt["parameters"][kb_09_key], group='KB settings')
    #plotting.line.plotHistory(evt["parameters"][kb_10_key], group='KB settings')
    plotting.line.plotHistory(evt["parameters"][kb_11_key], group='KB settings')
    plotting.line.plotHistory(evt["parameters"][kb_12_key], group='KB settings')
    plotting.line.plotHistory(evt["parameters"][kb_13_key], group='KB settings')
    plotting.line.plotHistory(evt["parameters"][kb_14_key], group='KB settings')

    # -----------------------------------------------------------------------

#    if data_mode != "simulation":
#        # Look at TOF detector
#        tof_trace = evt[acqiris_type][acqiris_key]
#        plotting.line.plotTrace(tof_trace, group="ToF detector")
#        plotting.correlation.plotScatter(evt["parameters"][injector_x_key], numpy.sum(evt[acqiris_type][acqiris_key]),
#                                         name='Injector x position vs. tof sum', history=int(1e6),
#                                         group="Injector")


    # -------------------------------------------------------------------------------------
    
    # Collecting background
    #if do_stacks and pulse and hit:
    #    # Update
    #    bg_front.add(front_no_geometry.data)
    #    # Reduce stack
    #    bg_front.reduce()
    #    # Write to file
    #    bg_front.write(evt, directory=bg_dir, verbose=True)
        
    # ------------------------------------------------------------------------------------ 

    # Gijs Analysis
    #glorious_hit = False
    #if hit and i_like_flamingos:
    #    glorious_hit = redflamingo.get_data(evt[front_type_s][front_key_s].data, hitcounter) == 1
    #if glorious_hit and False:
    #    print "Glory"
    #    plotting.image.plotImage(evt[front_type_s][front_key_s], group="Gijs", 
    #                             msg='', name="pnCCD front (glorious hit)")#, mask=mask_front) 

    # ---------------------------------------------------------------------------------------

    # Saving to file
    if do_save and hit:
        D = {}
        D['front']  = evt[front_type_s][front_key_s].data
        D['hitscore']  = evt["analysis"]["litpixel: hitscore"].data  
        D['timestamp'] = evt["eventID"]["Timestamp"].timestamp
        D['fiducial']  = evt["eventID"]["Timestamp"].fiducials
        W.write_slice(D)

    # ---------------------------------------------------------

def end_of_run():
    if do_save:
        tmpfile  = w_dir + "/r%04d.h5" %conf.run_nr
        donefile = w_dir + "/r%04d_done.h5" %conf.run_nr
        os.system('mv %s %s' %(tmpfile, donefile))
