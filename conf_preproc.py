import os,sys

this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(this_dir)

import time, numpy
import scipy.constants as constants
import logging
import h5py
import csv

import analysis.event
import analysis.pixel_detector
import analysis.beamline
import analysis.hitfinding
import analysis.patterson
import utils.reader
import utils.cxiwriter
utils.cxiwriter.h5writer.logger.setLevel("INFO")
import utils.array
import ipc.mpi
from backend import add_record

from quad_correction import fix_upper_right_quadrant

# Commandline arguments
from utils.cmdline_args import argparser, add_config_file_argument
add_config_file_argument('--out-dir', metavar='out_dir', nargs='?',
                         help="Output directory", required=True,
                         type=str)
add_config_file_argument('--hitscore-threshold', metavar='hitscore_threshold',
                         help="Hitscore threshold [if not provided read from CSV file]",
                         type=int)
add_config_file_argument('--multiscore-threshold', metavar='multiscore_threshold',
                         help="Multiscore threshold [if not provided read from CSV file]",
                         type=int)
add_config_file_argument('--photon-energy-ev', metavar='photon_energy_ev',
                         help="Manually set nominal photon energy in unit electron volts (used for example for pnCCD gain calculation and hit finding) [if not provided read from CSV file]",
                         type=float)
add_config_file_argument('--rescaling-factors-asics', metavar='rescaling_factors_asics',
                         help="Manually set the 4 rescaling factors for upper right quadrant (for example 2.0,1.0,1.0,1.0) [if not provided read from CSV file]",
                         type=str)
add_config_file_argument('--output-level',
                         help='Output level defines how much data per event will be stored (default=3, min=0, max=3)', 
                         type=int, default=3)
add_config_file_argument('--do-raw',
                         help="Only do pedestal correction, no other corrections",
                         type=int, default=0)
add_config_file_argument('--do-quad',
                         help="Correct artefacts in \"upper right\" quadrant",
                         type=int, default=1)
add_config_file_argument('--do-cmc',
                         help="Apply common mode correction along \"horizontal\" lines",
                         type=int, default=1)
add_config_file_argument('--do-metrology',
                         help="Move pixels to their physical locations",
                         type=int, default=1)


args = argparser.parse_args()

# The way how the output level argument is interpreted
save_anything = args.output_level > 0
save_tof = args.output_level >= 2
save_pnccd = args.output_level >= 3

# Skipping events with no FEL pulse
do_skipdark = True

# ------------------------------
# PSANA
# ------------------------------
import conf_amol3116 as conf
state = conf.state

# Switch between different datasets
data_mode = "amol3116"

# PNCCD Detector
front_type = conf.pnccd_type
front_key  = conf.pnccd_key

# Detector position (default values)
dx_front = conf.dx_front
dy_front = conf.dy_front
dy_gap   = conf.dy_gap

# Trigger for cluster source
cluster_trigger_key = conf.cluster_trigger_key
cluster_trigger_event_code_key = conf.cluster_trigger_event_code_key

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
        
# Tof detector
acqiris_type = conf.acqiris_type
acqiris_key  = conf.acqiris_key
        
if data_mode == "amol3116":
    # Trigger for cluster source
    cluster_trigger_key = conf.cluster_trigger_key
    # MCP Detector
    back_type = conf.mcp_type
    back_key  = conf.mcp_key

# ------------------------------
# PSANA
# ------------------------------
if data_mode == "amol3116":
    import conf_amol3116 as conf
    conf.state['LCLS/PsanaConf'] = this_dir + '/psana_cfg/pnccd_front.cfg'
else:
    raise Exception("%s is not a valid data_mode." % data_mode)
    sys.exit(1)
state = conf.state

# For CXI writer we reduce the number of event readers by one
state['reduce_nr_event_readers'] = 1

# Gain
gain = None
gain_mode = None
gain_mode_key = conf.gain_mode_key

# Thresholds for hit finding and nominal photon energy for run
# If these parameters are not passed from the command line set to tabulated value for this run
hitscore_threshold = None
multiscore_threshold = None
nominal_photon_energy_eV = None
rescaling_factors_asics = None
run_type = "Unknown"
filename_csv = '%s_run_params.csv' % (data_mode)
with open(filename_csv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if numpy.int16(row["RunNr"]) == conf.run_nr:
            hitscore_threshold = numpy.int16(row["HitscoreThreshold"])
            multiscore_threshold = numpy.int16(row["MultiscoreThreshold"])
            nominal_photon_energy_eV = numpy.float64(row["PhotonEnergyEV"])
            rescaling_factors_asics = numpy.array([numpy.float64(row["PnCCDRescaleQ1AC0"]),
                                                   numpy.float64(row["PnCCDRescaleQ1AC1"]),
                                                   numpy.float64(row["PnCCDRescaleQ1AC2"]),
                                                   numpy.float64(row["PnCCDRescaleQ1AC3"])], dtype='float64')
            dx_front = numpy.int16(row["PnCCDDX"])
            dy_front = numpy.int16(row["PnCCDDY"])
            dy_gap   = numpy.int16(row["PnCCDYGAP"])
            run_type = row["RunType"]

if args.hitscore_threshold is not None:
    hitscore_threshold = args.hitscore_threshold
if hitscore_threshold is None:
    raise Exception("Could not find tabulated values for the hitscore threshold of this run (%i) in %s. Please specify the hitscore manually (--hitscore-threshold) or pre-process another run." % (conf.run_nr, filename_csv))
    sys.exit(1)

if args.multiscore_threshold is not None:
    multiscore_threshold = args.multiscore_threshold
if multiscore_threshold is None:
    raise Exception("Could not find tabulated values for the multiscore threshold of this run (%i) in %s. Please specify the multiscore manually (--multiscore-threshold) or pre-process another run." % (conf.run_nr, filename_csv))
    sys.exit(1)

if args.photon_energy_ev is not None:
    nominal_photon_energy_eV = args.photon_energy_ev
if nominal_photon_energy_eV is None:
    raise Exception("Could not find tabulated values for the nominal photon energy of this run (%i) in %s. Please specify the photon energy manually (--photon-energy-ev) or pre-process another run." % (conf.run_nr, filename_csv))
    sys.exit(1)

if args.rescaling_factors_asics is not None:
    rescaling_factors_asics = numpy.array([numpy.float64(f) for f in args.rescaling_factors_asics.split(',')], dtype='float64')
if rescaling_factors_asics is None:
    raise Exception("Could not find tabulated values for the asics rescaling factors of this run (%i) in %s. Please specify these values manually (--rescaling-factors-asics) or pre-process another run." % (conf.run_nr, filename_csv))
    sys.exit(1)

# Patterson
# ---------
patterson_threshold = 5.
patterson_floor_cut = 50.
patterson_mask_smooth = 5.
patterson_diameter = 60.

# Masks
# -----
mask_front_templ = utils.reader.MaskReader(this_dir + "/mask/mask_front_production.h5","/data/data").boolean_mask
(ny_front,nx_front) = mask_front_templ.shape
mask_front_ext_templ   = utils.reader.MaskReader(this_dir + "/mask/mask_front_ext4.h5","/data/data").boolean_mask

event_number = 0
hitcounter = 0

# Cluster source trigger
# ----------------------
cluster_trigger_treshold = 0

W = None
out_dir = args.out_dir
tmpfile  = '%s/.%s_r%04d_ol%i.h5' % (out_dir, data_mode, conf.run_nr, args.output_level)
donefile = '%s/%s_r%04d_ol%i.h5' % (out_dir, data_mode, conf.run_nr, args.output_level)
D_solo = {}

def beginning_of_run():
    global W
    W = utils.cxiwriter.CXIWriter(tmpfile, chunksize=10, compression=None)
    
# ---------------------------------------------------------
# E V E N T   C A L L
# ---------------------------------------------------------

def onEvent(evt):
    global event_number
    global hitcounter
    global gain
    global gain_mode

    ft = front_type
    fk = front_key
  
    event_number += 1

    # ------------------- #
    # INITIAL DIAGNOSTICS #
    # ------------------- #

    # Time measurement
    analysis.event.printProcessingRate()

    # Average pulse energies, in case there is no pulse energy it returns None
    analysis.beamline.averagePulseEnergy(evt, evt["pulseEnergies"]) 
    # Skip frames that do not have a pulse energy
    try:
        pulse_energy = evt['analysis']['averagePulseEnergy']
        pulse = True
    except (TypeError,KeyError):
        pulse = False
    if do_skipdark and not pulse:
        print "No pulse energy. Skipping event."
        return

    # Average pulse energies, in case there is no pulse energy it returns None 
    analysis.beamline.averagePhotonEnergy(evt, evt["photonEnergies"])
    # Skip frames that do not have a photon energy
    try:
        evt['analysis']['averagePhotonEnergy']
    except (TypeError,KeyError):
        print "No photon energy. Skipping event."
        return       

    # Set to nominal photon energy
    photon_energy_ev_avg = evt['analysis']['averagePhotonEnergy']
    photon_energy_ev_nom = add_record(evt["analysis"], "analysis", "nominalPhotonEnergy" , nominal_photon_energy_eV)
    
    # Gain of PNCCD
    try:
        gain_mode = evt['parameters'][gain_mode_key].data
    except:
        print "Could not read gain mode from data stream. Skipping event."
        return
    analysis.pixel_detector.pnccdGain(evt, photon_energy_ev_nom, gain_mode)
    gain = evt['analysis']['gain'].data

    # Skip frames that do not have the pnCCDs
    try:
        evt[ft][fk]
    except (TypeError,KeyError):
        print "No front pnCCD. Skipping event."
        return

    front = evt[ft][fk]
    mask_front = mask_front_templ.copy()
    mask_front_ext = mask_front_ext_templ.copy()

    front_is_read_only = True

    if len(front.data.shape) == 3:
        tmp = numpy.zeros(shape=(1024, 1024), dtype=front.data.dtype)
        # Left
        tmp[:512,:512] = front.data[0,:,:]
        tmp[512:,:512] = front.data[1,::-1,::-1]
        # Right
        tmp[:512,512:] = front.data[3,:,:]
        tmp[512:,512:] = front.data[2,::-1,::-1]
        front_is_read_only = False
    else:
        tmp = front.data
    ft = "analysis"
    fk = "preprocessed - " + fk
    add_record(evt[ft], ft, fk , tmp)
    front = evt[ft][fk]

    if not args.do_raw:        
        if args.do_quad:
            if front_is_read_only:
                front.data = front.data.copy()
            # Compensate for artefacts in "upper right" quadrant
            fix_upper_right_quadrant(img_quad=front.data[:512, 512:], msk_quad=mask_front[:512, 512:], rescaling_factors_asics=rescaling_factors_asics)

        # Common mode correction for PNCCD
        if args.do_cmc and not args.do_raw:
            analysis.pixel_detector.commonModePNCCD2(evt, ft, fk, signal_threshold=10, min_nr_pixels_per_median=100)
            ft = "analysis"
            fk  = "corrected - " + fk
            front = evt[ft][fk]
            
        if args.do_metrology:
            # Fix wrong wiring: Move "upper" quadrants by one pixel inwards
            # "Top left"
            #front.data[:512,1:512] = front.data[:512,:512-1]
            #mask_front[:512,0] = False
            # "Top right"
            #front.data[:512,512:-1] = front.data[:512,512+1:]
            #mask_front[:512,-1] = False
    
            # Move right detector half (assemble)
            front = analysis.pixel_detector.moveHalf(evt, front, horizontal=dx_front, vertical=dy_front, outkey='data_half-moved')
            ft = "analysis"
            fk  = "data_half-moved"
            mask_front = analysis.pixel_detector.moveHalf(evt, add_record(evt["analysis"], "analysis", "mask",mask_front), horizontal=dx_front, vertical=dy_front, outkey='mask_half-moved').data

            mask_front_ext = analysis.pixel_detector.moveHalf(evt, add_record(evt["analysis"], "analysis", "mask_ext", mask_front_ext), horizontal=dx_front, vertical=dy_front, outkey='mask_half-moved_ext' ).data
            d = 35
            mask_front_ext[512-d:512+d,:] = 0

            # Add vertical gap
            gap = 6
            front_cpy = numpy.zeros((front.data.shape[0]+dy_gap, front.data.shape[1]))
            front_cpy[:512-dy_front,:512] = front.data[:512-dy_front,:512]
            front_cpy[512-dy_front+dy_gap:,:512] = front.data[512-dy_front:,:512]
            front_cpy[:512,512+dx_front:] = front.data[:512,512+dx_front:]
            front_cpy[512+dy_gap:,512+dx_front:] = front.data[512:,512+dx_front:]
            front = add_record(evt[ft], ft,fk,front_cpy)

            mask_front_cpy = numpy.zeros((mask_front.shape[0]+dy_gap, mask_front.shape[1]), dtype=mask_front.dtype)
            mask_front_cpy[:512-dy_front,:512] = mask_front[:512-dy_front,:512]
            mask_front_cpy[512-dy_front+dy_gap:,:512] = mask_front[512-dy_front:,:512]
            mask_front_cpy[:512,512+dx_front:] = mask_front[:512,512+dx_front:]
            mask_front_cpy[512+dy_gap:,512+dx_front:] = mask_front[512:,512+dx_front:]
            mask_front = mask_front_cpy

            mask_front_ext_cpy = numpy.zeros((mask_front_ext.shape[0]+dy_gap, mask_front_ext.shape[1]), dtype=mask_front_ext.dtype)
            mask_front_ext_cpy[:512-dy_front,:512] = mask_front_ext[:512-dy_front,:512]
            mask_front_ext_cpy[512-dy_front+dy_gap:,:512] = mask_front_ext[512-dy_front:,:512]
            mask_front_ext_cpy[:512,512+dx_front:] = mask_front_ext[:512,512+dx_front:]
            mask_front_ext_cpy[512+dy_gap:,512+dx_front:] = mask_front_ext[512:,512+dx_front:]
            mask_front_ext = mask_front_cpy

    # Finding hits
    analysis.hitfinding.countLitPixels(evt, front, aduThreshold=gain, 
                                       hitscoreThreshold=hitscore_threshold, mask=mask_front)
    hit = evt["analysis"]["litpixel: isHit"].data
    if hit: hitcounter += 1

    # Finding multiple hits
    if hit:
        analysis.patterson.patterson(evt, ft, fk, mask_front_ext, 
                                     floor_cut=patterson_floor_cut, mask_smooth=patterson_mask_smooth,
                                     threshold=patterson_threshold ,diameter_pix=patterson_diameter,
                                     crop=512, full_output=True)
        multiple_hit = evt["analysis"]["multiple score"].data > multiscore_threshold

    # Check if cluster source is turned on/off
    try:
        event_code_number = evt['parameters'][cluster_trigger_event_code_key].data
    except KeyError:
        print "WARNING: %s cannot be found in parameters." % cluster_trigger_event_code_key
        event_code_number = None
    if event_code_number is None or event_code_number in evt['eventCodes']['EvrEventCodes'].data:
        cluster_source_on = True
    else:
        cluster_source_on = False

    # Saving to file
    if hit and save_anything:
        D = {}

        D["entry_1"] = {}

        if save_pnccd:
            D["entry_1"]["detector_1"] = {}
            D["entry_1"]["detector_2"] = {}

        D["entry_1"]["detector_3"] = {}
        D["entry_1"]["event"] = {}
        D["entry_1"]["injector"] = {}
        D["entry_1"]["FEL"] = {}
        D["entry_1"]["result_1"] = {}

        if save_pnccd:
            # PNCCD
            D["entry_1"]["detector_1"]["data"] = numpy.asarray(evt[ft][fk].data, dtype='float16')
            if ipc.mpi.is_main_event_reader() and len(D_solo) == 0:
                bitmask = numpy.array(mask_front, dtype='uint16')
                bitmask[bitmask==0] = 512
                bitmask[bitmask==1] = 0
                D_solo["entry_1"] = {}
                D_solo["entry_1"]["detector_1"] = {}
                D_solo["entry_1"]["detector_1"]["mask"] = bitmask
            D["entry_1"]["detector_1"]["gain"] = gain
            D["entry_1"]["detector_1"]["gain_mode"] = gain_mode
            D["entry_1"]["detector_1"]["rescaling_factors_asics"] = rescaling_factors_asics
            D["entry_1"]["detector_1"]["patterson"] = numpy.asarray(evt["analysis"]["patterson"].data, dtype='float16')

        if save_tof:
            # TOF
            D["entry_1"]["detector_2"]["TOF"] = evt[acqiris_type][acqiris_key].data
        
        # GMD
        D["entry_1"]["detector_3"]["pulse_energy_mJ"] = pulse_energy.data
        
        # HIT PARAMETERS
        D["entry_1"]["result_1"]["hitscore_litpixel"] = evt["analysis"]["litpixel: hitscore"].data
        D["entry_1"]["result_1"]["hitscore_litpixel_threshold"] = hitscore_threshold
        D["entry_1"]["result_1"]["multiscore_patterson"] = evt["analysis"]["multiple score"].data
        D["entry_1"]["result_1"]["multiscore_patterson_threshold"] = multiscore_threshold

        # EVENT IDENTIFIERS
        D["entry_1"]["event"]["timestamp"] = evt["eventID"]["Timestamp"].timestamp
        D["entry_1"]["event"]["timestamp2"] =  evt["eventID"]["Timestamp"].timestamp2
        D["entry_1"]["event"]["fiducial"] = evt["eventID"]["Timestamp"].fiducials
        D["entry_1"]["event"]["clusteron"] = cluster_source_on

        # INJECTOR
        D["entry_1"]["injector"]["injector_position_x"] = evt["parameters"][injector_x_key].data
        D["entry_1"]["injector"]["injector_position_y"] = evt["parameters"][injector_z_key].data
        D["entry_1"]["injector"]["injector_pressures_1"] = evt["parameters"][injector_pressure_1_key].data
        D["entry_1"]["injector"]["injector_pressures_2"] = evt["parameters"][injector_pressure_2_key].data
        D["entry_1"]["injector"]["nozzle_pressures_1"] = evt["parameters"][nozzle_pressure_1_key].data
        D["entry_1"]["injector"]["nozzle_pressures_2"] = evt["parameters"][nozzle_pressure_2_key].data

        # FEL
        D["entry_1"]["FEL"]["photon_energy_eV_nominal"] = photon_energy_ev_nom.data
        D["entry_1"]["FEL"]["wavelength_nm_nominal"] = constants.c*constants.h/(photon_energy_ev_nom.data*constants.e)/1E-9
        D["entry_1"]["FEL"]["photon_energy_eV_SLAC"] = photon_energy_ev_avg.data
        D["entry_1"]["FEL"]["wavelength_nm_SLAC"] = constants.c*constants.h/(photon_energy_ev_avg.data*constants.e)/1E-9

        W.write_slice(D)

def end_of_run():
    if ipc.mpi.is_main_event_reader():
        if "entry_1" not in D_solo:
            D_solo["entry_1"] = {}
        D_solo["entry_1"]["run_type"] = run_type
        W.write_solo(D_solo)

    W.close(barrier=True)

    if ipc.mpi.is_main_event_reader():
        if save_anything:
            with h5py.File(tmpfile, "r+") as f:
                f["/entry_1/detector_3/data"] = h5py.SoftLink('/entry_1/detector_3/pulse_energy_mJ')
                if save_pnccd:
                    f["/entry_1/data_1"] = h5py.SoftLink('/entry_1/detector_1')
                if save_tof:
                    f["/entry_1/data_2"] = h5py.SoftLink('/entry_1/detector_2')
                    f["/entry_1/detector_2/data"] = h5py.SoftLink('/entry_1/detector_2/TOF')
                print "Successfully wrote soft links."
        os.system('mv %s %s' %(tmpfile, donefile))
        print "Moved temporary file %s to %s." % (tmpfile, donefile)

    print "Clean exit."
