# imports
import os, sys

# Commandgline arguments
from hummingbird import parse_cmdline_args
args = parse_cmdline_args()

# Path to current directory
this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(this_dir)

# Run online
do_online = False
do_autoonline = True

# --------------------------------------
# P S A N A (AMOC6914) - Use for testing
# --------------------------------------
state = {}
state['Facility'] = 'LCLS'

idxfile = None
amoopr = False
if do_autoonline:
    import getpass
    if getpass.getuser() == "amoopr":
        do_online=True
        amoopr=True

if do_online:    
    state['LCLS/DataSource'] = 'shmem=psana.1:stop=no'
    do_save = False
    do_stack = False
else:
    run_nr = args.lcls_run_number
    if run_nr is None:
        run_nr = 2
    state['LCLS/DataSource'] = 'exp=amol3116:run=%d' %run_nr

state['LCLS/PsanaConf'] = this_dir + '/psana_cfg/pnccd_front.cfg'

if idxfile is not None:
    import h5py
    idxdata = h5py.File(idxfile, "r")
    fiducials = idxdata['/fiducial'][()]
    times = idxdata['/timestamp2'][()]
    len = fiducials.shape[0]
    rankrange = numpy.arange(ipc.mpi.slave_rank(), len, ipc.mpi.nr_workers())
    state['fiducials'] = fiducials[rankrange]
    state['times'] = times[rankrange]

# INJECTOR MOTORS
# ---------------
injector_x_key = "AMO:PPL:MMS:01.RBV"
injector_z_key = "AMO:PPL:MMS:03.RBV"

# NOZZLE PRESSURES
#-----------------
nozzle_pressure_1_key = "AMO:SDS:PAREG:02:PRESS"
nozzle_pressure_2_key = "AMO:SDS:PAREG:03:PRESS"

# INJECTOR PRESSURES
#-------------------
injector_pressure_1_key = "AMO:LMP:VG:43:PRESS"
injector_pressure_2_key = "AMO:LMP:VG:40:PRESS"
injector_pressure_3_key = ""

# TRIGGER FOR CLUSTER SOURCE
# --------------------------
cluster_trigger_key = "AMO:EVR:01:TRIG0:CNT"
cluster_trigger_event_code_key = 'AMO:EVR:01:TRIG0:EC_RBV'

# THERMO COUPLE
# -------------
thermo_couple_1_key = "AMO:PPL:STC:05"
thermo_couple_2_key = "AMO:PPL:STC:06"

# PNCCD DETECTOR
# --------------
pnccd_type = "image"
pnccd_key  = "pnccdFront[%s]" % pnccd_type
gain_mode_key = "PNCCD:FRONT:GAIN"

dx_front = 42
dy_front = -3
dy_gap   = 0

# TOF DETECTOR
# ------------
acqiris_type = "ionTOFs"
acqiris_key  = "Acqiris 2 Channel 0"

# MCP DETECTOR
mcp_type = "camera"
mcp_key = "mcp"
