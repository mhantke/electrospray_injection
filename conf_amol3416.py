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

# PSANA
state = {}
state['Facility'] = 'LCLS'

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
    state['LCLS/DataSource'] = 'exp=amol3416:run=%d' %run_nr

state['LCLS/PsanaConf'] = this_dir + '/psana_cfg/pnccd_front.cfg'

# Read events from file
idxfile = None

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

# INJECTOR VOLTAGE
#-----------------
injector_voltage_key = "IOC:AMO:HV1:VHS2:CH0:VoltageMeasure"
injector_current_key = "IOC:AMO:HV1:VHS2:CH0:CurrentMeasure"

# ELECTROSPRAY
#----------------
electrospray_key = "AMO:PPL:MMS:18.RBV"

# PULSE LENGTH
pulselength_key = "SIOC:SYS0:ML00:AO820"

# PNCCD DETECTOR
# --------------
pnccd_type = "image"
pnccd_key  = "pnccdFront[%s]" % pnccd_type
gain_mode_key = "PNCCD:FRONT:GAIN"

dx_front  = 73
dy_front  = -3

# PNCCD STAGE
# -----------
pnccd_motorx = "AMO:LMP:MMS:09.RBV"
pnccd_motory_top = "AMO:LMP:MMS:07.RBV"
pnccd_motory_bottom = "AMO:LMP:MMS:08.RBV"
pnccd_motorz = "AMO:LMP:MMS:10.RBV"

# TOF DETECTOR
# ------------
acqiris_type = "ionTOFs"
acqiris_key  = "Acqiris 2 Channel 0"

# MCP DETECTOR
# ------------
mcp_type = "camera"
mcp_key = "mcp"

# KB MIRRORS
# ----------
kb_04_key = "AMO:KBO:MMS:04.RBV"
kb_05_key = "AMO:KBO:MMS:05.RBV"
kb_06_key = "AMO:KBO:MMS:06.RBV"
kb_07_key = "AMO:KBO:MMS:07.RBV"
kb_08_key = "AMO:KBO:MMS:08.RBV"
kb_09_key = "AMO:KBO:MMS:09.RBV"
kb_10_key = "AMO:KBO:MMS:10.RBV"
kb_11_key = "AMO:KBO:MMS:11.RBV"
kb_12_key = "AMO:KBO:MMS:12.RBV"
kb_13_key = "AMO:KBO:MMS:13.RBV"
kb_14_key = "AMO:KBO:MMS:14.RBV"

