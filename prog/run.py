import subprocess, os

INTERSECTION_CONES_TOP_DIRECTORY = os.path.expanduser("~/Dropbox/IntersectionCones")

GAP_FILENAME             = os.path.expanduser("~/Gap/gap4r8/bin/gap.sh")
SAGE_FILENAME            = os.path.expanduser("~/Sage/SageMath/sage")
GAP_PROGRAM_FILENAME     = os.path.expanduser("%s/prog/prog.g" % INTERSECTION_CONES_TOP_DIRECTORY)
SAGE_PROGRAM_FILENAME    = os.path.expanduser("%s/prog/sageprog.sage" % INTERSECTION_CONES_TOP_DIRECTORY)

MK_RC_DIRECTORY          = True # Change to false when directory exists
RAW_CONES_DIRECTORY      = os.path.expanduser("%s/raw_cones/" % INTERSECTION_CONES_TOP_DIRECTORY)

MK_TC_DIRECTORY          = True # Change to false when directory exists
TRUE_CONES_DIRECTORY     = os.path.expanduser("%s/true_cones/" % INTERSECTION_CONES_TOP_DIRECTORY)

MK_BC_DIRECTORY          = True # Change to false when directory exists
BLOCKED_CONES_DIRECTORY  = os.path.expanduser("%s/blocked_cones/" % INTERSECTION_CONES_TOP_DIRECTORY)

MK_IC_DIRECTORY          = True # Change to false when directory exists
INT_CONES_DIRECTORY      = os.path.expanduser("%s/intersected_cones/" % INTERSECTION_CONES_TOP_DIRECTORY)





subprocess.run(["sh", GAP_FILENAME, GAP_PROGRAM_FILENAME])
subprocess.run([SAGE_FILENAME])


# Sage part:

# load("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/prog/sageprog.sage")

