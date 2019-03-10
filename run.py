import subprocess, os

############ HEADER --- these attributes may require modifying. ############

INTERSECTION_CONES_TOP_DIRECTORY = os.path.expanduser(
	"~/Dropbox/IntersectionCones") # Exclude trailing slash.

GAP_FILENAME             = os.path.expanduser("~/Gap/gap-4.9.1/bin/gap.sh")
SAGE_FILENAME            = os.path.expanduser("~/Sage/SageMath/sage")
# Where you see "%s", read the top directory in its place.
GAP_PROGRAM_FILENAME     = os.path.expanduser("%s/prog/prog.g"
	% INTERSECTION_CONES_TOP_DIRECTORY)
SAGE_PROGRAM_FILENAME    = os.path.expanduser("%s/prog/sageprog.sage"
	% INTERSECTION_CONES_TOP_DIRECTORY)

MK_RC_DIRECTORY          = True # Change to false when directory exists
RAW_CONES_DIRECTORY      = os.path.expanduser("%s/raw_cones/"
	% INTERSECTION_CONES_TOP_DIRECTORY)

MK_TC_DIRECTORY          = True # Change to false when directory exists
TRUE_CONES_DIRECTORY     = os.path.expanduser("%s/true_cones/"
	% INTERSECTION_CONES_TOP_DIRECTORY)

MK_BC_DIRECTORY          = True # Change to false when directory exists
BLOCKED_CONES_DIRECTORY  = os.path.expanduser("%s/blocked_cones/"
	% INTERSECTION_CONES_TOP_DIRECTORY)

MK_IC_DIRECTORY          = True # Change to false when directory exists
INT_CONES_DIRECTORY      = os.path.expanduser("%s/intersected_cones/"
	% INTERSECTION_CONES_TOP_DIRECTORY)

############                        END HEADER                  ############



# I'm sure there's a better way to communicate these settings forward,
# but this is what I came up with.
f = open("%s/settings.txt" % INTERSECTION_CONES_TOP_DIRECTORY, "w")
f.write(
	"""return rec(
GAP_FILENAME             := "%s",
SAGE_FILENAME            := "%s",

GAP_PROGRAM_FILENAME     := "%s",
SAGE_PROGRAM_FILENAME    := "%s",  

MK_RC_DIRECTORY          := %s, 
RAW_CONES_DIRECTORY      := "%s", 

MK_TC_DIRECTORY          := %s, 
TRUE_CONES_DIRECTORY     := "%s", 

MK_BC_DIRECTORY          := %s, 
BLOCKED_CONES_DIRECTORY  := "%s", 

MK_IC_DIRECTORY          := %s, 
INT_CONES_DIRECTORY      := "%s", 
);
""" % (GAP_FILENAME, SAGE_FILENAME, GAP_PROGRAM_FILENAME, SAGE_PROGRAM_FILENAME, MK_RC_DIRECTORY,
RAW_CONES_DIRECTORY, MK_TC_DIRECTORY, TRUE_CONES_DIRECTORY, MK_BC_DIRECTORY, BLOCKED_CONES_DIRECTORY,
MK_IC_DIRECTORY, INT_CONES_DIRECTORY))
f.flush()
f.close()

subprocess.run(["sh", GAP_FILENAME, "-l", ";%s"%INTERSECTION_CONES_TOP_DIRECTORY, GAP_PROGRAM_FILENAME])
subprocess.run([SAGE_FILENAME, SAGE_PROGRAM_FILENAME])


# Sage part:

# load("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/prog/sageprog.sage")

