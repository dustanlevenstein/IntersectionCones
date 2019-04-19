############ HEADER --- these attributes may require modifying. ############

REWRITE_SETTINGS         = True  # Change to False if you want to keep
                                 # settings.txt as is.

INTERSECTION_CONES_TOP_DIRECTORY = \
	"~/Dropbox/IntersectionCones" # Exclude trailing slash.

GAP_FILENAME             = "~/Gap/gap-4.9.1/bin/gap.sh"
SAGE_FILENAME            = "~/Sage/SageMath/sage"
                           # Where you see "%s", read the top directory in
                           # its place.
GAP_PROGRAM_FILENAME     = "%s/prog/prog.g"
SAGE_PROGRAM_FILENAME    = "%s/prog/sageprog.sage"

MK_RC_DIRECTORY          = False # Change to False when directory exists
RAW_CONES_DIRECTORY      = "%s/raw_cones/"
FILL_RAW_CONES           = True 
RC_MIN_PRIME             = 2
RC_MAX_PRIME             = 5
RC_MIN_SYMMETRIC_GP_N    = 1
RC_MAX_SYMMETRIC_GP_N    = 6

MK_TC_DIRECTORY          = False # Change to False when directory exists
TRUE_CONES_DIRECTORY     = "%s/true_cones/"
FILL_TRUE_CONES          = False
TC_MIN_PRIME             = 2
TC_MAX_PRIME             = 5
TC_MIN_SYMMETRIC_GP_N    = 1
TC_MAX_SYMMETRIC_GP_N    = 20

MK_BC_DIRECTORY          = False # Change to False when directory exists
BLOCKED_CONES_DIRECTORY  = "%s/blocked_cones/"
FILL_BC_DIRECTORY        = False # This will fill in blocked cones for every
                                 # file which exists in the raw cones
                                 # directory.
OVERWRITE_OLD_BCS        = False # Overwrite already existing blocked cones.

MK_IC_DIRECTORY          = False # Change to False when directory exists
INT_CONES_DIRECTORY      = "%s/intersected_cones/"
FILL_IC_DIRECTORY        = False # This will fill in the intersection of
                                 # every blocked cone which exists in the
                                 # blocked cones directory.
OVERWRITE_OLD_ICS        = False # Overwrite already existing intersections
                                 # of cones.

MK_DC_DIRECTORY          = True  # Change to False when directory exists
DUAL_CONES_DIRECTORY     = "%s/dual_cones/"

############    END HEADER - EDIT FOLLOWING AT YOUR OWN RISK    ############

import subprocess, os

if REWRITE_SETTINGS:
    INTERSECTION_CONES_TOP_DIRECTORY = os.path.expanduser(INTERSECTION_CONES_TOP_DIRECTORY)

    GAP_FILENAME             = os.path.expanduser(GAP_FILENAME )
    SAGE_FILENAME            = os.path.expanduser(SAGE_FILENAME)

    GAP_PROGRAM_FILENAME     = os.path.expanduser(GAP_PROGRAM_FILENAME     % INTERSECTION_CONES_TOP_DIRECTORY)
    SAGE_PROGRAM_FILENAME    = os.path.expanduser(SAGE_PROGRAM_FILENAME    % INTERSECTION_CONES_TOP_DIRECTORY)
    RAW_CONES_DIRECTORY      = os.path.expanduser(RAW_CONES_DIRECTORY      % INTERSECTION_CONES_TOP_DIRECTORY)
    TRUE_CONES_DIRECTORY     = os.path.expanduser(TRUE_CONES_DIRECTORY     % INTERSECTION_CONES_TOP_DIRECTORY)
    BLOCKED_CONES_DIRECTORY  = os.path.expanduser(BLOCKED_CONES_DIRECTORY  % INTERSECTION_CONES_TOP_DIRECTORY)
    INT_CONES_DIRECTORY      = os.path.expanduser(INT_CONES_DIRECTORY      % INTERSECTION_CONES_TOP_DIRECTORY)

    if MK_RC_DIRECTORY:
        subprocess.run(["mkdir", RAW_CONES_DIRECTORY])
    if MK_TC_DIRECTORY:
        subprocess.run(["mkdir", TRUE_CONES_DIRECTORY])
    if MK_BC_DIRECTORY:
        subprocess.run(["mkdir", BLOCKED_CONES_DIRECTORY])
    if MK_IC_DIRECTORY:
        subprocess.run(["mkdir", INT_CONES_DIRECTORY])
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
FILL_RAW_CONES           := %s,
RC_MIN_PRIME             := %s,
RC_MAX_PRIME             := %s,
RC_MIN_SYMMETRIC_GP_N    := %s,
RC_MAX_SYMMETRIC_GP_N    := %s,
MK_TC_DIRECTORY          := %s,
TRUE_CONES_DIRECTORY     := "%s",
FILL_TRUE_CONES          := %s,
TC_MIN_PRIME             := %s,
TC_MAX_PRIME             := %s,
TC_MIN_SYMMETRIC_GP_N    := %s,
TC_MAX_SYMMETRIC_GP_N    := %s,
MK_BC_DIRECTORY          := %s,
BLOCKED_CONES_DIRECTORY  := "%s",
FILL_BC_DIRECTORY        := %s,
OVERWRITE_OLD_BCS        := %s,
MK_IC_DIRECTORY          := %s,
INT_CONES_DIRECTORY      := "%s",
FILL_IC_DIRECTORY        := %s,
OVERWRITE_OLD_ICS        := %s  
);
""" % (
GAP_FILENAME            ,
SAGE_FILENAME           ,
GAP_PROGRAM_FILENAME    ,
SAGE_PROGRAM_FILENAME   ,
MK_RC_DIRECTORY         ,
RAW_CONES_DIRECTORY     ,
FILL_RAW_CONES          ,
RC_MIN_PRIME            ,
RC_MAX_PRIME            ,
RC_MIN_SYMMETRIC_GP_N   ,
RC_MAX_SYMMETRIC_GP_N   ,
MK_TC_DIRECTORY         ,
TRUE_CONES_DIRECTORY    ,
FILL_TRUE_CONES         ,
TC_MIN_PRIME            ,
TC_MAX_PRIME            ,
TC_MIN_SYMMETRIC_GP_N   ,
TC_MAX_SYMMETRIC_GP_N   ,
MK_BC_DIRECTORY         ,
BLOCKED_CONES_DIRECTORY ,
FILL_BC_DIRECTORY       ,
OVERWRITE_OLD_BCS       ,
MK_IC_DIRECTORY         ,
INT_CONES_DIRECTORY     ,
FILL_IC_DIRECTORY       ,
OVERWRITE_OLD_ICS       ,
    ))
    f.flush()
    f.close()
else:
    pass # TODO load the relevant settings.

subprocess.run(["sh", GAP_FILENAME, "-l", ";%s"%INTERSECTION_CONES_TOP_DIRECTORY, GAP_PROGRAM_FILENAME])
subprocess.run([SAGE_FILENAME, SAGE_PROGRAM_FILENAME])


# Sage part:

# load("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/prog/sageprog.sage")

