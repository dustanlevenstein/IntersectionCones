# IntersectionCones

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

MK_RC_DIRECTORY          = True  # Change to False when directory exists
RAW_CONES_DIRECTORY      = "%s/raw_cones/"
FILL_RAW_CONES           = True
RC_MIN_PRIME             = 2
RC_MAX_PRIME             = 7
RC_MIN_SYMMETRIC_GP_N    = 1
RC_MAX_SYMMETRIC_GP_N    = 6

MK_TC_DIRECTORY          = True  # Change to False when directory exists
TRUE_CONES_DIRECTORY     = "%s/true_cones/"
FILL_TRUE_CONES          = True
TC_MIN_PRIME             = 2
TC_MAX_PRIME             = 7
TC_MIN_SYMMETRIC_GP_N    = 1
TC_MAX_SYMMETRIC_GP_N    = 20

MK_BC_DIRECTORY          = True  # Change to False when directory exists
BLOCKED_CONES_DIRECTORY  = "%s/blocked_cones/"
FILL_BC_DIRECTORY        = True  # This will fill in blocked cones for every
                                 # file which exists in the raw cones
                                 # directory.
OVERWRITE_OLD_BCS        = False # Overwrite already existing blocked cones.

MK_IC_DIRECTORY          = True  # Change to False when directory exists
INT_CONES_DIRECTORY      = "%s/intersected_cones/"
FILL_IC_DIRECTORY        = True  # This will fill in the intersection of
                                 # every blocked cone which exists in the
                                 # blocked cones directory.
OVERWRITE_OLD_ICS        = False # Overwrite already existing intersections
                                 # of cones.


