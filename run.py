import os
if os.path.isfile("paths.txt"):
    f = open("paths.txt", "r")
    directories = eval(f.read())
    f.close()
    INTERSECTION_CONES_TOP_DIRECTORY = directories["INTERSECTION_CONES_TOP_DIRECTORY"]
    GAP_FILENAME  = directories["GAP_FILENAME"]
    SAGE_FILENAME = directories["SAGE_FILENAME"]
else:
    INTERSECTION_CONES_TOP_DIRECTORY = input("Enter the path for the IntersectionCones directory. ")
    


                              # Where you see "%s", read the top directory in
                              # its place.
GAP_PROGRAM_FILENAME         = "%s/prog/prog.g"
SAGE_PROGRAM_FILENAME        = "%s/prog/sageprog.sage"

RAW_CONES_DIRECTORY          = "%s/raw_cones/"

RC_MIN_PRIME                 = input("Indicate the smallest prime number for which"
                                " you wish to generate Hecke algebra cones."
                                " (2 would be a good choice.) ")
RC_MAX_PRIME                 = input("Indicate the largest prime number for which"
                                " you wish to generate Hecke algebra cones."
                                " (e.g., 5.) ")
RC_MIN_SYMMETRIC_GP_N        = input("Indicate the smallest symmetric group for"
                                 " which you wish to generate cones.")
RC_MAX_SYMMETRIC_GP_N        = input("Indicate the largest symmetric group for"
                                 " which you wish to generate cones.")

TRUE_CONES_DIRECTORY         = "%s/true_cones/"
TC_MIN_PRIME                 = input("Indicate the smallest prime number for which"
                             =  " you wish to generate Hecke algebra cones."
                             =  " (2 would be a good choice.) ")
TC_MAX_PRIME                 = input("Indicate the largest prime number for which"
                             =  " you wish to generate Hecke algebra cones."
                             =  " (e.g., 5.) ")
TC_MIN_SYMMETRIC_GP_N        = input("Indicate the smallest symmetric group for"
                             =   " which you wish to generate cones.")
TC_MAX_SYMMETRIC_GP_N        = input("Indicate the largest symmetric group for"
                             =   " which you wish to generate cones.")


BLOCKED_CONES_DIRECTORY      = "%s/blocked_raw_cones/"
BLOCKED_TRUE_CONES_DIRECTORY = "%s/blocked_true_cones/"

DUAL_CONES_DIRECTORY         = "%s/dual_cones/"
INT_CONES_DIRECTORY          = "%s/intersected_cones/"

F_I_MACHINE_DIRECTORY        = "%s/f_i_machine/"
F_I_HUMAN_DIRECTORY          = "%S/f_i_human/"
IND_HUMAN_DIRECTORY          = "%s/ind_human/"
IND_MACHINE_DIRECTORY        = "%s/ind_machine/"

POSS_MACHINE_DIRECTORY       = "%s/poss_adj_machine/"
POSS_HUMAN_DIRECTORY         = "%s/poss_adj_human/"


import subprocess

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
GAP_FILENAME                 := "%s",
SAGE_FILENAME                := "%s",
GAP_PROGRAM_FILENAME         := "%s",
SAGE_PROGRAM_FILENAME        := "%s",
RAW_CONES_DIRECTORY          := "%s",
RC_MIN_PRIME                 := %s,
RC_MAX_PRIME                 := %s,
RC_MIN_SYMMETRIC_GP_N        := %s,
RC_MAX_SYMMETRIC_GP_N        := %s,
TRUE_CONES_DIRECTORY         := "%s",
TC_MIN_PRIME                 := %s,
TC_MAX_PRIME                 := %s,
TC_MIN_SYMMETRIC_GP_N        := %s,
TC_MAX_SYMMETRIC_GP_N        := %s,
BLOCKED_CONES_DIRECTORY      := "%s",
BLOCKED_TRUE_CONES_DIRECTORY := "%s",
DUAL_CONES_DIRECTORY         := "%s",
INT_CONES_DIRECTORY          := "%s",
F_I_MACHINE_DIRECTORY        := "%s",
F_I_HUMAN_DIRECTORY          := "%s",
IND_HUMAN_DIRECTORY          := "%s",
IND_MACHINE_DIRECTORY        := "%s",
POSS_MACHINE_DIRECTORY       := "%s",
POSS_HUMAN_DIRECTORY         := "%s",
);
""" % (
GAP_FILENAME                ,
SAGE_FILENAME               ,
GAP_PROGRAM_FILENAME        ,
SAGE_PROGRAM_FILENAME       ,
RAW_CONES_DIRECTORY         ,
RC_MIN_PRIME                ,
RC_MAX_PRIME                ,
RC_MIN_SYMMETRIC_GP_N       ,
RC_MAX_SYMMETRIC_GP_N       ,
TRUE_CONES_DIRECTORY        ,
TC_MIN_PRIME                ,
TC_MAX_PRIME                ,
TC_MIN_SYMMETRIC_GP_N       ,
TC_MAX_SYMMETRIC_GP_N       ,
BLOCKED_CONES_DIRECTORY     ,
BLOCKED_TRUE_CONES_DIRECTORY,
DUAL_CONES_DIRECTORY        ,
INT_CONES_DIRECTORY         ,
F_I_MACHINE_DIRECTORY       ,
F_I_HUMAN_DIRECTORY         ,
IND_HUMAN_DIRECTORY         ,
IND_MACHINE_DIRECTORY       ,
POSS_MACHINE_DIRECTORY      ,
POSS_HUMAN_DIRECTORY        
    ))
f.flush()
f.close()

subprocess.run(["sh", GAP_FILENAME, "-l", ";%s"%INTERSECTION_CONES_TOP_DIRECTORY, GAP_PROGRAM_FILENAME])
subprocess.run([SAGE_FILENAME, SAGE_PROGRAM_FILENAME])


# Sage part:

# load("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/prog/sageprog.sage")

