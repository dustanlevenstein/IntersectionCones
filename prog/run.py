import subprocess, os

GAP_FILENAME = os.path.expanduser("~/Gap/gap4r8/bin/gap.sh")
GAP_PROGRAM_FILENAME = os.path.expanduser("~/Dropbox/MyFuture/ProgrammingProjects/IntersectionCones/IntersectionCones/prog/prog.g")
SAGE_FILENAME = os.path.expanduser("~/Sage/SageMath/sage")
SAGE_PROGRAM_FILENAME = os.path.expanduser("~/Dropbox/MyFuture/ProgrammingProjects/IntersectionCones/IntersectionCones/prog/sageprog.sage")


subprocess.run(["sh", GAP_FILENAME, GAP_PROGRAM_FILENAME])
subprocess.run([SAGE_FILENAME])


# Sage part:

# load("../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/prog/sageprog.sage")

