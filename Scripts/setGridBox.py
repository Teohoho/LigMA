import numpy as np
import ls_parsepdb
import sys, glob, math
import mdtraj as md
# 
#       A
#      /|\
#     / | \
#    /  |  \
#   /   |   \
#  /____|____\
# -------------
# B     H     C
#
#	USAGE:
#	python setGridBox receptor.pdb ligandDir
def getDist(P1,P2):
#	print(P1,P2)
	dist = math.sqrt((P1[0]-P2[0])**2 + (P1[1]-P2[1])**2 + (P1[2]-P2[2])**2)
	return dist


def getGridFromTriangle(A, B, C, ligandSize):
	# Check
	if A.shape != (3,):
		print("Point A should be of shape (3).")
	if B.shape != (3,):
		print("Point B should be of shape (3).")
	if C.shape != (3,):
		print("Point C should be of shape (3).")

	# Get vectors
	BA = A - B
	BC = C - B

	# Get lengths (didn't trust Numpy)
	ba = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
	bc = np.sqrt(BC[0]**2 + BC[1]**2 + BC[2]**2)

	# Get theta
	cosTheta = np.dot(BA, BC) / (ba * bc)
	#theta = np.arccos(cosTheta)
	bh = cosTheta * ba
	
	# Get H and HA
	BH = (bh / bc) * BC
	H = B + BH
	HA = A - H
	ha = np.sqrt(HA[0]**2 + HA[1]**2 + HA[2]**2)

	# Grids centers
	xoffVec = (0.25*BC) - (np.array([ligandSize, 0.0, 0.0]) * 0.25)
	X1 = B + xoffVec
	X2 = C - xoffVec

	yoffVec = (HA / 2.0) + np.array([0.0, ligandSize/2.0, 0.0])
	Y1 = X1 + yoffVec
	Y2 = X2 + yoffVec

	gridCenters = np.array([Y1, Y2])

	# Get grids frame
	XVec = BH / bh
	YVec = HA / ha
	ZVec = np.cross(XVec, YVec)
	frame = np.array([XVec, YVec, ZVec])

	# Grids size
	size = np.array([(bc / 2.0) + ligandSize, ha + ligandSize, ha + ligandSize])
	size /= 0.375
 
	return (gridCenters, frame, size)

#

# Main
#A = np.array([35.471 , 34.934 , 26.984])
#B = np.array([4.101 , 7.737 , 38.821])
#C = np.array([42.418 , -2.093 , -4.947])

#A = np.array([60.470, -32.199, 42.585]) # index 3864
#B = np.array([-53.000, -2.109, 22.064]) # index 472 
#C = np.array([11.880, 0.927, 6.025]) # index 2637

if (len(sys.argv) != 4):
	print ("Usage: setGridBox.py PDBIn ligandDir outDir")
	print (sys.argv)
	sys.exit()

pdbParser = ls_parsepdb.ParsePdb()
pdbParser.Read(sys.argv[1])

#A_atom = pdbParser.getAtomByIndex(4307)
# This is Index, not Serial (in VMD)
# These atoms were chosen, I assume, so as to form a plane
# that encompasses the whole NS2 dimer.

A_atom = pdbParser.getAtomByIndex(1420)
#print(A_atom)
A = np.array([A_atom[0]['x'], A_atom[0]['y'], A_atom[0]['z']])

B_atom = pdbParser.getAtomByIndex(471)
#print(B_atom)
B = np.array([B_atom[0]['x'], B_atom[0]['y'], B_atom[0]['z']])

C_atom = pdbParser.getAtomByIndex(2636)
#print(C_atom)
C = np.array([C_atom[0]['x'], C_atom[0]['y'], C_atom[0]['z']])

# I assume this is the size of the largest ligand
maxDist = 0
fileList = glob.glob(sys.argv[2])
#print (fileList)
for fileix in range(len(fileList)):
	#print (fileList[fileix])
	MDTrajObj = md.load(fileList[fileix])
	for ix in range(MDTrajObj.n_atoms):
		for jx in range(ix,MDTrajObj.n_atoms):
			if (getDist(MDTrajObj.xyz[0][ix],MDTrajObj.xyz[0][jx]) > maxDist):
				maxDistID = fileix
				maxDist = getDist(MDTrajObj.xyz[0][ix],MDTrajObj.xyz[0][jx])
print ("Max ligand distance: {}, ligand: {}".format(maxDist,fileList[maxDistID]))
ligandSize = maxDist

gridCenters, frame, sizes = getGridFromTriangle(A, B, C, ligandSize)

#rotationMatrix = frame
#pdbParser.Rotate(rotationMatrix)
#moveVector = -1.0 * gridCenters[0]
#pdbParser.Move(moveVector)
#pdbParser.PrintPdb()

rootName = sys.argv[1].split("/")[-1].split(".pdb")[:-1]
for i in range(2):
	with open("{}/{}.Site{}.GPFRefference.txt".format(sys.argv[3], rootName[0], i+1), "w") as f:
		f.write("npts {} {} {}\n".format(int(sizes[0]),
													int(sizes[1]),
													int(sizes[2])))
		f.write("spacing 0.375\n")
		f.write("gridcenter {} {} {}\n".format(gridCenters[i][0],
															gridCenters[i][1],		
															gridCenters[i][2]))

#print("Grid Centers")
#print(gridCenters)
#print("Frame")
#print(frame)
#print("Sizes")
#print(sizes)
