import mdtraj as md
import numpy as np
import glob, sys
from scipy.spatial.transform import Rotation as R

if (len(sys.argv)) < 3:
    print ("USAGE: \n{} PDBIn DirContainingLigs".format(sys.argv[0]))
    sys.exit()

ligandList = glob.glob("{}/*pdb".format(sys.argv[2]))

## Get longest ligand
maxDist = 0
for fileix in range(len(ligandList)):
    MDTrajObj = md.load(ligandList[fileix])
    for ix in range(MDTrajObj.n_atoms):
        for jx in range(ix,MDTrajObj.n_atoms):
            if (np.linalg.norm(MDTrajObj.xyz[0][ix] - MDTrajObj.xyz[0][jx]) > maxDist):
                maxDistID = fileix
                maxDist = np.linalg.norm(MDTrajObj.xyz[0][ix] - MDTrajObj.xyz[0][jx])
print ("Max ligand distance: {}, ligand: {}".format(maxDist,ligandList[maxDistID]))

MDTrajObj = md.load(sys.argv[1])

## Compute inertia tensor (inTens)
inTens = md.compute_inertia_tensor(MDTrajObj)

## Compute eigenvalue/eigenvector of inTens
eVal,eVec = np.linalg.eig(inTens)

## Translate the protein to (0, 0, 0)
newPositions = MDTrajObj.xyz[0] - md.compute_center_of_mass(MDTrajObj) 

## Apply the rotation object on the positions of the protein
## Multiplied by 10 since it's Angstroms, not nm
rotObj = R.from_matrix(-eVec.T.reshape(1,3,3))
newPositions = rotObj.apply(newPositions)
newPositions = newPositions * 10

## Save the new coordinates into a new PDB
outputFile = "{}.Aligned.pdb".format(sys.argv[1].replace(".pdb", ""))
newPDB = md.formats.PDBTrajectoryFile(outputFile, "w")
newPDB.write(newPositions, MDTrajObj.topology)

#####################################
## Blind Docking
#####################################

if (len(sys.argv) == 3):

    ## Compute minimum/maximum for all coords
    Xmin = np.min(newPositions[:,0])
    Xmax = np.max(newPositions[:,0])
    Ymin = np.min(newPositions[:,1])
    Ymax = np.max(newPositions[:,1])
    Zmin = np.min(newPositions[:,2])
    Zmax = np.max(newPositions[:,2])
    
    print ("Max coords: {} {} {} \nMin coords: {} {} {}".format(Xmax, Ymax, Zmax, Xmin, Ymin, Zmin))
    spacing = 0.375
    offSet = maxDist * 10 * 0.5
    
    ## If you know what you're doing, you can modify these individual offsets
    XPlusOffset = 0
    XMinusOffset = 0
    YPlusOffset = 0
    YMinusOffset = 0
    ZPlusOffset = 0
    ZMinusOffset = 0

    Xmax = Xmax + XPlusOffset
    Xmin = Xmin + XMinusOffset
    Ymax = Ymax + YPlusOffset
    Ymin = Ymin + YMinusOffset
    Zmax = Zmax + ZPlusOffset
    Zmin = Zmin + ZMinusOffset
    
    ## Box Size
    XSize = (Xmax + offSet) - (Xmin - offSet)
    YSize = (Ymax + offSet) - (Ymin - offSet)
    ZSize = (Zmax + offSet) - (Zmin - offSet)
    XPts = XSize / spacing
    YPts = YSize / spacing 
    ZPts = ZSize / spacing  
    center = [(Xmax+Xmin)/2, (Ymax+Ymin)/2, (Zmax+Zmin)/2]

    print ("Grid Size:\n{} {} {}".format(XSize, YSize, ZSize))
    print ("Points: {} {} {}".format(XPts, YPts, ZPts))
    print ("Grid Center\n {} {} {}".format(center[0],center[1],center[2]))
    
    ## Write GPFRefference File
    with open("{}.GPFRefference.txt".format(sys.argv[1].replace(".pdb", "")), "w") as f:
        f.write("npts {} {} {}\n".format(int(XPts), int(YPts), int(ZPts)))
        f.write("spacing 0.375\n")
        f.write("gridcenter {} {} {}\n".format(center[0],center[1],center[2]))

#####################################
## Targeted Docking 
#####################################

else:

    ## Get indices of Binding Site selection
    BSIndices = MDTrajObj.topology.select(sys.argv[3])

    ## Compute minimum/maximum of chosen binding site
    BSPositions = np.zeros((len(BSIndices),3))
    for atomix in range(len(BSIndices)):
        BSPositions[atomix] = newPositions[atomix]
    
    ## Compute minimum/maximum for all coords
    Xmin = np.min(BSPositions[:,0])
    Xmax = np.max(BSPositions[:,0])
    Ymin = np.min(BSPositions[:,1])
    Ymax = np.max(BSPositions[:,1])
    Zmin = np.min(BSPositions[:,2])
    Zmax = np.max(BSPositions[:,2])

    print ("Max coords: {} {} {} \nMin coords: {} {} {}".format(Xmax, Ymax, Zmax, Xmin, Ymin, Zmin))
    spacing = 0.375
    offSet = maxDist * 10
    
    ## If you know what you're doing, you can modify these individual offsets
    XPlusOffset = 0
    XMinusOffset = 0
    YPlusOffset = 0
    YMinusOffset = 0
    ZPlusOffset = 0
    ZMinusOffset = 0
    
    Xmax = Xmax + XPlusOffset
    Xmin = Xmin + XMinusOffset
    Ymax = Ymax + YPlusOffset
    Ymin = Ymin + YMinusOffset
    Zmax = Zmax + ZPlusOffset
    Zmin = Zmin + ZMinusOffset
    
    ## Box Size
    XSize = (Xmax - Xmin) + offSet
    YSize = (Ymax - Ymin) + offSet
    ZSize = (Zmax - Zmin) + offSet
    print ("Size: {} {} {}".format(XSize, YSize, ZSize))
    XPts = XSize / spacing
    YPts = YSize / spacing 
    ZPts = ZSize / spacing
    
    ## Box Size
    XSize = (Xmax + offSet) - (Xmin - offSet)
    YSize = (Ymax + offSet) - (Ymin - offSet)
    ZSize = (Zmax + offSet) - (Zmin - offSet)
    XPts = XSize / spacing
    YPts = YSize / spacing 
    ZPts = ZSize / spacing  
    center = [np.average(BSPositions[:,0]),np.average(BSPositions[:,1]),np.average(BSPositions[:,2])]

    print ("Grid Size:\n{} {} {}".format(XSize, YSize, ZSize))
    print ("Points: {} {} {}".format(XPts, YPts, ZPts))
    print ("Grid Center\n {} {} {}".format(center[0], center[1], center[2]))
    
    ## Write GPFRefference File
    with open("{}.GPFRefference.txt".format(sys.argv[1].replace(".pdb", "")), "w") as f:
        f.write("npts {} {} {}\n".format(int(XPts), int(YPts), int(ZPts)))
        f.write("spacing 0.375\n")
        f.write("gridcenter {} {} {}\n".format(center[0], center[1], center[2]))

