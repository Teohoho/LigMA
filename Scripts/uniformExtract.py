import mdtraj as md
import sys, os

## USAGE:
##	uniformExtract.py NAME.prmtop NAME.dcd STRIDE OUTDIR 
###DEBUG###
#print (sys.argv[1])
#print (sys.argv[2])
#print (sys.argv[3])
#print (sys.argv[4])
###DEBUG###

MDTrajObj = md.load(sys.argv[2], top=sys.argv[1], stride=int(sys.argv[3]))
MDTrajObj.superpose(MDTrajObj, frame=0)
os.mkdir (sys.argv[4])
rootName = sys.argv[1].split("/")[-1].split(".")[:-1]
print (rootName)
for frameIx in range(MDTrajObj.n_frames):
	os.mkdir("{}/{}.frame{}".format(sys.argv[4],
					rootName[0],
					frameIx*int(sys.argv[3])))
	currentFrame = MDTrajObj[frameIx]
	fileOut = "{0}/{1}.frame{2}/{1}.frame{2}.pdb".format(
					sys.argv[4],rootName[0],
					frameIx*int(sys.argv[3]))
	print (fileOut)
	currentFrame.save_pdb(fileOut, force_overwrite=False)

print ("PDBs saved!")
