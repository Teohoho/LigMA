import sys

with open(sys.argv[1], "r") as f:
	lines = f.readlines()
	for i in range(len(lines)-1):
		if (" H2 " not in lines[i+1]):
			print(lines[i],end="")
		if ("OXT" in lines[i]):
			print("TER")
