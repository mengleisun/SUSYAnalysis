import os
from os import system, environ
path ="/eos/uscms/store/user/tmishra/elefakepho/DYResult18/"
outDir="DY18new"

system("cat {}*Bw-expo-pt*.txt >> {}/EleFakeRate-Data-Bw-expo-pt-60-120.txt".format(path,outDir))
system("cat {}*Bw-ker-pt*.txt >> {}/EleFakeRate-Data-Bw-ker-pt-60-120.txt".format(path,outDir))
system("cat {}*DY-ker-pt*.txt >> {}/EleFakeRate-Data-DY-ker-pt-60-120.txt".format(path,outDir))

system("cat {}*Bw-expo-eta*.txt >> {}/EleFakeRate-Data-Bw-expo-eta-60-120.txt".format(path,outDir))
system("cat {}*Bw-ker-eta*.txt >> {}/EleFakeRate-Data-Bw-ker-eta-60-120.txt".format(path,outDir))
system("cat {}*DY-ker-eta*.txt >> {}/EleFakeRate-Data-DY-ker-eta-60-120.txt".format(path,outDir))

system("cat {}*Bw-expo-vtx*.txt >> {}/EleFakeRate-Data-Bw-expo-vtx-60-120.txt".format(path,outDir))
system("cat {}*Bw-ker-vtx*.txt >> {}/EleFakeRate-Data-Bw-ker-vtx-60-120.txt".format(path,outDir))
system("cat {}*DY-ker-vtx*.txt >> {}/EleFakeRate-Data-DY-ker-vtx-60-120.txt".format(path,outDir))
