import os
import sys
import re
import shutil
import fileinput
import subprocess
import math
import time
from telnetlib import theNULL
from os.path import exists
from os.path import isfile

# ----------------------------------------------------------------------------------------- #

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# ----------------------------------------------------------------------------------------- #

def removeFile(pathToRemove):
    if (os.path.exists(pathToRemove) or pathToRemove[len(pathToRemove) -1] == "*"):
        subprocess.run(["rm","-r",pathToRemove])

# ----------------------------------------------------------------------------------------- #

def copyFile(pathToCopy,pathToPaste):
    if (os.path.exists(pathToCopy) or pathToCopy[len(pathToCopy) -1] == "*"):
        subprocess.run(["scp","-r",pathToCopy,pathToPaste])
    else:
        sys.exit("Error: "+pathToCopy+" does not exist!")

# ----------------------------------------------------------------------------------------- #

def moveFile(pathToOrigin,pathToDestination):
    subprocess.run(["mv",pathToOrigin,pathToDestination])

# ----------------------------------------------------------------------------------------- #

def getArguments():
    
    global maxAllowedImbalance
    global weightField
    global logFile
    global solver

    for arg in sys.argv:
        if "-i" in arg:
            splitArg=arg.split("=")
            maxAllowedImbalance=float(splitArg[1])
        if "-w" in arg:
            splitArg=arg.split("=")
            weightField=str(splitArg[1])
        if "-l" in arg:
            splitArg=arg.split("=")
            logFile=str(splitArg[1])
        if "-s" in arg:
            splitArg=arg.split("=")
            solver=str(splitArg[1])
            
    if (weightField==""):
        print("Error load balancing weight field must be given as argument...")
        exit()
                        
# ----------------------------------------------------------------------------------------- #

def restoreDictionaries():

    pathToDictOrig="./system/controlDict.orig"
    pathToDict="./system/controlDict"
    if (os.path.exists(pathToDictOrig)):
        removeFile(pathToDict)
        copyFile(pathToDictOrig,pathToDict)
    else:
        copyFile(pathToDict,pathToDictOrig)

    pathToDictOrig="./system/decomposeParDict.orig"
    pathToDict="./system/decomposeParDict"
    if (os.path.exists(pathToDictOrig)):
        removeFile(pathToDict)
        copyFile(pathToDictOrig,pathToDict)
    else:
        copyFile(pathToDict,pathToDictOrig)

# ----------------------------------------------------------------------------------------- #

def getRunInfo():

    global deltaT
    global writeTime
    global endTime
    global nProcs
    global weightField

    #get deltaT, writeTime, endTime
    pathToFile="./system/controlDict"
    with open(pathToFile, "r") as inFile:
        lines = inFile.readlines()

    for line in lines:
        if "deltaT" in line:
            lineSplit = re.split(' |;|\n',line)
            for string in lineSplit:
                if (is_number(string)):
                    deltaT=float(string)
                    break          
        if "writeInterval" in line:
            lineSplit = re.split(' |;|\n',line)
            for string in lineSplit:
                if (is_number(string)):
                    writeTime=float(string)
                    break                      
        if "endTime" in line:
            lineSplit = re.split(' |;|\n',line)
            for string in lineSplit:
                if (is_number(string)):
                    endTime=float(string)
                    break
                
    #get number of cfd and usp processors
    pathToFile="./system/decomposeParDict"
    with open(pathToFile, "r") as inFile:
            lines = inFile.readlines()

    for line in lines:
        if "numberOfSubdomains" in line:
            lineSplit = re.split(' |;|\n',line)
            for string in lineSplit:
                if (is_number(string)):
                    nProcs=int(string)
                    break        

    #avoid writing weightField more than one times
    addLine = True  
    
    #add weightField entry to decomposeParDict
    with open(pathToFile, "w") as outFile:
        for line in lines:
            outFile.write(line)
            if "method" in line and addLine:
                addLine = False
                outFile.write("weightField     "+weightField+";\n")

# ----------------------------------------------------------------------------------------- #

def checkImbalance():

    global isEndTime
    global writeTime
    global endTime

    #check if end time is reached
    isEndTime = False
    for dir in os.listdir("./processor0/"):
        if (is_number(dir) and float(dir) >= endTime-writeTime/2.e0 and float(dir) <= endTime+writeTime/2.e0):
            isEndTime = True

    with open(logFile, "r") as inFile:
        lines = inFile.readlines()

    #find maximum imbalance
    maxImbalance = -1 
    for lineI in range(max(0,len(lines)-20),len(lines)):
        if "Maximum imbalance" in lines[lineI]:
            lineSplit = re.split(' |;|\n',lines[lineI])
            for string in lineSplit:
                if (is_number(string)):
                    maxImbalance = float(string)

    if (maxImbalance >= maxAllowedImbalance):

        time.sleep(5)

        print("Found imbalance "+"{:8.2f}".format(maxImbalance)+"% larger than allowed imbalance "+"{:8.2f}".format(maxAllowedImbalance)+"%.")

        #find last time
        for lineI in range(max(0,len(lines)-20),len(lines)):
            if lines[lineI].startswith("Time"):
                lineSplit = re.split(' |;|\n',lines[lineI])
                for string in lineSplit:
                    if (is_number(string)):
                        currentTime = float(string)

        nextWriteTime = (int(currentTime/writeTime)+1)*writeTime

        print("Rebalancing simulation at time "+str(nextWriteTime)+".")
        print("")

        #stop simulation by changing controlDict
        pathToFile="./system/controlDict"
        with open(pathToFile, "r") as inFile:
            lines = inFile.readlines()

        with open(pathToFile, "w") as outFile:
            for line in lines:
                #if "stopAt" in line:
                #    outFile.write("stopAt writeNow;\n")
                if line.startswith("endTime"):
                    outFile.write("endTime "+str(nextWriteTime)+";\n")
                else:
                    outFile.write(line)

        #check if simulation has stopped by reading logFile
        simulationEnded = False
        while (not simulationEnded):
            with open(logFile, "r") as inFile:
                    lines = inFile.readlines()
            for lineI in range(max(0,len(lines)-5),len(lines)):
                if "Finalising parallel run" in lines[lineI]:
                    simulationEnded = True
                    break
            time.sleep(1)

        #reconstruct latest time directory
        reconstructLatestTime()

        #save any additional result directories that exist in processor0/
        dirToSave = list()
        for dir in os.listdir("./processor0/"):
            if (not is_number(dir) and dir != "constant" and dir != "system"):
                dirToSave.append(dir)
        for dir in dirToSave:
            copyFile("./processor0/"+dir,"./")
               
        #decompose latest time directory
        decomposeLatestTime()

        #restore any additional result directories that existed in processor0/
        for dir in dirToSave:
            copyFile("./"+dir,"./processor0/")
            
        #restart simulation from latest time
        pathToFile="./system/controlDict"
        with open(pathToFile, "r") as inFile:
            lines = inFile.readlines()

        with open(pathToFile, "w") as outFile:
            for line in lines:
                #if "stopAt" in line:
                #    outFile.write("stopAt endTime;\n")
                if line.startswith("endTime"):
                    outFile.write("endTime "+str(endTime)+";\n")             
                elif  "startFrom" in line:
                    outFile.write("startFrom latestTime;\n")
                else:
                    outFile.write(line)

        runSolver()

        time.sleep(1)

# ----------------------------------------------------------------------------------------- #

def reconstructLatestTime():
    runCommand = "reconstructPar -latestTime >> "+logFile
    subprocess.call([runCommand],shell=True)

# ----------------------------------------------------------------------------------------- #

def decomposeLatestTime():
    runCommand = "decomposePar -latestTime -force >> "+logFile
    subprocess.call([runCommand],shell=True)

# ----------------------------------------------------------------------------------------- #

def runSolver():
    runCommand = "mpirun -n "+str(nProcs)+" "+solver+" -parallel >> "+logFile+" &"
    subprocess.call([runCommand],shell=True)

# ----------------------------------------------------------------------------------------- #

#user defined variables
maxAllowedImbalance = 10.0    #maximum allowed imbalance between processors in percententage
weightField = ""              #load balancing weight field name
logFile = "log.uspFoam"      #output logfile name
solver = "uspFoam"           #solver 

#initialize global variables
isEndTime=False
deltaT=0
writeTime=0
endTime=0
nProcs=0

#get intitial time
tStart = time.time()

#read user defined variables from arguments
getArguments()

#save controlDict as controlDict.orig and decomposeParDict as decomposeParDict.orig
restoreDictionaries()

#get parallel run info
getRunInfo()

print("------------------------------------------------------")
print("Running "+solver+" on "+str(nProcs)+" cores.")
print("Maximum allowed imbalance is "+str(maxAllowedImbalance)+"%.")
print("Load balancing weight field is "+weightField+".")
print("------------------------------------------------------")
print("")

#start parallel run
runSolver()

#check if logFile is beign written
while not os.path.exists(logFile):
    time.sleep(1)

#load balance simulation
while (not isEndTime):

    checkImbalance()

#restore original controlDict
restoreDictionaries()

tEnd = time.time()
print("cpu time: "+"{:10.6f}".format((tEnd-tStart)/3600e0)+" hours.")