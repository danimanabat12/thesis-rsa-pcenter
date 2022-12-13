#!/opt/homebrew/bin/python3
#  Dili siya simple usr/bin since daghan kog python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys, os

def main():
  if (len(sys.argv) < 3 or len(sys.argv) > 3): 
    print("Wrong argument. Try again")
    exit()
  else:
    currPC = sys.argv[1]
    currDataset = sys.argv[2]

  path = os.getcwd()+f"/_Output/{currPC}"
  os.makedirs(path)

  generateConvergenceMap(currPC, currDataset)
  generateDatasetMap(currPC, currDataset)
  generatePlottedFacilities(currPC, currDataset)
  generateGenerationFitnessPlot(currPC, currDataset)
  generateRuntime(currPC, currDataset)
  generateBestFitTime(currPC, currDataset)
  generatePhasePerformance(currPC, currDataset)

def generateConvergenceMap(currPC, currDataset):
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  num_iters = len(runFitness)
  x = np.linspace(0, num_iters,num_iters)

  for i in range(0,len(runFitness.iloc[0])):
      plt.plot(x,runFitness[i],'g-')

  plt.title(f'Convergence Map of 6-center Problem on {currDataset} City Household Data\nUsing RSA with T={num_iters-1}')
  plt.xlabel("Iteration Number")
  plt.ylabel("Fitness Value")
  plt.savefig(f'_Output/{currPC}/{currPC}_ConvergenceMap.png')
  plt.clf()

def generateDatasetMap(currPC, currDataset):
  DX = pd.read_csv(f"{currPC}/{currDataset}_X.csv",header=None, skiprows=1)
  DY = pd.read_csv(f"{currPC}/{currDataset}_Y.csv",header=None, skiprows=1)
  plt.plot(DX, DY, 'x', color='black',markersize=.2);
  plt.title(f'Demand points based on {currPC}\n City, Philippines Household Data.')
  plt.xlabel('x-coordinate')
  plt.ylabel('y-coordinate')
  plt.savefig(f'_Output/{currPC}/{currPC}_DatasetMap.png')
  plt.clf()

def generatePlottedFacilities(currPC, currDataset):
  DX = pd.read_csv(f"{currPC}/{currDataset}_X.csv",header=None, skiprows=1)
  DY = pd.read_csv(f"{currPC}/{currDataset}_Y.csv",header=None, skiprows=1)
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  initXValues = pd.read_csv(f"{currPC}/{currPC}_RSA_InitBestX.csv",header=None)
  initYValues = pd.read_csv(f"{currPC}/{currPC}_RSA_InitBestY.csv",header=None)
  bestXValues = pd.read_csv(f"{currPC}/{currPC}_RSA_FinalBestX.csv",header=None)
  bestYValues = pd.read_csv(f"{currPC}/{currPC}_RSA_FinalBestY.csv",header=None)

  A = runFitness.iloc[len(runFitness)-1]
  index_min = np.argmin(A)

  spx=initXValues.iloc[index_min]
  spy=initYValues.iloc[index_min]
  bestX=bestXValues.iloc[index_min]
  bestY=bestYValues.iloc[index_min]

  plt.gca().set_aspect('equal', adjustable='box')

  ax = plt.gca()

  if (currDataset == 'Digos'):
    ax.set_ylim([min(bestYValues.iloc[index_min])-7000,max(bestYValues.iloc[index_min])+9000])
    ax.set_xlim([min(bestXValues.iloc[index_min])-7000,max(bestXValues.iloc[index_min])+12000])
  elif (currDataset == 'Davao'):
    ax.set_ylim([min(bestYValues.iloc[index_min])-13000,max(bestYValues.iloc[index_min])+10000])
    ax.set_xlim([min(bestXValues.iloc[index_min])-11000,max(bestXValues.iloc[index_min])+13000])
  else:
    ax.set_ylim([min(bestYValues.iloc[index_min])-10000,max(bestYValues.iloc[index_min])+12000])
    ax.set_xlim([min(bestXValues.iloc[index_min])-13000,max(bestXValues.iloc[index_min])+13000])

  # Plotting of demand poiints
  plt.scatter(DX,DY,color='black', marker='o',s=1)

  # Plotting of best facilities
  plt.scatter(bestX,bestY,color='red', marker='^')

  # Plotting of initial facilities
  plt.scatter(spx,spy,color='green')

  # # plot the radius/fitness of the best
  circle1=plt.Circle(( bestX[0],bestY[0]), runFitness[index_min][len(runFitness)-1], fill=False,color='r')
  circle2=plt.Circle((bestX[1],bestY[1]), runFitness[index_min][len(runFitness)-1], fill=False,color='r')
  circle3=plt.Circle((bestX[2],bestY[2]), runFitness[index_min][len(runFitness)-1], fill=False,color='r')
  circle4=plt.Circle(( bestX[3],bestY[3]), runFitness[index_min][len(runFitness)-1], fill=False,color='r')
  circle5=plt.Circle((bestX[4],bestY[4]), runFitness[index_min][len(runFitness)-1], fill=False,color='r')
  circle6=plt.Circle((bestX[5],bestY[5]), runFitness[index_min][len(runFitness)-1], fill=False,color='r')
  ax.add_artist(circle1)
  ax.add_artist(circle2)
  ax.add_artist(circle3)
  ax.add_artist(circle4)
  ax.add_artist(circle5)
  ax.add_artist(circle6)

  # # #plot the radius/fitness of the initial points
  circle1=plt.Circle((spx[0], spy[0]), runFitness[index_min][0], fill=False,color='g')
  circle2=plt.Circle((spx[1], spy[1]), runFitness[index_min][0], fill=False,color='g')
  circle3=plt.Circle((spx[2], spy[2]), runFitness[index_min][0], fill=False,color='g')
  circle4=plt.Circle((spx[3], spy[3]), runFitness[index_min][0], fill=False,color='g')
  circle5=plt.Circle((spx[4], spy[4]), runFitness[index_min][0], fill=False,color='g')
  circle6=plt.Circle((spx[5], spy[5]), runFitness[index_min][0], fill=False,color='g')
  ax.add_artist(circle1)
  ax.add_artist(circle2)
  ax.add_artist(circle3)
  ax.add_artist(circle4)
  ax.add_artist(circle5)
  ax.add_artist(circle6)

  # # # plt.scatter(bestX,bestY,color='red', marker='^')

  plt.title(f'Initial (green) and final locations (red) from \nRSA for 6-center on {currDataset} City data \nmaxiter = {len(runFitness)-1}')

  xmin, xmax, ymin, ymak = plt.axis()
  plt.text(xmax-100, ymin+400,"End Fitness ={:05.3f}".format(runFitness[index_min][len(runFitness)-1]), horizontalalignment='right', color='red') 
  plt.text(xmax-100, ymin+1600,"Init Fitness ={:05.3f}".format(runFitness[index_min][0]), horizontalalignment='right', color='green') 

  plt.xticks(rotation = 35)
  plt.xlabel('x coordinate')
  plt.ylabel('y coordinate')

  # # plt.show
  plt.savefig(f'_Output/{currPC}/{currPC}_6CenterPlot.jpg',dpi=80,bbox_inches='tight')
  plt.clf()

def generateGenerationFitnessPlot(currPC, currDataset):
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  finalFitness=[]
  maxiter=len(runFitness)-1
  for i in range(0,len(runFitness.iloc[0])):
      finalFitness.append(runFitness[i][maxiter])
  x = np.linspace(1, len(finalFitness),len(finalFitness))

  plt.bar(x, finalFitness)
  plt.title(f'Final Fitness Values for 6-center problem on {currDataset} City\nHousehold data using RSA with T = {maxiter}')
  plt.ylabel('Maximum Distance')
  plt.xlabel('Run Number')

  plt.savefig(f'_Output/{currPC}/{currPC}_FinalFitVal.png',dpi=300)
  plt.clf()

def generateRuntime(currPC, currDataset):
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  maxiter=len(runFitness)-1
  runTime=[]
  runTimes = pd.read_csv(f"{currPC}/{currPC}_RSA_RunTime.csv",header=None)
  for i in range(0,len(runTimes)):
      runTime.append(runTimes.iloc[i][0])
  x = np.linspace(1, len(runTimes), len(runTimes))
  plt.bar(x, runTime)
  plt.ylabel('Time (seconds)')
  plt.xlabel('Run Number')
  plt.title(f'Final Runtime for 6-center problem on {currDataset} City \n data using RSA with T = {maxiter}')
  plt.savefig(f'_Output/{currPC}/{currPC}_Runtimes.png',dpi=300)
  plt.clf()

def generateBestFitTime(currPC, currDataset):
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  num_iters=len(runFitness)

  bestFitTime = pd.read_csv(f"{currPC}/{currPC}_RSA_BestFitTime.csv",header=None)
  x = np.linspace(0, num_iters,num_iters)

  for i in range(0,len(bestFitTime.iloc[0])):
      plt.plot(x,bestFitTime[i],'g-')

  plt.title(f"Time when Best Solution is found for the 6-center problem on\n{currDataset} City Household Data Using RSA with configuration \n{currPC}")
  plt.xlabel("Iteration Number")
  plt.ylabel("Time (s)")
  plt.savefig(f'_Output/{currPC}/{currPC}_BestFitTime.png',dpi=300)
  plt.clf()

def generatePhasePerformance(currPC, currDataset):
  highWalking = []
  bellyWalking = []
  huntingCoordination = []
  huntingCooperation = []

  df = pd.read_csv(f'{currPC}/{currPC}_RSA_PhasesPerformance.csv', header=None)

  for i in df.iterrows(): 
    highWalking.append(i[1][0].item())
    bellyWalking.append(i[1][1].item())
    huntingCoordination.append(i[1][2].item())
    huntingCooperation.append(i[1][3].item())

  phasesAverage = []
  phasesAverage.append(sum(highWalking)/len(highWalking))
  phasesAverage.append(sum(bellyWalking)/len(bellyWalking))
  phasesAverage.append(sum(huntingCoordination)/len(huntingCoordination))
  phasesAverage.append(sum(huntingCooperation)/len(huntingCooperation))

  x = ["High\nWalking", "Belly\nWalking", "Hunting\nCoordination", "Hunting\nCooperation"]

  plt.bar(x, phasesAverage)
  plt.ylabel('Number of times')
  plt.xlabel('Phases')
  plt.title(f'Average number of times the solutions changed their values\nin each RSA phase using the RSA configuration {currPC}')
  plt.savefig(f'_Output/{currPC}/{currPC}_PhasePerformance.png',dpi=300)
  plt.clf()

if __name__ == "__main__":
  main()