# Usage: python generatePlots.py <configuration-folder> <configuration-dataset> <alpha-value> <beta-value>
# Note: Cannot actually be used here since the configuration folders as well as the results are hidden.
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys, os
import tempfile
import shutil

mainFont = {'fontname': 'Times New Roman'}
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

def main():
  if (len(sys.argv) < 5 or len(sys.argv) > 5): 
    print("Wrong argument. Try again")
    exit()
  else:
    currPC = sys.argv[1]
    currDataset = sys.argv[2]
    alpha = sys.argv[3]
    beta = sys.argv[4]

  path = os.getcwd()+f"/_Output/{currPC}"
  if (os.path.exists(path)):
    tmp = tempfile.mktemp(dir=os.path.dirname(path))
    shutil.move(path, tmp)
    shutil.rmtree(tmp)
  os.makedirs(path)
  arcgis = os.getcwd()+f"/_Output/ArcGIS/{currPC}"
  if (os.path.exists(arcgis)):
    tmp = tempfile.mktemp(dir=os.path.dirname(arcgis))
    shutil.move(arcgis, tmp)
    shutil.rmtree(tmp)
  os.makedirs(arcgis)

  generateConvergenceMap(currPC, currDataset, alpha, beta)
  generateDatasetMap(currPC, currDataset, alpha, beta)
  generatePlottedFacilities(currPC, currDataset, alpha, beta)
  generateGenerationFitnessPlot(currPC, currDataset, alpha, beta)
  generateRuntime(currPC, currDataset, alpha, beta)
  generateBestFitTime(currPC, currDataset, alpha, beta)
  generatePhasePerformance(currPC, currDataset, alpha, beta)

def generateConvergenceMap(currPC, currDataset, alpha, beta):
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  runFitness['average'] = runFitness.mean(axis=1)
  num_iters = len(runFitness)
  x = np.linspace(0, num_iters,num_iters)

  for i in range(0,len(runFitness.iloc[0])-1):
      plt.plot(x,runFitness[i],'g-')
  plt.plot(x, runFitness['average'], 'r-')    

  plt.title(f'Convergence Map of 6-center Problem on {currDataset} City Household Data\n for {currPC} ($\\alpha = {alpha}, \\beta$ = {beta}) using T={num_iters-1}')
  plt.xlabel("Iteration Number")
  plt.ylabel("Fitness Value")
  plt.savefig(f'_Output/{currPC}/{currPC}_ConvergenceMap.png')
  plt.clf()

def generateDatasetMap(currPC, currDataset, alpha, beta):
  DX = pd.read_csv(f"{currPC}/{currDataset}_X.csv",header=None, skiprows=1)
  DY = pd.read_csv(f"{currPC}/{currDataset}_Y.csv",header=None, skiprows=1)
  plt.plot(DX, DY, 'x', color='black',markersize=.2);
  plt.title(f'Demand points based on {currDataset} City,\nPhilippines Household Data.')
  plt.xlabel('x-coordinate')
  plt.ylabel('y-coordinate')
  plt.savefig(f'_Output/{currPC}/{currPC}_DatasetMap.png')
  plt.clf()

def generatePlottedFacilities(currPC, currDataset, alpha, beta):
  DX = pd.read_csv(f"{currPC}/{currDataset}_X.csv",header=None, skiprows=1)
  DY = pd.read_csv(f"{currPC}/{currDataset}_Y.csv",header=None, skiprows=1)
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  timePerRun = pd.read_csv(f"{currPC}/{currPC}_RSA_BestFitTime.csv", header=None)
  runTimes = pd.read_csv(f"{currPC}/{currPC}_RSA_RunTime.csv", header=None)
  initXValues = pd.read_csv(f"{currPC}/{currPC}_RSA_InitBestX.csv",header=None)
  initYValues = pd.read_csv(f"{currPC}/{currPC}_RSA_InitBestY.csv",header=None)
  bestXValues = pd.read_csv(f"{currPC}/{currPC}_RSA_FinalBestX.csv",header=None)
  bestYValues = pd.read_csv(f"{currPC}/{currPC}_RSA_FinalBestY.csv",header=None)

  A = runFitness.iloc[len(runFitness)-1]
  B = timePerRun.iloc[-1]
  index_min = np.argmin(A)
  table = pd.concat([A,B,runTimes], axis=1, ignore_index=True)
  std = [table[0].std(), table[1].std(), table[2].std()]
  mean = [table[0].mean(), table[1].mean(), table[2].mean()]
  table.loc[-1] = mean
  table.index = table.index+1
  table.sort_index()
  table.loc[-1] = std
  table.index = table.index+1
  table.sort_index()
  table[0]=['%.4f'%row for row in table[0]]
  table[1]=['%.4f'%row for row in table[1]]
  table[2]=['%.4f'%row for row in table[2]]

  table.to_csv(f'_Output/ArcGIS/{currPC}/{currPC}_table.csv', index=False, mode='w+')
  spx=initXValues.iloc[index_min]
  spy=initYValues.iloc[index_min]
  bestX=bestXValues.iloc[index_min]
  bestY=bestYValues.iloc[index_min]

  init_coords = pd.concat([spx, spy], axis=1)
  init_coords = init_coords.set_axis(['X', 'Y'], axis=1)
  final_coords = pd.concat([bestX, bestY], axis=1)
  final_coords = final_coords.set_axis(['X', 'Y'], axis=1)
  dataset = pd.concat([DX, DY], axis=1)
  dataset = dataset.set_axis(['X', 'Y'], axis=1)

  init_coords.to_csv(f'_Output/ArcGIS/{currPC}/{currPC}_initcoords_{runFitness.iloc[0][index_min]}.csv', index=False, mode='w+')
  final_coords.to_csv(f'_Output/ArcGIS/{currPC}/{currPC}_finalcoords_{runFitness.iloc[len(runFitness)-1][index_min]}.csv', index=False, mode='w+')
  dataset.to_csv(f"_Output/ArcGIS/{currPC}/{currDataset}.csv", index=False, mode='w+')

  ax = plt.gca()

  ax.set_ylim([min(spy)-runFitness[index_min][0]-runFitness[index_min][0]*0.5, max(spy)+runFitness[index_min][0]+runFitness[index_min][0]*0.5])
  ax.set_xlim([min(spx)-runFitness[index_min][0], max(spx)+runFitness[index_min][0]])
  ax.set_aspect('equal', adjustable='datalim')

  # Plotting of demand points
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

  plt.title(f'Initial (green) and Final locations (red) of the 6-center \non {currDataset} City dataset using RSA\'s {currPC} ($\\alpha = {alpha}, \\beta$ = {beta})\nwith T={len(runFitness)-1}')

  xmin, xmax, ymin, ymak = plt.axis()
  ax.text(0.99, 0.005,"End Fitness ={:05.3f}".format(runFitness[index_min][len(runFitness)-1]), transform=ax.transAxes,horizontalalignment='right', color='red', va='bottom') 
  ax.text(0.99, 0.04,"Init Fitness ={:05.3f}".format(runFitness[index_min][0]), transform=ax.transAxes,horizontalalignment='right', color='green', va='bottom') 
  # ax.text(0.99, 0.01,"Init Fitness ={:05.3f}".format(runFitness[index_min][0]), horizontalalignment='right', color='green') 

  plt.xticks(rotation = 35)
  plt.xlabel('x-coordinate')
  plt.ylabel('y-coordinate')

  plt.savefig(f'_Output/{currPC}/{currPC}_6CenterPlot.jpg',dpi=300,bbox_inches='tight')
  plt.clf()

def generateGenerationFitnessPlot(currPC, currDataset, alpha, beta):
  runFitness = pd.read_csv(f"{currPC}/{currPC}_RSA_GenerationBestFitnesses.csv",header=None)
  finalFitness=[]
  maxiter=len(runFitness)-1
  mean = 0
  for i in range(0,len(runFitness.iloc[0])):
      finalFitness.append(runFitness[i][maxiter])
      mean += runFitness[i][maxiter]
  x = np.linspace(1, len(finalFitness),len(finalFitness))
  mean = mean/30

  plt.bar(x, finalFitness)
  plt.axhline(y = mean, color = 'r', linestyle = '-')
  plt.title(f'Final Fitness Values for 6-center problem on {currDataset} City Household data\nusing RSA\'s {currPC} ($\\alpha = {alpha}, \\beta$ = {beta}) with T = {maxiter}')
  plt.ylabel('Maximum Distance')
  plt.xlabel('Run Number')

  plt.savefig(f'_Output/{currPC}/{currPC}_FinalFitVal.png',dpi=300)
  plt.clf()

def generateRuntime(currPC, currDataset, alpha, beta):
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
  plt.title(f'Final Runtime for 6-center problem on {currDataset} City household data\nusing RSA\'s {currPC} ($\\alpha = {alpha}, \\beta$ = {beta}) with T = {maxiter}')
  plt.savefig(f'_Output/{currPC}/{currPC}_Runtimes.png',dpi=300)
  plt.clf()

def generateBestFitTime(currPC, currDataset, alpha, beta):
  bestFitTime = pd.read_csv(f"{currPC}/{currPC}_RSA_BestFitTime.csv",header=None)
  bestFitTime['average'] = bestFitTime.mean(axis=1)
  num_iters=len(bestFitTime)
  x = np.linspace(0, num_iters,num_iters)

  for i in range(0,len(bestFitTime.iloc[0])-1):
      plt.plot(x,bestFitTime[i],'g-')
  plt.plot(x, bestFitTime['average'], 'r-')    

  plt.title(f"Time when best solution is found for the 6-center problem\non {currDataset} City Household Data using RSA\'s {currPC} ($\\alpha = {alpha}, \\beta$ = {beta})")
  plt.xlabel("Iteration Number")
  plt.ylabel("Time (s)")
  plt.tight_layout()
  plt.savefig(f'_Output/{currPC}/{currPC}_BestFitTime.png',dpi=300)
  plt.clf()

def generatePhasePerformance(currPC, currDataset, alpha, beta):
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
  plt.xlabel('Strategy')
  plt.title(f'Average number of times the solutions changed their values\nin each RSA strategy for {currPC} ($\\alpha = {alpha}, \\beta$ = {beta})')
  plt.savefig(f'_Output/{currPC}/{currPC}_PhasePerformance.png',dpi=300)
  plt.clf()

if __name__ == "__main__":
  main()