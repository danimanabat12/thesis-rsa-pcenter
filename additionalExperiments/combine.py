#!/opt/homebrew/bin/python3
import sys
import pandas as pd

def main():
  if (len(sys.argv) < 2 or len(sys.argv) > 2): 
    print("Wrong argument. Try again")
    exit()
  else:
    currPC = sys.argv[1]

  for i in range(1,10):
    interface(currPC+str(i))

def interface(currPC): 
  concatCsvs(currPC, "RSA_AverageFitnesses")
  concatCsvs(currPC, "RSA_BestFitTime", type=1)
  concatCsvs(currPC, "RSA_FinalBestX")
  concatCsvs(currPC, "RSA_FinalBestY")
  concatCsvs(currPC, "RSA_GenerationBestFitnesses", type=1)
  concatCsvs(currPC, "RSA_InitBestX")
  concatCsvs(currPC, "RSA_InitBestY")
  concatCsvs(currPC, "RSA_PhasesPerformance")
  concatCsvs(currPC, "RSA_RunTime")

def concatCsvs(currPC, fileName, type=0):
  file1 = pd.read_csv(f"{currPC}/run_1-10/{currPC}_{fileName}.csv", header=None)
  file2 = pd.read_csv(f"{currPC}/run_11-20/{currPC}_{fileName}.csv", header=None)
  file3 = pd.read_csv(f"{currPC}/run_21-30/{currPC}_{fileName}.csv", header=None)

  if type == 1:
    finaldf = pd.concat([file1, file2, file3], axis=1)
  else: 
    finaldf = pd.concat([file1, file2, file3])
  path = f"{currPC}/{currPC}_{fileName}.csv"
  finaldf.to_csv(path, index=False, header=False)
  with open(path) as f:
    lines = f.readlines()
    last = len(lines) - 1
    lines[last] = lines[last].replace('\r','').replace('\n','')
  with open(path, 'w') as wr:
    wr.writelines(lines)
if __name__ == "__main__":
  main()