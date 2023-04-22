
#include "headers.h"

using namespace std; 

// Extracted from a file.
int numberOfReptiles; 
int facilities; 
int numberOfRuns;
double upperBoundX;
double lowerBoundX; 
double upperBoundY;
double lowerBoundY;
float alpha;
float beta;
int dataset;
int T;
// Candidate solution vector
vector < vector< vector<double> > > candidateSolution;
// Demand points vector 
vector < vector< double> > demandPointsVec;
// Vector containing the best obtained solution so far
vector < vector< double> > bestSoln;
// Epsilon value
double epsilon = 0.0001;
// Denotes the current run
int currRun = 0;

// Save data mechanism
string paramsFileName;
vector <double> runTime;
vector <double> averageFitnesses;
vector < vector< double> > bestFitTime;
vector < vector< double> > generationBestFitnesses;
vector < vector< vector <double> > > initBests;
vector < vector< vector <double> > > finalBests;
vector < vector< int> > rsaPhasesIncrementation;


// Function prototype
void initializeParameters(char* parameterSettingFile);
void initializeDataset (int dataset);
double generateRandom();
void initializeSolution();
void printSolution();
double calculateDistance(double x1, double y1, double x2, double y2);
double calculateFitness(int solutionIndex);
double calculateFitness2(vector<vector<double> > thisSoln);
void findBest();
double calculateES (double t);
vector<double> calculateEta(int currentDim, vector<double> P);
int generateRandomSolIndex();
vector<double> calculateReduceFunction(int currentDim);
vector<double> calculateXYAverage(int currentSol);
vector<double> calculatePercentageDifference(int currentSol, int currentDim);
void rsaPcenter();
void resultsToCsv();

// Main function
int main(int argc, char** argv) {
  // Read from command line. 
  if (argc < 2 || argc > 2) {
    cout << "Usage: " << argv[0] << " <Parameter-settings.txt>" << endl;
    return 0;
  }
  
  paramsFileName = argv[1];
  paramsFileName = paramsFileName.substr(0, paramsFileName.find('.'));
  cout << paramsFileName << endl;
  // Initialize parameters and datasets
  initializeParameters(argv[1]);
  initializeDataset(dataset);
  
  for (int i = 1; i <= numberOfRuns; i++) {
    cout << "RUN " << i << endl;
    // Start time per run
    auto start = chrono::steady_clock::now();
    initializeSolution();
    rsaPcenter();
    auto end = chrono::steady_clock::now();
    auto duration = end - start;
    candidateSolution.clear();
    candidateSolution.resize(numberOfReptiles); 
    for(int j = 0; j < candidateSolution.size(); j++) {
        candidateSolution[j].resize(facilities);
    }
    currRun++;
    runTime.push_back(chrono::duration <double> (duration).count());
  }

  resultsToCsv();
}

void initializeDataset (int dataset) {
	string parameterSettingFile = (dataset == 1 ? "Davao" : (dataset == 2 ? "Digos" : "Tagum"));
  string demandPointX = parameterSettingFile + "_X.csv";
  string demandPointY = parameterSettingFile + "_Y.csv";
  vector<double> temp;
  ifstream datasetX;
  ifstream datasetY; 
  string line, word, extractedX, extractedY; 
  datasetX.open(demandPointX, ios::in);
  datasetY.open(demandPointY, ios::in);

  // Remove header line 
  getline(datasetX, extractedX, '\n'); 
  getline(datasetY, extractedY, '\n'); 

  // First iteration sa labas so that we can initialize upperbound and lowerbound or each coordinate sa value ng first iteration.
  getline(datasetX, extractedX, '\n'); 
  getline(datasetY, extractedY, '\n'); 
  temp.push_back(atof(extractedX.c_str()));
  temp.push_back(atof(extractedY.c_str()));
	// printf("[%f, %f]\n", temp.at(0), temp.at(1));
  demandPointsVec.push_back(temp);
  upperBoundX = temp.at(0);
  lowerBoundX = temp.at(0);
  upperBoundY = temp.at(1);
  lowerBoundY = temp.at(1);
  temp.clear();
  // Assuming that both file has the same number of values.
  while (!datasetX.eof()) {
    getline(datasetX, extractedX, '\n'); 
    getline(datasetY, extractedY, '\n'); 
    temp.push_back(atof(extractedX.c_str()));
    temp.push_back(atof(extractedY.c_str()));
    demandPointsVec.push_back(temp);

    upperBoundX = (temp.at(0) > upperBoundX) ? temp.at(0) : upperBoundX; 
    lowerBoundX = (temp.at(0) < lowerBoundX) ? temp.at(0) : lowerBoundX;
    upperBoundY = (temp.at(1) > upperBoundY) ? temp.at(1) : upperBoundY; 
    lowerBoundY = (temp.at(1) < lowerBoundY) ? temp.at(1) : lowerBoundY;
    temp.clear();
  }
  datasetX.close();
  datasetY.close();
  printf("X-coordinate bounds: [%f, %f]\n", lowerBoundX, upperBoundX);
  printf("Y-coordinate bounds: [%f, %f]\n", lowerBoundY, upperBoundY);
}

void initializeParameters(char* parameterSettingFile) {
  fstream params; 
  string fileDirectory = string(parameterSettingFile);
  string line, word; 
  vector<double> fileContent;
  params.open(fileDirectory, ios::in);
  if (params.is_open()) {
    while(getline(params, line)) {
      stringstream s(line);
      while(getline(s, word, ' ')) {}
      fileContent.push_back(stod(word));    
    }
    params.close();
  }
  // parameter file should follow a certain ordering.
  numberOfReptiles = int(fileContent.at(0)); 
  facilities = int(fileContent.at(1)); 
  T = fileContent.at(2);
  numberOfRuns = fileContent.at(3);
  alpha = fileContent.at(4);
  beta = fileContent.at(5);
  dataset = fileContent.at(6);
  // Resize vector to numberOfReptiles and facilities.
  candidateSolution.resize(numberOfReptiles); 
  for(int i = 0; i < candidateSolution.size(); i++) {
    candidateSolution[i].resize(facilities);
  }
  bestSoln.resize(facilities);
}

double generateRandom() {
	random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> dist(0.0, 1.0+epsilon);
	return dist(mt);
}

void initializeSolution() {
  for (int i = 0; i < numberOfReptiles; i++) {
    for (int j = 0; j < facilities; j++) {
      double x = generateRandom() * (upperBoundX - lowerBoundX) + lowerBoundX;
      double y = generateRandom() * (upperBoundY - lowerBoundY) + lowerBoundY;
      candidateSolution[i][j].push_back(x);
      candidateSolution[i][j].push_back(y);
    }
  }
}

void printSolution() {
  cout << "\nCANDIDATE SOLUTION MATRIX\n==============================" << endl;
  for (int i = 0; i < numberOfReptiles; i++) {
    cout << "repPopn[" << i+1 << "] = {";
    for (int j = 0; j < facilities; j++)  {
      printf("(%f, %f) ", candidateSolution[i][j].at(0), candidateSolution[i][j].at(1));
    }
    cout << "}" << endl;
  }
}

double calculateDistance(double x1, double y1, double x2, double y2) {
  return sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
}

double calculateFitness(int solutionIndex) {
  double distance;
  double min; 
  double max = 0;

  for (int i = 0; i < demandPointsVec.size(); i++) {
    min = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], candidateSolution[solutionIndex][0].at(0), candidateSolution[solutionIndex][0].at(1));
    for (int j = 1; j < facilities; j++) {
      distance = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], candidateSolution[solutionIndex][j].at(0), candidateSolution[solutionIndex][j].at(1));
      if (min > distance) {
        min = distance; 
      }
    }
    if (min > max) {
      max = min; 
    }
  }
  return max;
}

// Another variant of calculateFitness where the parameter is a vector instead of an index.
double calculateFitness2(vector<vector<double> > thisSoln) {
  double distance;
  double min; 
  double max = 0;
  for (int i = 0; i < demandPointsVec.size(); i++) {
    min = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], thisSoln[0].at(0), thisSoln[0].at(1));
    for (int j = 1; j < facilities; j++) {
      distance = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], thisSoln[j].at(0), thisSoln[j].at(1));
      if (min > distance) {
        min = distance; 
      }
    }
    if (min > max) {
      max = min; 
    }
  }
  return max;
}

void findBest() {
  double bestFitness; 
  int bestFitnessSolutionIndex;
  double currentFitness;

  bestFitness = calculateFitness(0);
  bestFitnessSolutionIndex = 0;

  for (int i = 1; i < numberOfReptiles; i++) {
    currentFitness = calculateFitness(i);
    if (bestFitness > currentFitness) {
      bestFitness = currentFitness; 
      bestFitnessSolutionIndex = i;
    }
  }
  bestSoln.clear();
  bestSoln.resize(facilities);
  for (int j = 0; j < facilities; j++)  {
    bestSoln[j].push_back(candidateSolution[bestFitnessSolutionIndex][j].at(0));
    bestSoln[j].push_back(candidateSolution[bestFitnessSolutionIndex][j].at(1));
  }
}

void rsaPcenter() {
  int t = 0, r1; 
  double ES, averageFitness = 0;
  vector<double> eta, R, P;
  vector<vector <double > > xNew; 
  vector<double> generationBestFitnessValues;
  vector<double> thisBestFitTime;
  vector<int> phasesPerformance(4, 0);
  
  double generationBestFitness;
  double xNewBestFitness;
  xNew.resize(facilities);

  auto start = chrono::steady_clock::now();
  findBest();
  auto end = chrono::steady_clock::now();
  double timestamp = chrono::duration <double> (end-start).count();

  initBests.push_back(bestSoln);

  cout << "Generation ";
  while (t <= T) {
    cout << t << " ";
    cout.flush();
    findBest();
    ES = calculateES(t);
    generationBestFitness = calculateFitness2(bestSoln);

    for (int i = 0; i < numberOfReptiles; i++) {
      xNew.clear();
      xNew.resize(facilities);
      for (int j = 0; j < facilities; j++) {
        // 4. Calculate hunting op, reduce function, percentage diff
        R = calculateReduceFunction(j);
        P = calculatePercentageDifference(i, j);
        eta = calculateEta(j, P);

        // High walking
        if (t <= T/4) {
          xNew[j].push_back(bestSoln[j].at(0) * -eta.at(0) * beta - R.at(0) * generateRandom());
          xNew[j].push_back(bestSoln[j].at(1) * -eta.at(1) * beta - R.at(1) * generateRandom());
        } 
        // Belly walking
        else if (t <= (2*T/4) && t > T/4) {
          r1 = generateRandomSolIndex();
          xNew[j].push_back(bestSoln[j].at(0) * candidateSolution[r1][j].at(0) * ES * generateRandom());
          xNew[j].push_back(bestSoln[j].at(1) * candidateSolution[r1][j].at(1) * ES * generateRandom());
        }
        // Hunting coordination
        else if (t <= (3*T/4) && t > (2*T/4)) {
          xNew[j].push_back(bestSoln[j].at(0) * P.at(0) * generateRandom());
          xNew[j].push_back(bestSoln[j].at(0) * P.at(1) * generateRandom());
        } 
        // Hunting cooperation
        else if (t <= T && t > (3*T/4)) {
          double var = generateRandom(); 
          xNew[j].push_back(bestSoln[j].at(0) - eta.at(0) * epsilon - R.at(0) * generateRandom());
          xNew[j].push_back(bestSoln[j].at(1) - eta.at(1) * epsilon - R.at(1) * generateRandom());
        }
        R.clear();
        P.clear();
        eta.clear();
      }
      xNewBestFitness = calculateFitness2(xNew);
      if (calculateFitness(i) > xNewBestFitness) {
        // Replace itself if new solution is better.
        if (t <= T/4) phasesPerformance.at(0) += 1;
        else if (t <= (2*T/4) && t > T/4) phasesPerformance.at(1) += 1;
        else if (t <= (3*T/4) && t > (2*T/4)) phasesPerformance.at(2) += 1;
        else if (t <= T && t > (3*T/4)) phasesPerformance.at(3) += 1;
        candidateSolution[i].clear();
        candidateSolution[i].resize(facilities);
        candidateSolution[i] = xNew;			
        // Replace current generation best fitness if new reptile if new solution is better 
        if (generationBestFitness > xNewBestFitness) {
          end = chrono::steady_clock::now();
          timestamp = chrono::duration <double> (end-start).count();
          bestSoln.clear();
          bestSoln.resize(facilities);
          bestSoln = xNew;
          generationBestFitness = xNewBestFitness;
        }
      }
      // cout << "Generation best fitness: " << generationBestFitness << endl;
    }
    thisBestFitTime.push_back(timestamp);
    generationBestFitnessValues.push_back(generationBestFitness);
    averageFitness += generationBestFitness;
    t++;
  }
  cout << "\n";
  bestFitTime.push_back(thisBestFitTime);
  generationBestFitnesses.push_back(generationBestFitnessValues);
  finalBests.push_back(bestSoln);
  averageFitnesses.push_back(averageFitness/(T+1));
  rsaPhasesIncrementation.push_back(phasesPerformance);
} 

double calculateES (double t) {
  double ES;
  random_device rd;
  mt19937 mt(rd());
  uniform_int_distribution<int> dist(-1, 1);
  ES = 2 * dist(mt) * (1 - (1/T)); 
  return ES;
}

vector<double> calculateEta(int currentDim, vector<double> P) {
  vector<double> eta;
  eta.push_back(bestSoln[currentDim].at(0) * P.at(0));
  eta.push_back(bestSoln[currentDim].at(1) * P.at(1));
  return eta;
}

int generateRandomSolIndex() {
  random_device rd;
  mt19937 mt(rd());
  uniform_int_distribution<int> dist(0, numberOfReptiles-1);
  return dist(mt);
}

vector<double> calculateReduceFunction(int currentDim) {
  int r2 = generateRandomSolIndex();
  vector<double> R;
  R.push_back((bestSoln[currentDim].at(0) - candidateSolution[r2][currentDim].at(0)) / (bestSoln[currentDim].at(0) + epsilon));
  R.push_back((bestSoln[currentDim].at(1) - candidateSolution[r2][currentDim].at(1)) / (bestSoln[currentDim].at(1) + epsilon));
  return R;
}

vector<double> calculateXYAverage(int currentSol) {
  vector<double> M;
  double xAve = 0, yAve = 0;
  for (int i = 0; i < facilities; i++) {
    xAve += candidateSolution[currentSol][i].at(0);
    yAve += candidateSolution[currentSol][i].at(1);
  }
  M.push_back(xAve/facilities);
  M.push_back(yAve/facilities);
  return M;
}

vector<double> calculatePercentageDifference(int currentSol, int currentDim) {
  vector<double> P, M;
  M = calculateXYAverage(currentSol);
  P.push_back(alpha + ((candidateSolution[currentSol][currentDim].at(0) - M.at(0)) / ((bestSoln[currentDim].at(0) * (upperBoundX - lowerBoundX)) + epsilon)));
  P.push_back(alpha + ((candidateSolution[currentSol][currentDim].at(1) - M.at(1)) / ((bestSoln[currentDim].at(1) * (upperBoundY - lowerBoundY)) + epsilon)));
  return P;
}

void resultsToCsv() {
	ofstream runtime(paramsFileName + "_RSA_RunTime.csv");
	ofstream averagefitness(paramsFileName + "_RSA_AverageFitnesses.csv");
	ofstream phasesPerformance(paramsFileName + "_RSA_PhasesPerformance.csv");
	ofstream initBestX(paramsFileName + "_RSA_InitBestX.csv");
	ofstream initBestY(paramsFileName + "_RSA_InitBestY.csv");
	ofstream finalBestX(paramsFileName + "_RSA_FinalBestX.csv");
	ofstream finalBestY(paramsFileName + "_RSA_FinalBestY.csv");
	ofstream bestTime(paramsFileName + "_RSA_BestFitTime.csv");
	ofstream bestFitnesses(paramsFileName + "_RSA_GenerationBestFitnesses.csv");

	// Output file for run time, initial best solutions, and final best solutions
	for (int i = 0; i < runTime.size(); i++) {
		runtime << runTime.at(i);
		averagefitness << averageFitnesses.at(i);
		for (int j = 0; j < finalBests.at(i).size(); j++) {
      initBestX << fixed << initBests[i][j][0];
      initBestY << fixed << initBests[i][j][1];
      finalBestX << fixed << finalBests[i][j][0];
      finalBestY << fixed << finalBests[i][j][1];
      if (j != finalBests.at(i).size()-1) {
        initBestX << ",";
        initBestY << ",";
        finalBestX << ",";
        finalBestY << ",";
      }
		}

		for (int j = 0; j < rsaPhasesIncrementation.at(i).size(); j++) {
			phasesPerformance << rsaPhasesIncrementation.at(i).at(j);
			if (j != rsaPhasesIncrementation.at(i).size()-1) {
				phasesPerformance << ",";
			}
		}

		if (i != runTime.size()-1) {
			initBestX << "\n";
			initBestY << "\n";
			finalBestX << "\n";
			finalBestY << "\n";
			runtime << "\n";
			phasesPerformance << "\n";
			averagefitness << "\n";
		}
	}

	// Output file for generation best fitnesses and best fitness time
	for (int i = 0; i < bestFitTime.at(0).size(); i++) {
		for (int j = 0; j < bestFitTime.size(); j++) {
			bestTime << bestFitTime[j][i];
			bestFitnesses << generationBestFitnesses[j][i];
			// Resolved the error from the last meeting where the first value repeats. 
			if (j != bestFitTime.size()-1) {
				bestTime << ",";
				bestFitnesses << ",";
			}
		}
		if (i != bestFitTime.at(0).size() - 1) {
			bestTime << "\n";
			bestFitnesses << "\n";
		}
	}
}