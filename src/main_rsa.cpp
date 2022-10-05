
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
vector<vector<double> > demandPointsVec;
// Vector containing the best obtained solution so far
vector<vector<double> > bestSoln;
// 
double epsilon = 0.0001;

// Function prototype
void initializeParameters(char* parameterSettingFile);
void initializeDataset (int dataset);
double generateRandom();
void initializeSolution();
void printSolution();
double calculateDistance(double x1, double y1, double x2, double y2);
double calculateFitness(int solutionIndex);

double calculateFitness2(vector<vector<double> > thisSoln);
// Temporary, for checking purposes only.
double calculateFitness3();
void findBest();
double calculateES (double t);
vector<double> calculateEta(int currentDim, vector<double> P);
int generateRandomSolIndex();
vector<double> calculateReduceFunction(int currentDim);
vector<double> calculateXYAverage(int currentSol);
vector<double> calculatePercentageDifference(int currentSol, int currentDim);
void rsaPcenter ();

// Main function
int main(int argc, char** argv) {

		// Read from command line. 
		if (argc < 2 || argc > 2) {
				cout << "Usage: " << argv[0] << " <Parameter-settings.txt>" << endl;
				return 0;
		}
		initializeParameters(argv[1]);
		initializeDataset(dataset);
		for (int i = 0; i < numberOfRuns; i++) {
			// calculateFitness3();
			initializeSolution();
			printSolution();
			// Moved find best sa sulod sa rsaPcenter(). 
			// findBest();
			rsaPcenter();
			// findBest();
			candidateSolution.clear();
			candidateSolution.resize(numberOfReptiles); 
			for(int i = 0; i < candidateSolution.size(); i++) {
					candidateSolution[i].resize(facilities);
			}
		}
}

// Location na lang ang ipasa. Conditional statement for X and Y
// Assignment: change initialize to read
void initializeDataset (int dataset) {
	string parameterSettingFile = (dataset == 1 ? "Davao" : (dataset == 2 ? "Digos" : "Tagum"));
	// Anything related to file path, dili mag work. I-setup ang code such that all files that will be called by your source code will be placed on the directory where you will invoke the executable file. 
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
		// printf("[%f, %f]\n", temp.at(0), temp.at(1));
    // Find upper and lowerbound: X-coordinate
    upperBoundX = (temp.at(0) > upperBoundX) ? temp.at(0) : upperBoundX; 
    lowerBoundX = (temp.at(0) < lowerBoundY) ? temp.at(0) : lowerBoundX;
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

		cout << "\nPARAMETERS\n============================"<< endl;
		cout << "Number of reptiles (solution): " << numberOfReptiles << endl;
		cout << "Number of facilities (dimension): " << facilities << endl;
		cout << "Max iteration: " << T << endl;
		cout << "Number of runs: " << numberOfRuns << endl;
		cout << "Alpha: " << alpha << endl;
		cout << "Beta: " << beta << endl;
		cout << "Dataset: " << (dataset == 1 ? "Davao" : (dataset == 2 ? "Digos" : "Tagum")) << endl;
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
						candidateSolution[i][j].push_back(generateRandom() * (upperBoundX - lowerBoundX) + lowerBoundX);
						candidateSolution[i][j].push_back(generateRandom() * (upperBoundY - lowerBoundY) + lowerBoundY);
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
		// inner min - Buhaton is: per demand point, pangitaon nato ang nearest facility sa iyaha given a set of facility (located by index)
		// max - Pangitaon nato ang demand point with the farthest minimum distance jud. 
				// Approach: 
						// 1. Mag for loop ko, iterating each demand point, then i-calculate nako ang distance niya with all the given position sa solution
						// 2. Pangitaon nako ang min distance with its corresponding position index (facility)
						// 3. After ma-determine ang min distance, i-compare it with max para isahan na lang
						// 4. Return max as the fitness value.
						// Note: no need na mag store sa list ug i-take note ang index kay ang concern ra man is ang max jud which we can already optimize naman
		for (int i = 0; i < demandPointsVec.size(); i++) {
				// Initialize min as the distance between the current demand point and the first facility coordinate in the solution.
				// cout << "Facility 1" << endl;
				min = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], candidateSolution[solutionIndex][0].at(0), candidateSolution[solutionIndex][0].at(1));
				// printf("%f\n", min); //CONTINUE HERE!!!!
				for (int j = 1; j < facilities; j++) {
						// cout << "Facility " << j+1 << endl;
						distance = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], candidateSolution[solutionIndex][j].at(0), candidateSolution[solutionIndex][j].at(1));
						// printf("%f\n", distance);
						if (min > distance) {
								min = distance; 
								// printf("MINIMUM HAS BEEN REPLACED, NEW VAL: %f\n", min);
						}
				}
				// printf("MIN: %f\n", min);
				// Max holds the maximum of all distances from a demand point to a facility.
				if (min > max) {
					max = min; 
					// printf("MAX HAS BEEN REPLACED, NEW VAL: %f\n", max);
				}
		}
		return max;
		// outer min - dili na siya included sa calculate fitness, pero bale pangitaon nato ang set of facility nga naka-produce og smallest distance/fitness value.
}

// Another variant of calculateFitness where the parameter is a vector instead of an index.
double calculateFitness2(vector<vector<double> > thisSoln) {
		double distance;
		double min; 
		double max = 0;
		// inner min - Buhaton is: per demand point, pangitaon nato ang nearest facility sa iyaha given a set of facility (located by index)
		// max - Pangitaon nato ang demand point with the farthest minimum distance jud. 
				// Approach: 
						// 1. Mag for loop ko, iterating each demand point, then i-calculate nako ang distance niya with all the given position sa solution
						// 2. Pangitaon nako ang min distance with its corresponding position index (facility)
						// 3. After ma-determine ang min distance, i-compare it with max para isahan na lang
						// 4. Return max as the fitness value.
						// Note: no need na mag store sa list ug i-take note ang index kay ang concern ra man is ang max jud which we can already optimize naman
		for (int i = 0; i < demandPointsVec.size(); i++) {
				// Initialize min as the distance between the current demand point and the first facility coordinate in the solution.
				// cout << "Facility 1" << endl;
				min = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], thisSoln[0].at(0), thisSoln[0].at(1));
				// printf("%f\n", min); //CONTINUE HERE!!!!
				for (int j = 1; j < facilities; j++) {
						// cout << "Facility " << j+1 << endl;
						distance = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], thisSoln[j].at(0), thisSoln[j].at(1));
						// printf("%f\n", distance);
						if (min > distance) {
								min = distance; 
								// printf("MINIMUM HAS BEEN REPLACED, NEW VAL: %f\n", min);
						}
				}
				// printf("MIN: %f\n", min);
				// Max holds the maximum of all distances from a demand point to a facility.
				if (min > max) {
					max = min; 
					// printf("MAX HAS BEEN REPLACED, NEW VAL: %f\n", max);
				}
		}
		return max;
		// outer min - dili na siya included sa calculate fitness, pero bale pangitaon nato ang set of facility nga naka-produce og smallest distance/fitness value.
}

// Temporary, delete function once done.
double calculateFitness3() {
		double distance;
		double min; 
		double max = 0;
		double tempDemandPoints[10][2] = {10,76,30,76,40,70,21,24,65,90,38,19,20,70,13,65,10,63,60,55}; 
		double facility[3][2]={20,80,59,75,48,21};

		cout << "FACILITY[7]: " << facility[0][0] << " " << facility[0][1] << endl;
		for (int i = 0; i < 10; i++) {
				min = calculateDistance(tempDemandPoints[i][0], tempDemandPoints[i][1], facility[0][0], facility[0][1]);
				for (int j = 1; j < 3; j++) {
						distance = calculateDistance(tempDemandPoints[i][0], tempDemandPoints[i][1], facility[j][0], facility[j][1]);
						if (min > distance) {
								min = distance; 
						}
				}
				printf("MIN: %f\n", min);
				if (min > max) max = min; 
		}
		printf("%f", max);
		return max;
}

//I edit such that mag match siya sa pseudocode.
void findBest() {
		double bestFitness; 
		int bestFitnessSolutionIndex;
		double currentFitness;
		// Diri na tong outer min. 
				// Approach: 
				// 1. Basically, I will iterate through the solution vector, tapos ipasa ko yung index sa calculateFitness na function (index lang kay global naman ang vector)
				// 2. Take note the best solution (index na lang) and best fitness.
		// Initialize bestFitness with the fitness value of the first solution.
		bestFitness = calculateFitness(0);
		bestFitnessSolutionIndex = 0;
		// cout << "\n===================================" << endl;
		printf("FITNESS VALUE OF REPTILE[1]: %f\n", bestFitness);
		for (int i = 1; i < numberOfReptiles; i++) {
				currentFitness = calculateFitness(i);
				printf("FITNESS VALUE OF REPTILE[%i]: %f\n", i+1, currentFitness);
				if (bestFitness > currentFitness) {
						bestFitness = currentFitness; 
						bestFitnessSolutionIndex = i;
						// cout << "BEST FITNESS UPDATED" << endl;
				}
		}
		// cout << "\n===========" << endl; 
		// cout << "BEST FITNESS VALUE: ";
		// printf("%f \n", bestFitness);
		// cout << "BEST SOLUTION: ";
		// cout << "rPn[" << bestFitnessSolutionIndex+1 << "] = { ";  
		for (int j = 0; j < facilities; j++)  {
				// printf("(%f, %f) ", candidateSolution[bestFitnessSolutionIndex][j].at(0), candidateSolution[bestFitnessSolutionIndex][j].at(1));
        bestSoln[j].push_back(candidateSolution[bestFitnessSolutionIndex][j].at(0));
        bestSoln[j].push_back(candidateSolution[bestFitnessSolutionIndex][j].at(1));
		}
		cout << "}\n" << endl;
}

void rsaPcenter () {
		cout << "\nRSA START!" << endl;
		cout << "========================================" << endl;
		int t = 1, r1; 
		double ES;
		vector<double> eta, R, P;
		vector<vector <double > > xNew; 
		vector<double > generationBestFitnessValues;
		double generationBestFitness;
		double xNewBestFitness;
		xNew.resize(facilities);

		while (t <= T) {
				// 1. Calculate Fitness function of each candidate solution
				// 2. Find best solution and its fitness
				findBest();
				// 3. Calculate ES
				ES = calculateES(t);
				for (int i = 0; i < numberOfReptiles; i++) {
						cout << "Reptile #" << i+1 << endl;
						xNew.clear();
						xNew.resize(facilities);
						for (int j = 0; j < facilities; j++) {
								cout << "Facility #" << j+1 << endl;
								// 4. Calculate hunting op, reduce function, percentage diff
								R = calculateReduceFunction(j);
								P = calculatePercentageDifference(i, j);
								eta = calculateEta(j, P);
								// cout << setprecision(0) << P.at(0) << " " << P.at(1) << endl;
								printf("R: %f, %f\n, P: %e, %e\n, Eta: %f, %f\n", R.at(0), R.at(1), P.at(0), P.at(1), eta.at(0), eta.at(1));

								// CONTINUE HERE!!!
								// High walking
								if (t <= T/4) {
										double var1 = bestSoln[j].at(0);
										double var2 = eta.at(0);
										double var3 =  beta; 
										double var4 = R.at(0);
										double var5 = generateRandom();
										xNew[j].push_back(bestSoln[j].at(0) - eta.at(0) * beta - R.at(0) * var5);
										xNew[j].push_back(bestSoln[j].at(1) - eta.at(1) * beta - R.at(1) * var5);
								} 
								// Belly walking
								else if (t <= (2*T/4) && t > T/4) {
										r1 = generateRandomSolIndex();
										double var = generateRandom(); 
										xNew[j].push_back(bestSoln[j].at(0) * candidateSolution[r1][j].at(0) * ES * var);
										xNew[j].push_back(bestSoln[j].at(1) * candidateSolution[r1][j].at(1) * ES * var);
								}
								// Hunting coordination
								else if (t <= (3*T/4) && t > (2*T/4)) {
										double var = generateRandom(); 
										xNew[j].push_back(bestSoln[j].at(0) * P.at(0) * var);
										xNew[j].push_back(bestSoln[j].at(0) * P.at(1) * var);
								} 
								// Hunting cooperation
								else if (t <= T && t > (3*T/4)) {
										double var = generateRandom(); 
										xNew[j].push_back(bestSoln[j].at(0) - eta.at(0) * epsilon - R.at(0) * var);
										xNew[j].push_back(bestSoln[j].at(1) - eta.at(1) * epsilon - R.at(1) * var);
								}
								// cout << t << " " << i << " " << j << endl;
								R.clear();
								P.clear();
								eta.clear();
								// cout << "CLEARED" << endl;
						}
						generationBestFitness = calculateFitness2(bestSoln);
						xNewBestFitness = calculateFitness2(xNew);
						if (calculateFitness(i) > xNewBestFitness) {
							candidateSolution[i] = xNew;				
							if (generationBestFitness > xNewBestFitness) {
								bestSoln = xNew;
								generationBestFitness = xNewBestFitness;
							}
						}
				}
				generationBestFitnessValues.push_back(generationBestFitness);
				t++;
		}

		// cout << "GENERATION BEST FITNESSES" << endl;
		// for (int k = 0; k < T; k++) {
		// 	printf("%d: %f\n", k, generationBestFitnessValues.at(k));
		// }
}

double calculateES (double t) {
		double ES;
		random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> dist(-1, 1);
		ES = 2 * dist(mt) * (1 -  (t/T)); 
		return ES;
}

vector<double> calculateEta(int currentDim, vector<double> P) {
		vector<double> eta;
		eta.push_back(bestSoln[currentDim].at(0) * P.at(0));
		eta.push_back(bestSoln[currentDim].at(1) * P.at(1));
		return eta;
}

int generateRandomSolIndex() {
		// Change rand
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
		// Overhaul
		// Evaluate variables:
		// P is a vector jud. So is M. Kay duha man ang value (per coordinate)
		// I also think x_ij kay vector pud dapat, so with Bestj
		vector<double> P, x_ij, M, bestJ;
    x_ij.push_back(candidateSolution[currentSol][currentDim].at(0)); 
		x_ij.push_back(candidateSolution[currentSol][currentDim].at(1)); 
		// printf("X_ij: (%f, %f)\n", x_ij.at(0), x_ij.at(1));
    M = calculateXYAverage(currentSol);
		printf("M: %f, %f\n", M.at(0), M.at(1));
    bestJ.push_back(bestSoln[currentDim].at(0));
		bestJ.push_back(bestSoln[currentDim].at(1));
		// printf("Best_J: (%f, %f)\n", bestJ.at(0), bestJ.at(1));
		double numeratorX = (x_ij.at(0) - M.at(0));
		double numeratorY = (x_ij.at(1) - M.at(1));
		double denominatorX = ((bestJ.at(0) * (upperBoundX - lowerBoundX)) + epsilon);
		double denominatorY = ((bestJ.at(1) * (upperBoundY - lowerBoundY)) + epsilon);
		printf("Bounds of X: %f, %f\n", lowerBoundX, upperBoundX);
		printf("Bounds of Y: %f, %f\n", lowerBoundY, upperBoundY);
		printf("Best: %f, %f\n", bestJ.at(0), bestJ.at(1));
		printf("Numerator: %f, %f Denominator: %f, %f\n", numeratorX, numeratorY, denominatorX, denominatorY);
    P.push_back(alpha + ((x_ij.at(0) - M.at(0)) / ((bestJ.at(0) * (upperBoundX - lowerBoundX)) + epsilon)));
		P.push_back(alpha + ((x_ij.at(1) - M.at(1)) / ((bestJ.at(1) * (upperBoundY - lowerBoundY)) + epsilon)));
		// printf("Numerator: %f\n", numerator);
		// printf("Denominator: %f\n", denominator);
		// printf("Epsilon: %f\n", std::numeric_limits<double>::epsilon());
		// printf("X Bounds: [%f, %f]\n", lowerBoundX, upperBoundX);
		// printf("Y Bounds: [%f, %f]\n", lowerBoundY, upperBoundY);
		// printf("M: (%f, %f)\n", M.at(0), M.at(1));
		// printf("P: (%f, %f)\n", P.at(0), P.at(1));
		return P;
}