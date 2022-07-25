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

// Location na lang ang ipasa. Conditional statement for X and Y
void initializeDataset (string parameterSettingFile) {
  string demandPointX = "../p-center-datasets/" + parameterSettingFile + "_X.csv";
  string demandPointY = "../p-center-datasets/" + parameterSettingFile + "_Y.csv";
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
		string fileDirectory = "../parameter-setups/" + string(parameterSettingFile);
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

		cout << "\nPARAMETERS\n============================"<< endl;
		cout << "Number of reptiles (solution): " << numberOfReptiles << endl;
		cout << "Number of facilities (dimension): " << facilities << endl;
		cout << "Max iteration: " << T << endl;
		cout << "Number of runs: " << numberOfRuns << endl;
		cout << "Alpha: " << alpha << endl;
		cout << "Beta: " << beta << endl;
		cout << "Dataset: " << (dataset == 1 ? "Davao" : (dataset == 2 ? "Digos" : "Tagum")) << endl;

		// NOTE: Call initializeDataset afterwards
		initializeDataset((dataset == 1 ? "Davao" : (dataset == 2 ? "Digos" : "Tagum")));
}

double generateRandom() {
		return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

void initializeSolution() {
		srand(time(NULL));
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
		// max - Pangitaon nato ang demand point with the farthest maximum distance jud. 
				// Approach: 
						// 1. Mag for loop ko, iterating each demand point, then i-calculate nako ang distance niya with all the given position sa solution
						// 2. Pangitaon nako ang min distance with its corresponding position index (facility)
						// 3. After ma-determine ang min distance, i-compare it with max para isahan na lang
						// 4. Return max as the fitness value.
						// Note: no need na mag store sa list ug i-take note ang index kay ang concern ra man is ang max jud which we can already optimize naman
		for (int i = 0; i < demandPointsVec.size(); i++) {
				// Initialize min as the distance between the current demand point and the first facility coordinate in the solution.
				min = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], candidateSolution[solutionIndex][0].at(0), candidateSolution[solutionIndex][0].at(1));
				for (int j = 1; j < facilities; j++) {
						distance = calculateDistance(demandPointsVec[i][0], demandPointsVec[i][1], candidateSolution[solutionIndex][j].at(0), candidateSolution[solutionIndex][j].at(1));
						if (min > distance) {
								min = distance; 
						}
				}
				if (min > max) max = min; 
		}
		return max;
		// outer min - dili na siya included sa calculate fitness, pero bale pangitaon nato ang set of facility nga naka-produce og smallest distance/fitness value.
}

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
		cout << "\n===================================" << endl;
		printf("FITNESS VALUE OF REPTILE[1]: %f\n", bestFitness);
		for (int i = 1; i < numberOfReptiles; i++) {
				currentFitness = calculateFitness(i);
				printf("FITNESS VALUE OF REPTILE[%i]: %f\n", i+1, currentFitness);
				if (bestFitness > currentFitness) {
						bestFitness = currentFitness; 
						bestFitnessSolutionIndex = i;
				}
		}
		cout << "\n===========" << endl; 
		cout << "BEST FITNESS VALUE: ";
		printf("%f \n", bestFitness);
		cout << "BEST SOLUTION: ";
		cout << "rPn[" << bestFitnessSolutionIndex+1 << "] = { ";  
		for (int j = 0; j < facilities; j++)  {
				printf("(%f, %f) ", candidateSolution[bestFitnessSolutionIndex][j].at(0), candidateSolution[bestFitnessSolutionIndex][j].at(1));
		}
		cout << "}\n" << endl;
}

double calculateES () {
		double ES; 
		return ES;
}

double calculateEta() {
		double eta;
		return eta;
}

double calculateReduceFunction() {
		double R;
		return R;
}

double calculatePercentageDifference() {
		double P;
		return P;
}

void rsaPcenter () {
		int t = 1; 
		double ES, eta, R, P; 

		while (t <= T) {
				// 1. Calculate Fitness function of each candidate solution
				// 2. Find best solution and its fitness
				// 3. Calculate ES
				ES = calculateES();
				for (int i = 0; i < numberOfReptiles; i++) {
						for (int j = 0; j < facilities; j++) {
								// 4. Calculate hunting op, reduce function, percentage diff
								eta = calculateEta();
								R = calculateReduceFunction();
								P = calculatePercentageDifference();

								// High walking
								if (t <= T/4) {
										
								} 
								// Belly walking
								else if (t <= (2*T/4) && t > T/4) {

								}
								// Hunting coordination
								else if (t <= (3*T/4) && t > (2*T/4)) {

								} 
								// Hunting cooperation
								else if (t <= T && t > (3*T/4)) {

								}
						}
				}
		}
}

int main(int argc, char** argv) {
		// Read from command line. 
		if (argc < 2 || argc > 2) {
				cout << "Usage: " << argv[0] << " <Parameter-settings.txt>" << endl;
				return 0;
		}
		initializeParameters(argv[1]);
		initializeSolution();
		printSolution();
		findBest();
}