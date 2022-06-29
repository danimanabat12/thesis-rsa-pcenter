#include<iostream>
#include<random>
#include<vector>
#include<time.h>
#include<fstream>
#include<sstream>
#include<string>
#include<typeinfo>

using namespace std; 

double demandPoints[6][2] = {{10, 30}, {-10, 30}, {10, 30}, {5, 15}, {-5, 15}, {5, -15}};


int main(int argc, char** argv) {

  // string xd = "835671.5378";
  // // double xd2 = 835671.5378;
  // printf("%f", atof(xd.c_str()));

  vector<vector<double> > demandPointsVec;
  vector<double> temp;
  fstream dataset; 
  string line, word; 
  dataset.open(argv[1], ios::in);
  while (!dataset.eof()) {
    temp.clear();
    getline(dataset, line); 
    stringstream s(line);
    while(getline(s, word, ',')) {
      temp.push_back(atof(word.c_str()));
    }
    cout << "(";
    for (int i = 0; i < temp.size(); i++) {
      printf("%f ", temp.at(i));
    }
    cout << ")" << endl;
    demandPointsVec.push_back(temp);
  }
  dataset.close();


  // for (int i = 0; i < 6; i++) {
  //   for (int j = 0; j < 2; j++) {
  //     cout << demandPoints[i][j] << " "; 
  //   }
  //   cout << endl;
  // }
  // cout << "=====" << endl;
  // for (int i = 0; i < demandPointsVec.size(); i++) {
  //   for (int j = 0; j < demandPointsVec[i].size(); j++) {
  //     cout << demandPointsVec[i][j] << " "; 
  //   }
  //   cout << endl;
  // }
  // srand(time(NULL));
  // cout << (double) rand()/RAND_MAX << endl;

  // if (params.is_open()) {
  //   while(getline(params, line)) {
  //     cout << line << endl;
  //   }
  //   params.close();
  // }
  // int repNum = 30; 
  // int fac = 6; 

  // vector < vector< vector<double> > > candidateSolution;
  // candidateSolution.resize(30); 
  // for (int i = 0; i < candidateSolution.size(); i++) {
  //   cout << i+1 << endl;
  //   candidateSolution[i].resize(fac);
  //   for (int j = 0; j < candidateSolution[i].size(); j++) {
  //     // cout << ". ";
  //     cout << candidateSolution[i][j].size() << " ";
  //     candidateSolution[i][j].push_back(13);
  //     candidateSolution[i][j].push_back(15);
  //   }
  //   cout << candidateSolution[i][1].size();
  //   cout << endl;
  // }
}