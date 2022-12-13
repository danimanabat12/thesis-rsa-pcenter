#! /bin/bash
# We shall remove this in the future and perhaps add a workaround for each parameter setting texts?
g++ -std=c++11 -o rsa main_rsa.cpp
if [ $? -eq 0 ]; then
  ./rsa MPC2.txt
else
  echo "Compilation error. Check your code."
fi