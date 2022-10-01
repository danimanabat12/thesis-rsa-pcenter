#! /bin/bash
cd src
g++ -std=c++11 -o rsa main_rsa.cpp
if [ $? -eq 0 ]; then
  ./rsa PS1.txt
else
  echo "Compilation error. Check your code."
fi