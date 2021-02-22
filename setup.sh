#!/bin/bash

# Download and update all Git submodules
git submodule update --init --recursive

# Create required directories
mkdir simulation_results
mkdir bin

# Compile pandemic model
cd bin || { echo "Failed to change to bin directory"; exit; }
g++ -g -I../cadmium/include -I../json/include -std=c++17 -o CiSE-Pandemic ../model/main.cpp
