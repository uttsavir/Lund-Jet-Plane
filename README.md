# Lund Jet Plane: 
# Environment Setup

This guide outlines the steps to set up the Monte Carlo (MC) Generators required for this project. We will be installing MadGraph, Pythia 8, and Delphes.

## 1. Workspace Setup
First, create and navigate into a dedicated directory for your MC Generators:
```bash
mkdir -p MCGenerators
cd MCGenerators

## 2. Madgraph

## 3. Pythia

# Download Pythia 8.312
wget [https://www.pythia.org/download/pythia83/pythia8312.tgz](https://www.pythia.org/download/pythia83/pythia8312.tgz)

# Unzip the tarball
tar xvfz pythia8312.tgz

# Compile Pythia
cd pythia8312
make

## 3. Delphies

# Return to the MCGenerators directory
cd ..

# Download Delphes 3.5.0
wget [http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz](http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz)

# Unzip the tarball
tar -zxf Delphes-3.5.0.tar.gz

# Compile Delphes
cd Delphes-3.5.0
make
