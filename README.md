# Lund Jet Plane: 
## Environment Setup

This guide outlines the steps to set up the Monte Carlo (MC) Generators required for this project. We will be installing MadGraph, Pythia 8, and Delphes.

### 1. Workspace Setup
Create a dedicated directory for MC Generators in your workshpace and navigate in it:
```bash
mkdir MCGenerators
cd MCGenerators
```

### 2. Madgraph
MadGraph can be downloaded directly from the official website.
Note: I have used version 3.5.5.

```bash
# Download Pythia 8.312
$ wget [https://www.pythia.org/download/pythia83/pythia8312.tgz].(https://www.pythia.org/download/pythia83/pythia8312.tgz)

# Unzip the tarball
$ tar xvfz pythia8312.tgz
```
### 3. Pythia
Pythia can be downloaded and compiled from the official website.
Note: I have used version 8.312

```bash
# Return to the MCGenerators directory
$ cd ..

# Download Pythia 8.312
$ wget [https://www.pythia.org/download/pythia83/pythia8312.tgz].(https://www.pythia.org/download/pythia83/pythia8312.tgz)

# Unzip the tarball
$ tar xvfz pythia8312.tgz

# Compile Pythia
$ cd pythia8312
$ make
```
### 4. Delphies
```bash
# Return to the MCGenerators directory
$ cd ..

# Download Delphes 3.5.0
$ wget [http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz](http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz)

# Unzip the tarball
$ tar -zxf Delphes-3.5.0.tar.gz

# Compile Delphes
$ cd Delphes-3.5.0
make
```
### 5. Integrating Pythia 8 with Delphes

Important: Do not use Pythia standalone, as it can be quite unstable for this workflow. Instead, we will use it inside Delphes by creating an executable called DelphesPythia8.

** Make the Pythia library available globally:**
Edit your .bashrc file and add the following line:

``` bash
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rastogi/Uttsavir/Work/MCGenerators/pythia8312/lib
# /home/rastogi/Uttsavir/Work/MCGenerators/pythia8312 is the path to your Pythia folder
```
In the Delphes folder, declare a global variable at the terminal and recompile:

``` bash
$ export PYTHIA8=/home/rastogi/Uttsavir/Work/MCGenerators/pythia8312
$ make HAS_PYTHIA8=true
```    
This will compile an executable called DelphesPythia8, which you'll see in your Delphes folder.


## Workflow
### 1. Generating the events
First, produce NLO events in MadGraph:

``` bash
MG5_aMC> generate process_name [QCD]
```

Take the resulting LHE file from MadGraph. Use the `delphes_card_IDEA.tcl` detector card and your Pythia 8 configuration file to shower the events and run the fast detector simulation, generating a ROOT file:

```bash
$ ./DelphesPythia8 cards/delphes_card_IDEA.tcl pythia_configLHE.cmnd output_file.root
```
### 2. Running the code on the generated sample
Go to the directory where the analysis files are located. Compile them from the terminal prompt:
``` bash
$ source COMPILEIT
```
Next, execute the driver file to feed the appropriate values to the input-output parameters and run the analysis:
``` bash
$ ./ana 2
```

