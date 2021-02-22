# Cellular Pandemic Model
This repository contains a cellular pandemic model developed with [Cadmium Cell-DEVS](https://github.com/SimulationEverywhere/cadmium). The model enables the definition of multiple pandemic scenarios by tuning a set of configuration parameters.

## Software dependencies
This software relies on the following third-party tools:

1. [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
2. A C++17 compliant version of the `g++` compiler.
3. [Boost C++ libraries](https://www.boost.org).
4. [Jupyter](https://jupyter.org/install).
5. [Pandas](https://pandas.pydata.org).
6. [matplotlib](https://matplotlib.org/stable/index.html).


## Compilation
1. Clone this repository and move inside it:
  ```bash
  git clone https://github.com/SimulationEverywhere-Models/CiSE-Pandemic.git
  cd CiSE-Pandemic
  ```
2. Run the `setup.sh`script:
  ```bash
  ./setup.sh
  ```

## Usage
1. Define a scenario configuration. An example is shown in `config/scenario.json`. All the configuration parameters are described in the [Parameters README file](TODO).
2. Move to the repository's `bin/` folder and run the simulation passing the path to the scenario file and, optionally, the maximum number of time steps before the simulation stops (default is 500):
   ```bash
   cd bin
   ./CiSE-Pandemic ../config/scenario.json [MAX_TIME_STEPS]
   ```

## Visualization
We provide a Jupyter Notebook to visualize the simulation outcome. Instructions for running this notebook can be found in the [Visualization README file](TODO).
