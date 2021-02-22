# Configuration the simulation scenario
You don't need to re-compile the model every time you want to explore a new scenario. Instead, you can simply change the input configuration file.
An example configuration file can be found in the `scenario.json` file. In this `README` file, we describe all the configuration parameters.

## `"shape"` (array of integer numbers)
This parameter sets the dimensions of scenario. The shape of the scenario is an array of integer numbers.
In the example configuration file, `"shape"` is set to `[25, 25]`, which means that the scenario is a 2-dimensional grid of 25 by 25 cells.

##`"wrapped"` (boolean)
If set to `true`, the cell grid will be wrapped (i.e., cells in the limit of the scenario will be connected to the opposite limit cells, resulting into an spherical scenario).
By default, `"wrapped"` is set to `false`.

##`"cells"` (JSON object)
This parameter contains all the possible initial cell configurations.
By default, all the cells of the scenario will be initialized with the `"default"` configuration.
However, you can specify which cells start with an alternative configuration in the `"cell_map"` JSON object (more information below).

### `"default"` (JSON object): Default initial cell configuration
The `"default"` JSON object describes the default initial configuration of every cell in the scenario.
You can define the following parameters inside the `"default"` object:
- `"delay"` (string): delay buffer used by the cell. You can choose between `"inertial"`, `"transport"`, and `"hybrid"`.
  For the sake of this model, we recommend you to stick to the `"inertial"` delay, which is the default value.
- `"state"` (JSON object): this object defines the default initial state of the cells. A cell's state contains the following fields:
  - `"population"` (array of integer numbers): this array sets the number of people that lives inside a cell.
  You can divide the population into age segments by adding elements to the `"population"` array.
    For instance, if `"population"` is set to `[100]`, by default all cells will have 100 individuals.
    On the other hand, if `"population"` is set to `[20, 50, 50]`, by default all cells will have 120 individuals divided into 3 age segments (e.g., young, middle-age, and elder).
    - **IMPORTANT**: all cell state initial configurations must have the exact same number of age segments.
  - `"susceptible"` (array of decimal numbers): the ratio of the population of the cell that is susceptible to get infected at the beginning of the simulation. 
    Each number of the array corresponds to a different age group.
  - `"infected"` (array of decimal numbers): the ratio of the population of the cell that is infected at the beginning of the simulation.
    Each number of the array corresponds to a different age group.
  - `"recovered"` (array of decimal numbers): the ratio of the population of the cell that has recovered from the disease (and thus are immune) at the beginning of the simulation.
    Each number of the array corresponds to a different age group.
  - `"deceased"` (array of decimal numbers): the ratio of the population of the cell that has died due to the disease at the beginning of the simulation.
    Each number of the array corresponds to a different age group.
  - **IMPORTANT**: the number of elements of the `"population"`, `"susceptible"`, `"infected"`, `"recovered"`, and `"deceased"` arrays must be the same.
  - **IMPORTANT**: for each age segment, the initial susceptible, infected, recovered, and deceased ratios must be between 0 and 1 (including both limits).
  - **IMPORTANT**: for each age segment, the sum of the susceptible, infected, recovered, and deceased ratios must be equal to 1.
- `"config"` (JSON object): this object defines additional configuration parameters of the model. 
  You can define the following configuration parameters:
  - `"n_decimals"` (unsigned integer number): number of decimal digits. By default, it is set to 2.
  - `"susceptibility"` (array of decimals): for each age segment, it indicates how easily a susceptible person can be infected by the disease (from 0 to 1).
  - `"virulence"` (array of decimals): for each age segment, it indicates how easily an infected person can infect another person (from 0 to 1).
  - `"recovery"` (array of decimals):for each age segment, it indicates how easily an infected person can get recovered (from 0 to 1).
  - `"mortality"` (array of decimals): for each age segment, it indicates how easily an infected person can die due to the disease (from 0 to 1).
  - `"immunity"` (array of decimals): for each age segment, it indicates the probability a recovered person can maintain their immunity (from 0 to 1).
  - **IMPORTANT**: the number of elements of the `"susceptibility"`, `"virulence"`, `"recovery"`, `"mortality"`, and `"immunity"` arrays must be equal to the number of age segments in the scenario.
- `"neighborhood"` (array of JSON objects): this object defines the neighborhood (i.e., which other cells affect a cell's state).
  Each JSON object of the array must define the following parameters:
  - `"vicinity"` (JSON object): it defines the degree of dependency between neighboring cells and the cell.
    In this model, the vicinity is defined by the following fields:
    - `"connectivity"` (decimal number): it defines the degree of connectivity (e.g., transport means) between the cell and its neighboring cells.
    - `"mobility"` (array of decimal numbers): for each age segment, it defines the ratio of people that move from the neighboring cells to the cell as a daily basis.
      - **IMPORTANT**: the number of elements of the `"mobility"` array must be equal to the number of age segments in the scenario.
  - `"type"` (string): neighborhood type. Currently, Cadmium supports [`"moore"`](https://en.wikipedia.org/wiki/Moore_neighborhood), [`"von_neumann"`](https://en.wikipedia.org/wiki/Von_Neumann_neighborhood), and `"relative"` neighborhood types.
  - `"range"` (integer number, only for `"moore"` and `"von_neumann"` neighborhood types): it determines the range of the neighborhood.
    By default, it is set to 1 (i.e., regular neighborhood). Use numbers greater than 1 to use extended neighborhoods.
  - `"neighbors"` (array of arrays of integer numbers, only for `"relative"` neighborhood types):
    array containing the relative location of the cell that comprises the neighborhood.

# Other possible initial cell configurations (JSON object)
Usually, we define more than one initial configuration. To do so, you just have to add a new JSON object with the different configuration parameters.
Note that these additional configurations work as patches of the default configuration: if a parameter is not defined, then it will inherit the default value for that parameter.
In the example scenario, we define the `"infection_epicenter"` initial cell configuration.
We only define the new `"susceptible"` and `"infected"` ratios, as we want to keep the rest of the parameters as in the default configuration.

##`"cell_map"` (JSON object)
This section of the configuration file describes which cells are not initialized with the default configuration.
In the example configuration, we select the cell `[12, 12]` (i.e., the one in the middle of the scenario) to start with the `"infection_epicenter"` initial configuration.
