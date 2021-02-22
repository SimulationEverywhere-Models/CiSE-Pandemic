# Visualization
We provide a Jupyter Notebook to visualize the simulation outcome.
The jupyter notebook script generates different graphs of the susceptible, infected, recovered, and deceased ratio of the population over time.

In order to use the jupyter notebook file:
1. Assert that the `simulation_results/` folder contains the file `output_messages.txt` file (this is the outcome of the simulation).
2. Open a terminal in this directory and run the `launch_visualizer.sh` script:
```bash
./launch_visualizer.sh
```
3. A web browser tab will be automatically opened with the Jupyter Notebooks server. Open the `visualizer.ipynb` notebook.

Feel free to read through the notebook and play around with the parameters and plots.
All the plots shown in the Jupyter notebook will be stored in the `visualization/res` folder.
