# Hierarchical-hetero by Toviah Moldwin, Menachem Kalmenson, and Idan Segev, 2022.

All code to run the experiments found in "Asymmetric voltage attenuation in dendrites can enable hierarchical heterosynaptic plasticity" can be found here.

For help downloading and running NEURON, see Oran Amsalem's tutorial here: https://github.com/orena1/NEURON_tutorial

To run an experiment, start by downloading the folder titled 'experiment_scripts', and choose the experiment you want to run. 
Name of script will be the name of the figure it is presented in.

Example: to run the experiment that creates the traces seen in Figure 3A2-A5, open script 'fig_3a_2to5', choose the parameters you want, and run script. 

It is suggested to download all the scripts in the main file, as they are used in the running of the experiments.

## **_Script descriptions:_**

**cell_parent:** Script containing the class _cell_. Object functions include methods to add synapses, measure transfer resistances, etc.

**ball_and_stick:** Parent class for ball and stick model neuron (soma and single dendrite). Uses cell_parent as its parent class. You can see here a list of morphological and biophysical parameters relevant to this model, including nseg for number of segments on the stick (dendrite), lengths, diameters, etc.

**ball_and_many_sticks:** Like ball_and_stick, but for the many dendrite model.

**L5PC:** Parent class for the morphologically detailed layer 5 pyramidal neuron. Morphology comes from folder named _morphologies_.
