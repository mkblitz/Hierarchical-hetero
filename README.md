# Hierarchical-hetero by Toviah Moldwin, Menachem Kalmenson, and Idan Segev, 2022.

All code to run the experiments found in "Asymmetric voltage attenuation in dendrites can enable hierarchical heterosynaptic plasticity" can be found here.

For help downloading and running NEURON, see Oren Amsalem's tutorial here: https://github.com/orena1/NEURON_tutorial

With NEURON installed, navigate to the code folder and run 'nrnivmodl mods'.

To run an experiment, start by downloading the folder titled 'experiment_scripts', and choose the experiment you want to run. 
Name of script will be the name of the figure it is presented in.

Example: to run the experiment that creates the traces seen in Figure 3A2-A5, open script 'fig_3a_2to5', choose the parameters you want, and run script. 

It is suggested to download all the scripts in the main file, as they are used in the running of the experiments.

## **_Script descriptions:_**

**cell_parent:** Script containing the class _cell_. Object functions include methods to add synapses, measure transfer resistances, etc.

**ball_and_stick:** Parent class for ball and stick model neuron (soma and single dendrite). Uses cell_parent as its parent class. You can see here a list of morphological and biophysical parameters relevant to this model, including nseg for number of segments on the stick (dendrite), lengths, diameters, etc.

**ball_and_many_sticks:** Like ball_and_stick, but for the many dendrite model.

**L5PC:** Parent class for the morphologically detailed layer 5 pyramidal neuron. Morphology comes from folder named _morphologies_.

**model_run:** Script we used to run many of our experiments. You do not need to use this script, unless you want to change additional parameters not accessible by the _experiment_scripts_.

**plot_functions:** Script containing most of the plot functions used to create the figures. 

**utils:** Script containing general helper functions.

**calcium_plasticity_model_run, colored_branches, and peak_ca_plus_plasticity_for_three_locations:** Scripts used for figures 3,4 and parts of 5 and 6. You do not need 
to use this script, unless you want to change additional parameters not accessible by the _experiment_scripts_.

**fig_xx:** Scripts used to run actual experiments and create the figures. These are all found in the folder _experiment_scripts_, and give the option to change several - but not all - the parameters. To change additional parameters you must open scripts mentioned above. xx just a general placeholder for actual figure numbers, e.g. fig_2b_2c, fig_3a_2to5, etc.
