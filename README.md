# DREAM CHALLENGE: FINRISK

## How to participate:

We must register to the challenge and create a team with all of us (info [here!](https://www.synapse.org/#!Synapse:syn27130803/wiki/619282))

## Data
To run the code, download ```train/``` and ```test/``` folders into project root. 

* From the synapse if you are logged in, [here](https://www.synapse.org/#!Synapse:syn38067910)
* In Arthemysa are the ```train/``` and ```test/```

## Project workflow
All code must be included in ```src/``` folder. In this folder, we have an ```utils/``` repository to stored those function and/or scripts for general stuffs (such as data preprocessing, etc.).

Each new model/pipeline must be have a folder on ```src/``` folder (see [baseline]() as example). The new models/pipelines have to import data from ```train/``` and predict using ```test/``` data. The predicted results must be saved on ```output/``` as ```.csv``` (see [here](https://www.synapse.org/#!Synapse:syn27130803/wiki/619273) for more information)