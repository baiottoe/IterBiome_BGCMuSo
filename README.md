# Instructions

To use this model, you will need to pull the container using singularity. First, make sure:
* Your system is compatible (Linux OS)
* You have singularity installed (https://docs.sylabs.io/guides/latest/user-guide/)

## Opening the container
1. Pull the container to your machine using `singularity pull docker://uvarc/rhessys:leb3t`
2. Run the container using singularity run `docker://uvarc/rhessys:leb3t`

## Cloning this github repository
Now, you will want to install the contents of this github repository, such that it can be accessed within the container. 
1. Change directory to where you would like the contents of this repository to go using `cd <DIRECTORY>`
2. Clone this repository locally with `git clone git@github.com:baiottoe/BBGC_Muso_I.git`

## Start Jupyter Notebook session
On local machine: `jupyter notebook`
 
Or, if accessing system via ssh:
  1. `jupyter notebook --no-browser --port=8080`
  2. `ssh -L 8080:localhost:8080 <REMOTE_USER>@<REMOTE_HOST>`
  3. Navigate to `http://localhost:8080/`

## Running the Iterative BBGC-Muso Model
From the jupyter notebook window you opened, navigate to where you cloned this repository, open up the **Step_1** notebook file and proceed following these instructions. The notebooks are designed to walk you through the entire modelling process.
