# Instructions

To use this model, you will need to pull the container using singularity. First, make sure:
* Your system is compatible (Linux OS)
* You have singularity installed (https://docs.sylabs.io/guides/latest/user-guide/)

## Opening the container
1. Pull the container to your machine using `singularity pull docker://uvarc/rhessys:leb3t`
2. Run the container using singularity run `docker://uvarc/rhessys:leb3t`

## Start Jupyter Notebook session
On local machine: `jupyter notebook`
 
Or, if accessing system via ssh:
  1. `jupyter notebook --no-browser --port=8080`
  2. `ssh -L 8080:localhost:8080 <REMOTE_USER>@<REMOTE_HOST>`
  3. Navigate to `http://localhost:8080/`

