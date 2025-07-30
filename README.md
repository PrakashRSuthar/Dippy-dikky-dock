Create the core environment:
conda env create -f environment.yml

Activate the environment:
conda activate dippydock_env

Install the Python application libraries:
pip install -r requirements.txt

Install Vina:
pip install vina


conda create -n vina_env python=3.9  # Create a new environment (optional but recommended)
conda activate vina_env
conda install -c conda-forge vina