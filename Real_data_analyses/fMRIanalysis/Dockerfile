FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the environment.yml file into the container
COPY environment.yml /app/environment.yml

# Create the conda environment
RUN conda env create -f /app/environment.yml

# Activate the environment and set it as default
RUN echo "conda activate $(head -1 /app/environment.yml | cut -d' ' -f2)" >> ~/.bashrc

# Use the conda environment as default in the container
ENV PATH /opt/conda/envs/$(head -1 /app/environment.yml | cut -d' ' -f2)/bin:$PATH

# Copy the rest of the application files
COPY . /app

# Install gpid
WORKDIR /app/gpid
RUN pip install .
WORKDIR /app

# Set the default command
CMD ["/bin/bash"]
