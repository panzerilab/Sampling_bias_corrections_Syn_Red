# **README**

## **Contents**
- `main_code_bias_corrs_on_data.py` – Main function to correct the bias of Partial Information Decomposition (PID) measures in Gaussian data.
- `gpid/` – Package developed by Venkatesh et al.
- `utils_parallelization.py` – Utilities for parallelizing estimation methods.
- `Gaussian_bias_correction_methods_routines.py` – Collection of bias correction methods including 1) 'Venkatesh_bias_correction': the Venkatesh bias correction approach (Venkatesh et al., 2023 Nips), 2) 'informative_bias_correction': our best bias correction approach, 3) 'shuffsub_bias_correction', 4) 'zero_info_bias_correction'. The 'no_bias_correction' outputs the biased (naive) values of information without making any bias correction. 
- `environment.yml` – Dependencies for Conda environment.
- `requirements.txt` – Dependencies for virtual environments (venv).
- `Dockerfile` -  Defines a Docker image for running the project.

---

## **Setup**

### **Option 1: Conda (Recommended)**
```bash
conda env create -f environment.yml
conda activate gpid
```

### **Option 2: venv**
```bash
python -m venv myenv
source myenv/bin/activate  # macOS/Linux
myenv\Scripts\activate     # Windows
pip install -r requirements.txt
```
Deactivate with `deactivate`.

### **Option 3: Docker**
```bash
docker build -t my_python_env .
docker run -it -v "$(pwd)":/app --name my_container my_python_env /bin/bash
cd /app
```

---

## **Installing `gpid`**
```bash
cd gpid
pip install .
```

---

## **Usage: `execute_bias_correction_routines`**
This is the main function in the folder, used to apply bias correction methods to measured PID values in datasets.


### **Function Signature**
```python
execute_bias_correction_routines(input_data, bias_routines, num_iterations=10, parallel=False, save=False, save_format="matlab")
```

### **Parameters**
- **`input_data` (np.ndarray)** – Dataset in shape `(samples, dimensions)`, with variables ordered as `[M, X, Y]`:
  - **M** – Target variable  
  - **X, Y** – Source variables  

- **`bias_routines` (list of str)** – Available correction methods:
  - `"no_bias_correction"`
  - `"informative_bias_correction"`
  - `"zero_information_bias_correction"`
  - `"shuffle_subtraction_bias_correction"`
  - `"uniform_bias_correction"`

- **`num_iterations`** (`int`, default: `10`) – Number of randomization used in the bias correction methods.  
- **`parallel`** (`bool`, default: `False`) – Run in parallel (`True`) or sequentially (`False`).  
- **`save`** (`bool`, default: `False`) – Whether to save results.  
- **`save_format`** (`str`, default: `"matlab"`) – Save format:
  - `"matlab"` → `.mat`
  - `"python"` → `.pkl`
  - `"json"` → `.json`

### **Return Value**
A list of tuples `(bias_title, unbiased_results, error)`, where:
- `bias_title` – Name of applied bias correction.
- `unbiased_results` – Processed dataset.
- `error` – Any encountered error.

---

## **Example Usage**

### **1. Run in Parallel & Save as `.mat`**
```python
execute_bias_correction_routines(
    input_data, 
    ["plugin_values", "shuffle_subtraction_bias_correction"], 
    num_iterations=15, 
    parallel=True, 
    save=True, 
    save_format="MATLAB"
)
```

### **2. Run Sequentially & Save as `.pkl`**
```python
execute_bias_correction_routines(
    input_data, 
    ["informative_bias_correction"], 
    num_iterations=10, 
    parallel=False, 
    save=True, 
    save_format="Python"
)
```

---

## **Output Structure**
Each method returns a dictionary with the following quantities:

### **Information Measures**
- `"imxy"` – Mutual information: M & {X, Y}
- `"uix"` – Unique info: X → M
- `"uiy"` – Unique info: Y → M
- `"imx"` – Mutual info: X → M
- `"imy"` – Mutual info: Y → M
- `"ri"` – Redundant info: X & Y → M
- `"si"` – Synergistic info: X & Y → M

### **Entropy Measures**
- `"hm"` – Entropy of M
- `"hxy"` – Entropy of {X, Y}
- `"hmxy"` – Joint entropy: {M, X, Y}

### **Conditionally Independent Distribution (Pind)**
- `"imxy_ind"`
- `"hm_ind"`
- `"hxy_ind"`
- `"hmxy_ind"`

---

