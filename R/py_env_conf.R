# https://rstudio.github.io/renv/articles/python.html
library(tidyverse)
library(reticulate)

# Instalamos la versión de python que vamos a usar
version <- "3.10.8"
install_python(version)

# Creamos un entorno virtual, con el python que hemos instalado, en
# "~/.virtualenvs"
SYSTEM_VIRTUAL_ENV <- paste0("virt_env_", version)
if (!virtualenv_exists(SYSTEM_VIRTUAL_ENV)) {
  virtualenv_create(SYSTEM_VIRTUAL_ENV, version = version)
}
# use_virtualenv("cofidis_clv23")

# python.exe del entorno virtual creado
python_exe <- file.path(virtualenv_root(), SYSTEM_VIRTUAL_ENV, "Scripts", "python.exe")

# Creamos el entorno virtual del proyecto en 
renv::use_python(name = "tst_lifetimes", python = python_exe)

py_install(c("ipython", "jupyter", "jupyterhub", "jupyterlab", "ipykernel",
             "pandas", "numpy", "statsmodels", "scipy",
             "scikit-learn", "matplotlib", "seaborn", "lifetimes"))
