Palantir
------

Palantir is an algorithm to align cells along differentiation trajectories. Palantir models differentiation as a stochastic process where stem cells differentiate to terminally differentiated cells by a series of steps through a low dimensional phenotypic manifold. Palantir effectively captures the continuity in cell states and the stochasticity in cell fate determination. Palantir has been designed to work with multidimensional single cell data from diverse technologies such as Mass cytometry and single cell RNA-seq. 


#### Installation and dependencies
1. Palantir has been implemented in Python3 and can be installed using:

        $> git clone git://github.com/dpeerlab/Palantir.git
        $> cd Palantir
        $> sudo -H pip3 install .

2. Palantir depends on a number of `python3` packages available on pypi and these dependencies are listed in `setup.py`
All the dependencies will be automatically installed using the above commands

3. To uninstall:
		
		$> sudo -H pip3 uninstall palantir

4. If you would like to determine gene expression trends, please install <a href="https://cran.r-project.org"> R <a> programming language and the R package <a href="https://cran.r-project.org/web/packages/gam/">GAM </a>. You will also need to install the rpy2 module using 
	
		$> sudo -H pip3 install rpy2
		
5. Palantir can also be used with [**Scanpy**](https://github.com/theislab/scanpy). It is fully integrated into Scanpy, and can be found under Scanpy's external modules ([link](https://scanpy.readthedocs.io/en/latest/api/scanpy.external.html#external-api))


#### Usage

A tutorial on Palantir usage and results visualization for single cell RNA-seq data can be found in this notebook: http://nbviewer.jupyter.org/github/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb


#### Citations
Palantir manuscript is in submission: "Palantir characterizes cell fate continuities in human hematopoiesis". Preprint is available at https://www.biorxiv.org/content/early/2018/08/05/385328
