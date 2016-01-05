# Introduction #

The prediction of the behavior of mixtures in liquid or vapor phase, including the phase equilibrium, is interesting to Chemical Engineering. Usually this task is executed with models that rely heavily on experimental data. Recently, some models that are (almost) totally predictive became available. In these models the necessary information about the substances come from quantum mechanics calculations.

JCOSMO was initially based on the FORTRAN COSMO-SAC code from [Virginia Tech](http://www.design.che.vt.edu). Currently it contains improvements in code organization, in the method itself and in the global parameters (see the Citation section below).
Regarding the COSMO-SAC method, it is actually a variation of the COSMO-RS method (J. Phys. Chem., 1995, 99, 2224-2235).

The purpose of this study is to develop a computational tool capable of predicting properties of mixtures using these models.
Please note that results of COSMO-based methods are frequently semi-quantitative, for more precise results (but less general method) check the [F-SAC method](http://www.enq.ufrgs.br/labs/lvpp).

# Process Simulation #

As far as we know JCosmo is available in two process simulators:
  * [DWSIM](http://dwsim.inforside.com.br): an open-source sequential modular simulator
  * [iiSE](http://www.vrtech.com.br/pt_br/simulador-ise/simulador-de-processo-qu-micos-e-petroqu-micos-ise.html): a new equation based process simulator

# Features #

  * Database with 1432 (sigma profile database with 1432 developed by [Virginia Tech](http://www.design.che.vt.edu/VT-Databases.html) under the supervision of Prof. Y. A. Liu).
  * Can predict activity coefficients for liquid mixtures with any number of components using the COSMO-SAC model
  * Can predict the infinite dilution activity coefficients for binary mixtures
  * A graphical interface for plotting the activity coefficient profiles
  * A graphical interface for plotting the sigma profiles

Activity coefficient profiles of binary mixtures:

![http://jcosmo.googlecode.com/svn/wiki/gamma.png](http://jcosmo.googlecode.com/svn/wiki/gamma.png)

Sigma profiles for two compounds:

![http://jcosmo.googlecode.com/svn/wiki/sigma.png](http://jcosmo.googlecode.com/svn/wiki/sigma.png)

**NOTE**: Currently JCosmo **do not** calculate the sigma profiles, it uses the already prepared sigma profiles from the [VT-2005](http://www.design.che.vt.edu/VT-Databases.html) database. We cannot assure the quality of the sigma profiles used or if they were produced with the correct conformations.

# Documentation #

## Installation Notes ##

> First, download the distribution file on the [downloads tab](http://code.google.com/p/jcosmo/downloads/list).

  * On Windows: just extract the contents of the distribution package and execute jcosmo.exe to run the demonstration dialog.
  * On Linux: just extract the contents of the distribution package and execute the script ./jcosmo.sh to run the demonstration dialog.

**Note:** JCOSMO requires a **Java version 6** or higher.

## Tutorials ##

For some simple tutorials in how to use JCosmo to produce phase equilibrium diagrams and more check the [wiki pages](http://code.google.com/p/jcosmo/w/list).

# Authors #

  * Rafael de Pelegrini Soares
  * Renan Pereira Gerber

# License #

JCosmo is licensed under [LGPL](http://www.gnu.org/licenses/lgpl.html).

# Citation #

If using this code, please cite as:
  * Ind. Eng. Chem. Res., 2010, 49 (16), pp 7488–7496 DOI: [10.1021/ie901947m](http://pubs.acs.org/doi/abs/10.1021/ie901947m)
  * Ind. Eng. Chem. Res., 2011, 50 (5), pp 3060–3063 DOI: [10.1021/ie102087p](http://pubs.acs.org/doi/abs/10.1021/ie102087p)
  * Braz. J. Chem. Eng. vol.30 no.1 São Paulo Jan./Mar. 2013 DOI: [S0104-66322013000100002](http://dx.doi.org/10.1590/S0104-66322013000100002)