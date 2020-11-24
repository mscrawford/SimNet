# PPA

Implementation of the Perfect Plasticity Approximation model (Strigul et al. 2008) for tropical forest succession. Adapted from Rüger et al. 2020 (code written by Caroline Farrior, cfarrior@gmail.com; https://github.com/cfarrior/Ruger_etal_2020). Parameterized for the SimNet BEF experimental protocol, to test how local coexistence mechanisms interact with regional metacommunity processes to drive a plant community's emergent BEF patterns.

Can be run serially, or with `snow` in parallel on a local or cluster computer. On a local computer, it should run out-of-the-box in RStudio. The cluster submission script will have to be tailored to the given HPC.

## References 

Rüger, N. et al. 2020. Demographic trade-offs predict tropical forest dynamics. - Science (80-. ). 368: 165–168.

Strigul, N. et al. 2008. Scaling from Trees to Forests: Tractable Macroscopic Equations for Forest Dynamics. - Ecol. Monogr. 78: 523–545.

