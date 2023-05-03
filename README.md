# Code for the paper:
> **Proximity measures for orchard networks**  
> *Leo van Iersel, Mark Jones, Esther Julien, Yukihiro Murakami*

To run the code, the following packages have to be installed: `networkx`, `numpy`, and `PySCIPOpt`. 
The used version can be found in `requirements.txt`.

### Generate simulated networks
In the folder `realNetworks` and `simulatedNetworks`, the real and simulated networks that were used in the experiments 
of the paper are given, respectively. 

The simulated networks were generated with the code in `network_gen.py` by running the following: 
```bash
python network_gen.py -l -r -s
```
where `l` is the number of leaves, `r` the number of reticulations and `s` the instance number.

### Leaf distance MILP
The MILP for the leaf distance can be found in `leaf_distance_milp.py`.
To run real networks:
```bash
python leaf_distance_milp.py -real=1 -s
```
where `s` is the instance number.
To run simulated networks, the number of leaves (`l`) and the number of reticulations (`r`) have to be indicated:
```bash
python leaf_distance_milp.py -real=0 -l -r -s
```
In the experiments of the paper, 50 instances were made per pair (l, r) for l in {20, 50, 100, 150, 200} and for 
r in {5, 10, 20, 30, 40, 50, 100, 200}.

- If you have any questions or comments, please send an email to Esther Julien via e.a.t.julien@tudelft.nl.