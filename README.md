# Volterra-paper

This repo contains the code used to generate plots presented in the paper [Development of an Efficient Formulation for Volterra’s Equations of Motion for Multibody Dynamical Systems](https://scientiairanica.sharif.edu/article_23551.html). 

Three folders are provided that correspond to three case studies provided in the paper. Each folder contain subfolders named after the method used. Each subfolder contains a `RunMe.m` that should be used to run various approaches. These approaches are: 

1. Gibbs-Appel
2. Maggi
3. Lagrange
4. Efficient Volterra (proposed).

Also, to regenerate the exact plots shown in the paper, run the code inside `Graphs' subfolder in each folder. Run the `RunMe.m` in `Errors Graph' folder to generate the error plots throughout the paper, and run `RunMe.m` in `GCs graph' to plot the graph of generalized coordinates.

If you found this code useful and used it in your research code, please cite the paper as

> @article{yoosefian2024development,
> title={Development of an Efficient Formulation for Volterra’s Equations of Motion for Multibody Dynamical Systems},
> author={Yoosefian Nooshabadi, Mohammad Hussein and Nejat Pishkenari, Hossein},
> journal={Scientia Iranica},
> year={2024},
> publisher={Sharif University of Technology}
> }

# First case study
This example considers a cart and a 2-DOF pendulum, as shown below.
