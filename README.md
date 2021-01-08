# Exit_time
MATLAB codes to support Simpson et al. (2021)  Exit time for diffusion on irregular domains

These MATLAB codes require GMSH 4.7.1 to generate unstructured meshes for the finite volume and random walk calculations.  The software can be downloaded at https://gmsh.info/

The MATLAB codes also include ellipse fitting tools from Szpak et al. (2015) https://link.springer.com/article/10.1007%2Fs10851-014-0536-x

- The main MATLAB scripts are
  - *unperturbed_disc_main.m* for the unperturbed disc problem;
  -  *unperturbed_ellipse_main.m* for the unperturbed ellipse problem;
  - *perturbed_disc_main.m* for the perturbed disc problem; 
  - *perturbed_ellipse_main.m* for the perturbed ellipse problem;
  - *tasmania_analysis.m* for the Tasmania case study;
  - *taiwan_main.m* for the Taiwan case study.

- Some notes for each figure:
  - To produce **Figure 1**, run *unperturbed_disc_main.m* Lines 1 through 64.
  - To produce **Figure 2**, run *unperturbed_ellipse_main.m* Lines 1 through 64.
  - To produce **Figure 3**, run *perturbed_disc_main.m* Lines 1 through 78,
  - To produce **Figure 4**, run *perturbed_ellipse_main.m* Lines 1 through 82.
  - To produce **Figure 6**, run *tasmania_analysis.m* Lines 1 through 214. 
  - To produce **Figure 7**, run *taiwan_analysis.m* Lines 1 through 192. 
  - To produce **Figure 8**, run *perturbed_disc_main.m* Lines 1 through 29, Line 37, Lines 80 through 136.
  - To produce **Figure 9**, run *perturbed_ellipse_main.m* Lines 1 through 42, Lines 84 through 123.

- About inpoly2:
- In *perturbed_ellipse_walk_2.m*, the inpoly2 function from Engwirda (2020) [https://github.com/dengwirda/inpoly] with a slight modification to prevent errors in the case of a single remaining walk. This modification is given in the documentation of *perturbed_ellipse_walk_2.m*.


