# electromagnetic_FEM
Experiments with various FEM techniques in electromagnetism, starting with electrostatics (note analogy with diffusion problems).

# Electrostatics

The problem is $\nabla \cdot D = 0$, $E=-\nabla \phi$ for a continuous scalar potential function $\phi$.  $D= \epsilon E$.

$D_{\perp}$ and $E_{\parallel}$ are continuous at material boundaries.  

For a primal formulation of the problem of a 2D circular dielectric in a uniform applied electric field see electrostatic_primal.py.

For a more advanced FEEC formulation of the same problem, see electrostatic_FEEC.py.  This uses a choice of elements such that $\sigma=\nabla \cdot D$ is fully continuous (in fact it vanishes identically) and $u=E=-\grad{\phi}$ is continuous only in the component parallel to the element edges (this is clearly good since material boundaries lie along element edges in a conforming discretization, as here).

![E_magnitude_primal](../png/Ed_magnitude_primal.png "Magnitude of the electric field from the primal solution.")

![E_magnitude_FEEC](../png/Ed_magnitude_FEEC.png "Magnitude of the electric field from the FEEC solution as described in the text.")
