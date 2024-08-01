# electrostatic_primal.py
# a Firedrake attempt to do 2D electrostatics with piecewise constant dielectrics
# simple primal scalar Laplace equation
# this is dielectric version not hole
# as long as the dielectric constant is DG, this version works fine
# it works because it solves for potential phi which is fully continuous
# (so, FEEC is NOT needed for this problem)

from firedrake import *

mesh = Mesh("square_circ.msh")

x,y = SpatialCoordinate(mesh)

SS=FunctionSpace(mesh,"CG",1)

phi = TrialFunction(SS)
phi_test = TestFunction(SS)

circle_eps = 10.0
SS2 = FunctionSpace(mesh, "DG", 0)
eps = Function(SS2)
eps_rule = conditional(le(x*x+y*y,0.0625),circle_eps,1.0)  # circle
eps.interpolate(eps_rule)

a = (eps*dot(grad(phi),grad(phi_test)))*dx # dielectric

f=Function(SS)
f.interpolate(0.0*x+0.0)
L=inner(f,phi_test) *ds(13)

g=Function(SS)

# analytic perfect conducting metal, or dielectric, boundary condition, with circle
phi_boundary = Function(SS)
#phi_s.interpolate(x+0.0625*x/(x**2+y**2))  # perfect conductor
phi_boundary.interpolate(x+0.0625*((1.0-circle_eps)/(1.0+circle_eps))*x/(x**2+y**2))  # dielectric
bc_cond = DirichletBC(SS, phi_boundary, "on_boundary")

nullspace = VectorSpaceBasis(constant=True)

solve(a == L, g, bcs=[bc_cond], nullspace=nullspace)

g.rename("phi")

VV = VectorFunctionSpace(mesh,"DG",1)

grad_g = Function(VV)
grad_g.project(grad(g))
grad_g.rename("grad_phi_s")

VS = FunctionSpace(mesh,"DG",0)
Laplacian_g = Function(VS)
Laplacian_g.project(div(grad(g)))
Laplacian_g.rename("Laplacian_phi_s")

phi_boundary.rename("phi_exact")

File("electrostatic_primal.pvd").write(g, grad_g, Laplacian_g, phi_boundary)

