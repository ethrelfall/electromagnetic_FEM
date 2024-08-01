# electrostatic_FEEC.py
# Firedrake attempt to do 2D electrostatics with dielectric consts
# building up to aniso diffusion using FEEC
# dielectric version - harmonic removal not needed
# actually seems to work very well - handles discontinuities in dielectric const
# dielectric problem is analogous to diffusion problem, dielectric const analogous to diffusivity
# see Maxwell's eqs for electrostatics, div D=0, D=eps E = eps (-grad phi)

from firedrake import *

#import netgen  # needs installing in order to try out higher order geometry ...

from firedrake.petsc import PETSc
try:
	from slepc4py import SLEPc
except ImportError:
		import sys
		warning("Unable to import SLEPc, eigenvalue computation not possible (try firedrake-update --slepc)")
		sys.exit(0)


#mesh = Mesh("square_square.msh")
mesh = Mesh("square_circ.msh")
#mesh = Mesh("square_circhole.msh")

circle_eps = 10.0  # dielectric const in circle

x,y = SpatialCoordinate(mesh)

# there are two "good" choices of vector finite element space:

# 1) Raviart-Thomas edge elements - order should be same as that of space SS below
SV=FunctionSpace(mesh,"RTE",1)  # note element type N1curl is equivalent to RTE
#SV=VectorFunctionSpace(mesh,"DG",1)  # does not give converged solution ... reason = no Hilbert complex? 

# 2) Brezzi-Douglas-Marini edge elements - order should be one less than that of space SS below
#SV=FunctionSpace(mesh, "BDME",1)  # note element type N2curl is equivalent to BDME

# ... there is also the "obvious" choice of CG function space, which does not work well:
#SV=VectorFunctionSpace(mesh, "CG", 1)

SS=FunctionSpace(mesh,"CG",1)  # order should be 1, 2, 3, ...
V = SV*SS

# mixed formulation with sigma = div(eps u), u is E (=-grad phi), sigma is div D (should vanish identically)
u, sigma = TrialFunctions(V)
v, tau   = TestFunctions(V)

# spatial dependence of dielectric constant
SS2 = FunctionSpace(mesh,"DG",0)
eps = Function(SS2)
#eps_rule = conditional(And(And(le(x,0.25), ge(x,-0.25)), And(le(y,0.25), ge(y,-0.25))),2.0,1.0)  # square
eps_rule = conditional(le(x*x+y*y,0.0625), circle_eps, 1.0)  # circle
eps.interpolate(eps_rule)

a=( inner(sigma, tau) - inner(u*eps, grad(tau)) +inner(grad((sigma-dot(grad(eps),u))/eps),v)+inner(curl(u),curl(v))) *dx

x, y = SpatialCoordinate(mesh)

# here is attempt at dielectric circle solution BC on normal cpt of electric field
alpha = 0.0625*(1.0-circle_eps)/(1.0+circle_eps)
DielectricFlux = Function(SV)
DielectricFlux.interpolate(as_vector([(1+alpha*(y**2-x**2)/((x**2+y**2)**2)), -2*alpha*x*y/((x**2+y**2)**2)]))

L = (DielectricFlux[0]*tau)*ds(13) - (DielectricFlux[0]*tau)*ds(15) \
  + (DielectricFlux[1]*tau)*ds(16) - (DielectricFlux[1]*tau)*ds(14) \

g=Function(V)

solve(a == L, g)

u, sigma = g.split()
u.rename("E")
sigma.rename("divD")

File("electrostatic_FEEC.pvd").write(u, sigma)

