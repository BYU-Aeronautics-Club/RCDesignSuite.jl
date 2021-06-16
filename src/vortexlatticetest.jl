using VortexLattice

#MAIN WING
xle = [0.0; 0.0]
yle = [0.0; 750.0]
zle = [0.0; 0.0]

chord = [300.0; 300.0]
twist = [-3.0*pi/180; 0.0]
dihedral = [0.0; 0.0]

ns = 12
nc = 4

spacing_s = Sine()
spacing_c = Uniform()
mirror = true


# horizontal stabilizer
xle_h = [0.0; 0.0]
yle_h = [0.0; 250.0]
zle_h = [0.0; 0.0]
chord_h = [200.0; 200.0]
theta_h = [0.0; 0.0]
phi_h = [0.0; 0.0]
ns_h = 6
nc_h = 3
spacing_s_h = Uniform()
spacing_c_h = Uniform()
mirror_h = true



# vertical stabilizer
xle_v = [0.0; 0.0]
yle_v = [0.0; 0.0]
zle_v = [0.0; 300.0]
chord_v = [200.0; 150.0]
theta_v = [0.0; 0.0]
phi_v = [0.0; 0.0]
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false


# REFERENCE
cref = 300.0
bref = 1500.0
Sref = cref*bref
cgpose = [135.916; 0.0; -75.863]

ref = Reference(Sref, cref, bref, cgpose)


# FREESTREAM
alpha = 0.0
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(alpha, beta, Omega)

symmetric = false


# HORSESHOE SETUP
# horseshoe vortices
wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, twist, dihedral, ns, nc;
    mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

htail = wing_to_horseshoe_vortices(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(htail, [1000.0, 0.0, 0.0])

vtail = wing_to_horseshoe_vortices(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
translate!(vtail, [1000.0, 0.0, 0.0])

surfaces = [wing, htail, vtail]
surface_id = [1, 2, 3]


# SOLVE
AIC = influence_coefficients(surfaces, surface_id, symmetric)
b = normal_velocity(surfaces, ref, fs)
Gamma = circulation(AIC, b)
CF, CM, panelprops = near_field_forces(surfaces, surface_id, ref, fs, symmetric, Gamma; frame=Stability())
CDiff = far_field_drag(surfaces, ref, fs, symmetric, Gamma)


# FORCE/MOMENT COEFFS
CD, CY, CL = CF
Cl, Cm, Cn = CM


# DERIVATIVES
dCFb, dCMb = body_derivatives(surfaces, ref, fs, symmetric, AIC)

CXu, CYu, CZu = dCFb.u
CXv, CYv, CZv = dCFb.v
CXw, CYw, CZw = dCFb.w
CXp, CYp, CZp = dCFb.p
CXq, CYq, CZq = dCFb.q
CXr, CYr, CZr = dCFb.r

Clu, Cmu, Cnu = dCMb.u
Clv, Cmv, Cnv = dCMb.v
Clw, Cmw, Cnw = dCMb.w
Clp_b, Cmp_b, Cnp_b = dCMb.p
Clq_b, Cmq_b, Cnq_b = dCMb.q
Clr_b, Cmr_b, Cnr_b = dCMb.r

dCFs, dCMs = stability_derivatives(surfaces, ref, fs, symmetric, AIC)

CDa, CYa, CLa = dCFs.alpha
CDb, CYb, CLb = dCFs.beta
CDp, CYp, CLp = dCFs.p
CDq, CYq, CLq = dCFs.q
CDr, CYr, CLr = dCFs.r

Cla, Cma, Cna = dCMs.alpha
Clb, Cmb, Cnb = dCMs.beta
Clp_s, Cmp_s, Cnp_s = dCMs.p
Clq_s, Cmq_s, Cnq_s = dCMs.q
Clr_s, Cmr_s, Cnr_s = dCMs.r


# VIZUALIZE
write_vtk("wing", wing)
write_vtk("htail", htail)
write_vtk("vtail", vtail)