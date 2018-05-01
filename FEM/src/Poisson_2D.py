import dolfin
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import os
import subprocess
import mesh_writer
import numpy as np

def worker(amp_in, offset, timestep, step, params, outfile):
    #elec_length, elec_spacing use 10nm as conservative estimate of current generation feature size.
    #elec_height uses 20nm, shown possible for Cr here: https://www.nature.com/articles/srep23823.pdf
    #the voltage at the boundary increases with elec_length and elec_spacing.
    #the voltage at the boundary decreases with boundary_dielectric by the same amount. 
    #E.G for voltage at boundary V(elec_spacing, elec_length, boundary_dielectric), V(1, 1, 1) == V(2, 2, 2) == V(3, 3, 3)...
    #Also roughly linear with applied V.
    elec_length, elec_height, elec_spacing, \
    boundary_x_min, boundary_x_max, boundary_dielectric, \
    boundary_y_min, boundary_y_max, mid_x, mid_y = params

    if elec_height/2.0 >= boundary_dielectric:
        print "Electrode height too tall, boundary conditions will intersect with Air-Si boundary. Leaving..."
        return
    # Create classes for defining parts of the boundaries and the interior
    # of the domain
    class Left(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global boundary_x_min
            return dolfin.near(x[0], boundary_x_min)

    class Right(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global boundary_x_max
            return dolfin.near(x[0], boundary_x_max)

    class Bottom(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global boundary_y_min
            return dolfin.near(x[1], boundary_y_min)

    class Top(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global boundary_y_max
            return dolfin.near(x[1], boundary_y_max)

    # INTERNAL BOUNDARY CONDITION
    # ELECTRODE_Phase_shift
    class Electrode_0(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global elec_length, elec_spacing, mid_x
            return ( \
                (dolfin.between(x[0], (0.5*elec_spacing, 0.5*elec_spacing+elec_length))) and \
                dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height)) 
               # (dolfin.between(x[0], (mid_x-0.5*elec_length+2*(elec_length+elec_spacing), mid_x+1.999*(elec_length+elec_spacing))) and \
               #  dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) or \
               # (dolfin.between(x[0], (mid_x-1.999*(elec_length+elec_spacing), mid_x+0.5*elec_length-2*(elec_length+elec_spacing))) and \
               #  dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
            )

    class Electrode_90(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global elec_length, elec_spacing, mid_x
            return ( \
                (dolfin.between(x[0], (1.5*elec_spacing+elec_length, 1.5*elec_spacing+2*elec_length))) and \
                dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))
               # (dolfin.between(x[0], (mid_x-0.5*elec_length-(elec_length+elec_spacing), mid_x+0.5*elec_length-(elec_length+elec_spacing))) and \
               #  dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
            )

    class Electrode_180(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global elec_length, elec_spacing, mid_x
            return ( \
                (dolfin.between(x[0], (2.5*elec_spacing+2*elec_length, 2.5*elec_spacing+3*elec_length))) and \
                dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))
               # (dolfin.between(x[0], (mid_x-0.5*elec_length, mid_x+0.5*elec_length)) and \
               #  dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
            )

    class Electrode_270(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global elec_length, elec_spacing, mid_x
            return ( \
                (dolfin.between(x[0], (3.5*elec_spacing+3*elec_length, 3.5*elec_spacing+4*elec_length))) and \
                dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))
               # (dolfin.between(x[0], (mid_x-0.5*elec_length+(elec_length+elec_spacing), mid_x+0.5*elec_length+(elec_length+elec_spacing))) and \
               #  dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
            )

    # INTERNAL BOUNDARY CONDITION
    #DIELECTRIC
    class Air(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global boundary_dielectric, boundary_y_max
            return (dolfin.between(x[1], (boundary_dielectric, boundary_y_max)) )

    # Sub domain for Periodic boundary condition
    class PeriodicBoundary(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            # global boundary_x_min
            return bool(x[0] < dolfin.DOLFIN_EPS + boundary_x_min and x[0] > -(dolfin.DOLFIN_EPS + boundary_x_min) and on_boundary)
        def map(self, x, y):
            # global boundary_x_max
            #map the left domain to the right domain to boundary_x_max
            y[0] = x[0] - boundary_x_max 
            y[1] = x[1]

    # Create periodic boundary condition
    pbc = PeriodicBoundary()

    # Initialize sub-domain instances
    print "Initialising subdomain instances..."
    left = Left()
    top = Top()
    right = Right()
    bottom = Bottom()
    electrode_0 = Electrode_0()
    electrode_90 = Electrode_90()
    electrode_180 = Electrode_180()
    electrode_270 = Electrode_270()
    air = Air()

    #grab the mesh from file.
    mesh = dolfin.Mesh('../data/domain.xml')

    # Initialize mesh function for interior domains
    print "Marking domains..."
    domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    domains.set_all(0)
    air.mark(domains, 1)

    # Initialize mesh function for boundary domains
    print "Marking boundaries..."
    boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    boundaries.set_all(0)
    left.mark(boundaries, 1)
    top.mark(boundaries, 2)
    right.mark(boundaries, 3)
    bottom.mark(boundaries, 4)
    electrode_0.mark(boundaries, 5)
    electrode_90.mark(boundaries, 6)
    electrode_180.mark(boundaries, 7)
    electrode_270.mark(boundaries, 8)


    # Define function space and basis functions
    print "Defining function space and basis..."
    V = dolfin.FunctionSpace(mesh, "CG", 3, constrained_domain=pbc)
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    # Define Dirichlet boundary conditions at top, bottom, front, back
    print "Defining Dirichlet boundaries..."

    ############# This entire part is kinda ehhhhhhhhh...
    #Should pick a metal with LOWER WF than chi_si to get ohmic contact
    #Note: at 4K, there are no free carriers inside Si -> no depletion region -> no 
    chi_si = 4.05 #eV
    phi_gold = 5.1 #eV
    phi_bi = phi_gold - chi_si
    x = np.linspace(0, 2*np.pi, 4, endpoint=False)
    
    amp = amp_in
    elec_pot = amp*np.sin(x+timestep) + offset
    print "Potentials used:"
    for pot in elec_pot:
        print pot
    #############
    bcs = [dolfin.DirichletBC(V, elec_pot[0]+phi_bi, boundaries, 5),
           dolfin.DirichletBC(V, elec_pot[1]+phi_bi, boundaries, 6),
           dolfin.DirichletBC(V, elec_pot[2]+phi_bi, boundaries, 7),
           dolfin.DirichletBC(V, elec_pot[3]+phi_bi, boundaries, 8)]
           # dolfin.DirichletBC(V, np.abs(amp/boundary_y_max), boundaries, 2),
           # dolfin.DirichletBC(V, np.abs(amp/boundary_y_min), boundaries, 4),
           # dolfin.DirichletBC(V, 0, boundaries, 2)]
           # dolfin.DirichletBC(V, 0, boundaries, 4)]
           # dolfin.DirichletBC(V, 0, boundaries, 1),
           # dolfin.DirichletBC(V, 0, boundaries, 3)]

    # Define input data
    print "Defining inputs..."
    EPS_0 = 8.854e-12       #Absolute permittivity, [Farad/metre]
    Q_E = 1.602e-19         #Elementary charge, [Coulomb/charge]
    k = 8.617e-5            #Boltzmann constant in eV/K
    T = 4                   #Temperature in Kelvin
    N_c = 2.82e25           #N_c for silicon in 1/m^3
    N_v = 1.83E25           #N_v for silicon in 1/m^3
    n_i = np.sqrt(N_c*N_v)*np.power(T,1.5)*np.exp(-1.1/(2.0*k*T))
    print "n_i = ", n_i
    # a0 = dolfin.Constant(11.68*EPS_0) #Permittivity, Si
    # a1 = dolfin.Constant(1.0*EPS_0) #Permittivity, Air
    a0 = dolfin.Constant(11.6*EPS_0) #Permittivity, Si
    a1 = dolfin.Constant(1.0*EPS_0) #Permittivity, Air
    g_T = dolfin.Constant("0.0")
    g_B = dolfin.Constant("0.0")
    g_L = dolfin.Constant("0.0")
    g_R = dolfin.Constant("0.0")
    #with T = 4 kelvin, k = 8.617e-5 ev/K, carrier charge density is ~0 inside silicon.
    f0 = dolfin.Constant(n_i) #Temperature dependent carrier density
    f1 = dolfin.Constant(0.0) #Charge density, Air
    # h_T = dolfin.Constant(np.abs(1.0/boundary_y_max))
    # h_B = dolfin.Constant(np.abs(1.0/boundary_y_min))
    h_T = dolfin.Constant(0.0)
    h_B = dolfin.Constant(0.0)
    # h_L = dolfin.Constant(0.0)
    # h_R = dolfin.Constant(0.0)

    # Define new measures associated with the interior domains and
    # exterior boundaries
    print "Defining measures..."
    dx = dolfin.dx(subdomain_data=domains)
    ds = dolfin.ds(subdomain_data=boundaries)

    # 1 left, 2 top, 3 right, 4 bottom
    # Define variational form
    print "Defining variational form..."
    F = ( a0*dolfin.inner(dolfin.grad(u), dolfin.grad(v))*dx(0) \
        + a1*dolfin.inner(dolfin.grad(u), dolfin.grad(v))*dx(1) \
        # + h_L*u*v*ds(1) \
        + h_T*u*v*ds(2) \
        # + h_R*u*v*ds(3) \
        + h_B*u*v*ds(4) \
        # - g_L*v*ds(1) \
        # - g_T*v*ds(2) \
        # - g_R*v*ds(3) \
        # - g_B*v*ds(4) \
        - f0*v*dx(0) - f1*v*dx(1) )

    # Separate left and right hand sides of equation
    print "Separating LHS and RHS..."
    a, L = dolfin.lhs(F), dolfin.rhs(F)

    # Solve problem
    u = dolfin.Function(V)
    problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
    solver = dolfin.LinearVariationalSolver(problem)
    solver.parameters['linear_solver'] = 'gmres'
    solver.parameters['preconditioner'] = 'ilu'
    prm = dolfin.parameters.krylov_solver  # short form
    prm.absolute_tolerance = 1E-25
    prm.relative_tolerance = 1E-5
    prm.maximum_iterations = 10000
    print "Solving problem..."
    solver.solve()
    print "Solve finished"

    # #Print to PVD file
    # print "Printing solution to file..."
    # file_string = "../data/Potential.pvd"
    # dolfin.File(file_string) << u

    u_P1 = dolfin.project(u, V)
    u_nodal_values = u_P1.vector()
    u_array = u_nodal_values.get_local()
    coor = mesh.coordinates()
    boundary_vals = []
    #gather the values at the Si-Air
    for i in range(len(coor)):
        if np.abs(coor[i][1] - boundary_dielectric) < dolfin.DOLFIN_EPS:
            boundary_vals += [[coor[i][0],u(coor[i][0],coor[i][1])]]
    #sort by ascending x values.
    boundary_vals.sort(key=lambda x: x[0])
    data = np.array(boundary_vals)
    x, v = data.T
    print "Amp: %f Max: %f Min: %f"%(amp_in, max(v), min(v))

    ####Potential at boundary
    print "Plotting in matplotlib..."
    # plt.figure()
    # plt.title("Potential at Si-Air boundary vs x")
    # plt.xlim(boundary_x_min, boundary_x_max)
    # plt.ylim(-0.10, 0.10)
    # plt.plot(x, v)
    # savestring = '../data/Plots/SiAirBoundary%02d.png'%(step) 
    # plt.savefig(savestring)
    # ####Potential gradient
    # V_vec = dolfin.VectorFunctionSpace(mesh, "CG", 3)
    # gradu = dolfin.project(dolfin.grad(u),V_vec)
    # plt.figure()
    # plt.ylim(boundary_y_min, boundary_y_max)
    # plt.xlim(boundary_x_min, boundary_x_max)
    # dolfin.plot(gradu, title="gradient of u")
    # ####Potential
    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(1,1,1, adjustable='box')
    # ax1.set_title("Potential due to 4-phase clocking scheme")
    label_size = 18
    tick_size = 15
    bot_trunc = 0.1
    top_trunc = 0.3
    ax1.set_ylim(bot_trunc*boundary_y_min, top_trunc*boundary_y_max)
    ax1.set_xlim(boundary_x_min, boundary_x_max)
    ax1.set_ylabel("Depth (nm)", fontsize=label_size)
    ax1.set_xlabel("Position (nm)", fontsize=label_size)
    # p = dolfin.plot(u, scalarbar=True, title="u")
    p = dolfin.plot(u, alpha=1.0)
    cbar = plt.colorbar(p, fraction = 0.046*(1.75*(bot_trunc+top_trunc)))
    cbar.set_label(label="Potential (V)", fontsize=label_size)
    cbar.ax.tick_params(labelsize=tick_size)
    ##Mesh
    dolfin.plot(mesh,alpha=1.0)
    
    air_line_data = np.array([boundary_dielectric for i in xrange(len(x))])
    ax1.plot(x, air_line_data, 'k--', linewidth=2)
    
    x_left = 0.5*elec_spacing
    for i in range(4):
        ax1.add_patch(ptc.Rectangle((x_left, -0.5*elec_height), elec_length, elec_height, facecolor="white", zorder=10))
        x_left += elec_spacing + elec_length
    locs, labels = plt.yticks()
    locs += 0.5*elec_height
    labels = []
    for loc in locs:
        labels += [str(round(loc*1e9-112.5, 2))]
    plt.yticks(locs, labels)
    
    locs, labels = plt.xticks()
    labels = []
    for loc in locs:
        labels += [str(round(loc*1e9, 2))]
    plt.xticks(locs, labels)
    ax1.tick_params(labelsize=tick_size)
    # plt.yticks()
    plt.savefig("../data/snapshot.pdf")
    plt.show()
    t = timestep*np.ones(len(x))
    print "Saving numpy arrays to file."
    np.savetxt(outfile, np.c_[t,x,v], fmt='%1.4e')
    print "Ending."
    return amp_in, max(v)

def main():
    #ONLY NEED TO CREATE THE MESH ONCE, then save it to file.
    elec_length = 25.0e-9
    elec_height = 25.0e-9 #not much effect on V.
    elec_spacing = 50.0e-9 
    boundary_x_min = 0.0e-9 #periodic, wraps around to x_max
    boundary_x_max = (4*elec_length+4*elec_spacing)
    boundary_dielectric = elec_height/2.0+100.0e-9 #surface will be placed relative to top of electrode.
    boundary_y_min = -4*boundary_dielectric #V stops changing after setting simulation to boundary +/-5*boundary_dielectric
    boundary_y_max = 4*boundary_dielectric
    mid_x = (boundary_x_max+boundary_x_min)/2.0
    mid_y = (boundary_y_max+boundary_y_min)/2.0

    params = [elec_length, elec_height, elec_spacing,
            boundary_x_min, boundary_x_max, boundary_dielectric, 
            boundary_y_min, boundary_y_max, mid_x, mid_y]
            
    # Use in-house package to define mesh
    print "Initializing mesh with MeshWriter..."
    mw = mesh_writer.MeshWriter()
    mw.resolution = min((boundary_x_max-boundary_x_min)/10.0, (boundary_y_max-boundary_y_min)/10.0)*0.5
    mw.addOuterBound([boundary_x_min,boundary_y_min], [boundary_x_max,boundary_y_max], 1)
    mw.addSeam([0.001*boundary_x_max,boundary_dielectric],[0.999*boundary_x_max,boundary_dielectric],0.1)
    fields = []
    
    #Threshold fields use indices of existing lines. Add the field directly after adding the associated seam or line. 
    fields += [mw.addTHField(0.5, 1, 0.1*boundary_dielectric, 0.5*boundary_dielectric)]
    mw.addLine([0,0],[boundary_x_max,0],0.1)
    fields += [mw.addTHField(0.5, 1, 0.75*elec_height, 3.0*elec_height)]
    mw.addMinField(fields)

    x_left = 0.5*elec_spacing
    x_right = x_left+elec_length
    mw.addSeamBox([x_left,-0.5*elec_height],\
                   [x_right,0.5*elec_height],0.1)

    x_left += elec_spacing + elec_length
    x_right = x_left+elec_length
    mw.addSeamBox([x_left,-0.5*elec_height],\
                   [x_right,0.5*elec_height],0.1)

    x_left += elec_spacing + elec_length
    x_right = x_left+elec_length
    mw.addSeamBox([x_left,-0.5*elec_height],\
                   [x_right,0.5*elec_height],0.1)

    x_left += elec_spacing + elec_length
    x_right = x_left+elec_length
    mw.addSeamBox([x_left,-0.5*elec_height],\
                   [x_right,0.5*elec_height],0.1)

    dir = "../data/Plots"
    if not os.path.exists(dir):
        os.makedirs(dir)
    with open('../data/domain.geo', 'w') as f: f.write(mw.file_string)
    subprocess.call(['gmsh -2 ../data/domain.geo -string "General.ExpertMode=1;"'+\
                     ' -string "Mesh.CharacteristicLengthFromPoints=0;"'+\
                     ' -string "Mesh.CharacteristicLengthExtendFromBoundary=0;"'], shell=True) #Expert mode to suppress warnings about fine mesh
    subprocess.call(['dolfin-convert ../data/domain.msh ../data/domain.xml'], shell=True)
    
    results = []
    # testvals = np.linspace(0.1, 0.2, 11)
    # testvals = np.array([0.12]) #works fine for the 2_1_1 and 1_2_1 and 2_2_1 case.
    # testvals = np.array([0.16]) #works when elec_length, elec_spacing, dielectric_boundary are the same
    # testvals = np.array([0.32]) #works for 1_1_2
    testvals = np.array([0.6])
    offsets = np.array([-1.05])
    # timesteps = np.linspace(0,2*np.pi, 40, endpoint=False)
    timesteps = [0]
    outfile = file("../data/surfacepotentials.txt", "w")
    for i in range(len(testvals)):
        for j in range(len(timesteps)):
            results += [worker(testvals[i], offsets[i], timesteps[j],j, params, outfile)]
    res_x, res_v = zip(*results)    
    print max(res_v), min(res_v)
    print "Wobble: ", (max(res_v) - min(res_v))/((max(res_v)+min(res_v))/2.0)*100, "%"
    
    #read from file
    # outfile.close()
    # outfile = file("../data/surfacepotentials.txt", "r")
    # t, x, v = np.loadtxt(outfile, unpack=True)
    # print t
    # print x
    # print v
main()

#For making gifs
# subprocess.call(['convert -delay 10 -loop 0 ../data/Plots/*.png ../data/Plots/SiAirBoundary.gif'], shell=True)

