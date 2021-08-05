
#=

Authors: Judd Mehr,

Preliminary Design functions

=#

export addairfoilflaps, saveairfoilpolar, airfoilcomp


############################
#######    General    ######
############################

"""
    stabilityderivatives(a::prelim, ns, nc, Vinf)

Outputs all stability derivatives using vortex lattice

**Inputs**

A struct of type `prelim`
`ns::Float64` : number of panels along the span
`nc::Float64` : number of panels along the chord
`Vinf::Float64` : freestream speed
"""

function stabilityderivatives(a::prelim, ns, nc, Vinf)
# geometry (right half of the wing)
    xle = [0.0, a.thickness] #thickness ofr sweep?
    yle = [0.0, a.span/2]
    zle = [0.0, 0.0] #tan(sweep)*span/2?
    chord = [a.rootchord, a.tipchord]
    theta = [a.angleofattack*pi/180, a.twistangle*pi/180]
    phi = [a.dihedral, a.polyhedral]
    fc = fill((xc) -> 0, 2) # camberline function for each section  ??

    # discretization parameters
    spacing_s = Uniform()
    spacing_c = Uniform()

    # reference parameters
    Sref = .5*a.span*(chord[1] + chord[2])
    cref = .5*(chord[1] + chord[2])
    bref = a.span
    rref = [a.cg, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = a.angleofattack*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all surfaces
    surfaces = [surface]

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = true

    # perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

    dCF, dCM = stability_derivatives(system)

    CDa, CYa, CLa = dCF.alpha
    Cla, Cma, Cna = dCM.alpha
    CDb, CYb, CLb = dCF.beta
    Clb, Cmb, Cnb = dCM.beta
    CDp, CYp, CLp = dCF.p
    Clp, Cmp, Cnp = dCM.p
    CDq, CYq, CLq = dCF.q
    Clq, Cmq, Cnq = dCM.q
    CDr, CYr, CLr = dCF.r
    Clr, Cmr, Cnr = dCM.r

    println("for α:", alpha*180/pi)
    println("CLa:", CLa)
    println("CLb:", CLb)
    println("CYa:", CYa)
    println("CYb:", CYb)
    println("Cla:", Cla)
    println("Clb:", Clb)
    println("Cma:", Cma)
    println("Cmb:", Cmb)
    println("Cna:", Cna)
    println("Cnb:", Cnb)
    println("CLp:", CLp)
    println("CLq:", CLq)
    println("CLr:", CLr)
    println("CYp:", CYp)
    println("CYq:", CYq)
    println("CYr:", CYr)
    println("Clp:", Clp)
    println("Clq:", Clq)
    println("Clr:", Clr)
    println("Cmp:", Cmp)
    println("Cmq:", Cmq)
    println("Cmr:", Cmr)
    println("Cnp:", Cnp)
    println("Cnq:", Cnq)
    println("Cnr:", Cnr)

return dCF, dCM
end



############################
#######    Sizing    #######
############################



############################
####    Aerodynamics    ####
############################
"""
    addairfoilflaps(filename,datapath,xc,yt,angle;npan=140)

Add a flap to airfoil geometry at a given percent chord and percent thickness at that location.

*Inputs:*
- filename::String : file name of airfoil geometry (format: column 1 = x/c, column 2 = y/c)
- suffix::String : file type of file (e.g. .dat, .csv, etc.)
- datapath::String : the path where the file can be found
- xc::Float64 : x/c location where flap is to be defined (e.g. 0.7 is at 70% chord)
- yt::Float64 : y/t relative to thickness at x/c location (e.g. 0.5 is halfway between upper and lower surfaces).
- angle::Float64 : angle of flap deflection, in degrees. + is down, - is up.

*Keyword Arguments:*
- npan::Int64 : number of panels to use for geometry
- savefile::String : path and file name where to save output geometry.

*Outputs:*
- geometry file with flap added, default file name is datapath*filename*"$(angle)deg_$(xc*100)chord.dat"
"""
function addairfoilflaps(filename, suffix, datapath, xc, yt, angle;
                npan=140,
                savefile=nothing
                )

    #Get original airfoil data
    xy = convert(Array,DataFrames.DataFrame!(CSV.File(datapath*filename*suffix,header=false,datarow=2)))

    #Pane data so it's standardized
    Xfoil.set_coordinates(xy[:,1],xy[:,2])
    Xfoil.pane(npan=npan)
    x = Xfoil.get_globals().x[1:npan]
    y = Xfoil.get_globals().y[1:npan]

    #identify points to rotate
    flapidx = findall(x->x>=xc,x)
    flaparr = [x[flapidx] y[flapidx]]

    #calculate y-axis rotation point
    #find halfway point of flapidx
    halfway = Int(length(flapidx)/2)

    #then grab y points there and find distance between
    ylo = y[flapidx[halfway]]
    yhi = y[flapidx[halfway+1]]
    #then calculate actual y point based on relative thickness input
    ypt = ylo + yt*(yhi-ylo)

    #translate to origin from rotation point
    flaparr .-= [xc ypt]

    #rotate
    rotmat = [cosd(-angle) -sind(-angle);
              sind(-angle)  cosd(-angle)]

    rotarr = rotmat*flaparr'
    rotarr = rotarr'

    #translate back
    rotarr .+= [xc ypt]

    #replace previous flap locations with rotated ones
    xnew = copy(x)
    xnew[flapidx] = rotarr[:,1]
    ynew = copy(y)
    ynew[flapidx] = rotarr[:,2]

    #delete the points that dip below (or above) old y values to avoid kinks.
    if sign(angle) < 0 #if flap is up
        rmidx = findall(x->x<y[flapidx[halfway]],ynew[flapidx[1:halfway]])
        deleteat!(ynew,flapidx[rmidx])
        deleteat!(xnew,flapidx[rmidx])
    else
        rmidx = findall(x->x>y[flapidx[halfway+1]],ynew[flapidx[halfway+1]:end])
        deleteat!(ynew,flapidx[halfway.+rmidx])
        deleteat!(xnew,flapidx[halfway.+rmidx])
    end

    #savefile:
    if savefile === nothing
        savefile =  datapath*filename*"$(angle)deg_$(xc*100)chord.dat"
    end

    f = open(savefile, "w")

    string = "x/c y/c\n"

    write(f, string)

    for i=1:length(xnew)

        string = "$(xnew[i]) $(ynew[i])\n"

        write(f, string)

    end

    close(f)

end


"""
    saveairfoilpolar(filename, datapath, Re, alpha, cl, cd, cdp, cm, conv)

Save outputs from Xfoil.jl analysis.

*Inputs:*
- filename::String : file name of save file.
- datapath::String : path where to save file.
- Re::Float64 : Reynolds number.
- alpha::Array{Float64} : angles of attack.
- cl::Array{Float64} : lift coefficients.
- cd::Array{Float64} : drag coefficients.
- cdp::Array{Float64} : parasitic drag coefficients.
- cm::Array{Float64} : moment coefficients.
- conv::Array{Bool} : convergence flags.

*Keyword Arguments:*
- savefile::String : path and file where to save output.

*Outputs:*
- saved file with location and name. Default: datapath*filename*"$(@sprintf("%.2E", Re)).csv"  File will always be saved with "$(@sprintf("%.2E", Re)).csv" tagged on the end.
"""
function saveairfoilpolar(filename, datapath, Re, alpha, cl, cd, cdp, cm, conv;
                    savefile = nothing)

    if savefile === nothing
        savefile = datapath*filename*"$(@sprintf("%.2E", Re)).csv"
    else
        savefile *= "$(@sprintf("%.2E", Re)).csv"
    end

    f = open(savefile, "w")

    string = "alpha,cl,cd,cdp,cm,conv\n"

    write(f, string)

    for i=1:length(alpha)

        string = "$(alpha[i]),$(cl[i]),$(cd[i]),$(cdp[i]),$(cm[i]),$(conv[i])\n"

        write(f, string)

    end

    close(f)

end

"""
    runxfoil(filenames, datapath, Revalues, alpharange=collect(-5:20); savefile)

Run xfoil for give airfoil geometry and Reynolds numbers, and pass data to saveairfoilpolar().

*Inputs:*
- filenames::Array{String} : geometry file names of airfoils to analyze.
- datapath::String : path where files are located.
- alpharange::Array{Float64} : angles of attack (in degrees) to analyze.

*Keyword Arguments:*
- savefile::Array{String} : names of files (including path) to save the polar data. (File will always be saved with "$(@sprintf("%.2E", Re)).csv" tagged on the end.)
"""
function runxfoil(filenames, datapath, Revalues, alpharange=collect(-5:20);
            savefile = nothing)

    for i=1:length(filenames)

        #Load Airfoil
        xy = convert(Array,DataFrames.DataFrame!(CSV.File(datapath*filenames[i],header=false,datarow=2)))

        #Loop through Reynolds numbers
        for j=1:length(Revalues)

            #Run Xfoil
            cl, cd, cdp, cm, conv = Xfoil.alpha_sweep(xy[:,1],xy[:,2], alpharange, Revalues[j])

            #Save data so you don't have to re-run xfoil
            saveairfoilpolar(filenames[i],datapath,Revalues[j],alpharange,cl,cd,cdp,cm,conv)

        end

    end

end



############################
#####    Propulsion    #####
############################

"""
    motorprop(a::prelim, Kvr, i0, R, v, d, num)


*Inputs:*
A struct of type 'prelim'
`Kv::Float64` : Motor electricity constant
`i0::Float64` : no load current
`KR:Float64` : Resistance
`v::Float64` : battery voltage
`d::Float64` : propellar diameter
`num::Float64` : number of motors/propellars

*Outputs a graph of efficiency vs rotational speed and four other graphs of force, efficiency, RPMs, and current vs airspeed


"""


function motorprop(a::prelim, Kvr, i0, R, v, d, num)
    #plane
    b = a.span
    S = b*(a.rootchord + a.tipchord)/2
    println(S,a.weight)
    tc = a.tipchord/a.rootchord
    AR = b^2/S
    fr = a.finenessratio
    Λ = a.sweep

    conceptual = concept(S, a.weight, "conventional tail", "front propeller")  #object of type 'concept' from aircrafttypes.jl
    df = a.df
    # motor
    m_KvRPM = Kvr
    m_i0 = i0  # no load current (Amps)
    m_R = R # resistance (Ohms)
    m_vbatt = v  # battery voltage (Volts)

    # esc
    throttle = .7  # number between 0 and 1 - not exactly throttle depending on esc
                     # more precisely it is applied voltage relative to battery capacity: v/vbatt

    # prop
    # propeller data with columns J, CT, CP, eta (same format and normalization as from the UIUC site)
    p_data = [
    0.137   0.1426   0.0835   0.234
    0.173   0.1400   0.0834   0.290
    0.207   0.1368   0.0831   0.341
    0.245   0.1329   0.0828   0.393
    0.280   0.1288   0.0823   0.438
    0.319   0.1232   0.0811   0.485
    0.357   0.1173   0.0797   0.525
    0.397   0.1105   0.0776   0.565
    0.432   0.1046   0.0758   0.597
    0.466   0.0986   0.0737   0.623
    0.508   0.0908   0.0706   0.653
    0.539   0.0848   0.0680   0.672
    0.579   0.0776   0.0649   0.692
    0.610   0.0725   0.0627   0.706
    0.653   0.0652   0.0590   0.722
    0.685   0.0590   0.0557   0.725
    0.726   0.0520   0.0519   0.727
    ]

    p_D = d * 0.0254  # prop diamter (m)
    p_num = num  # number of motor/props

    # atmosphere
    rho = 1.225  # atmospheric density (kg/m^3)
    μ = 1.81e-5
    # other aircraft inputs

    Vdesign = designspeed(conceptual)  # design flight speed. Function from conceptiual.jl
    println(Vdesign)
    q = .5*rho*Vdesign^2

    ac_CDp =  dragparasitic(S, Vdesign, tc, AR, fr, Λ, rho, μ)/(q*S) # parasitic drag coefficient. Function from conceptual.jl
    ac_einv = .98(1-2(df/b)^2)  # inviscid span efficiency (NOT Oswald efficiency factor, that is computed internally)
    ac_Sref = S  # reference area for CDp (m^2)
    ac_b = b  # wing span (m)
    ac_mass = (a.weight/9.81)  # mass of airplane (kg)



    # ---------------------------------------



    m_Kv = m_KvRPM * pi/30
    p_rho = rho
    ac_rho = rho


    # -------- functions ------------

    function motor(m_vbatt, m_Kv, m_i0, m_R, Omega, throttle)

        vb = m_vbatt*throttle  # assumes linear
        i = (vb .- Omega/m_Kv)/m_R  # current
        Q = (i .- m_i0)/m_Kv  # torque
        Pout = Q.*Omega  # power out
        eta = Pout/(i*vb)  # efficiency

        return i, Q, eta

    end


    function prop(p_data, p_D, p_rho, Omega, V)
        # unit conversion
        n = Omega/(2*pi)

        # extract data
        J_data = p_data[:, 1]
        CT_data = p_data[:, 2]
        CP_data = p_data[:, 3]
        eta_data = p_data[:, 4]

        # interpolate off of charts
        J = V ./(n*p_D)  # advance ratio
        CT = linear(J_data, CT_data, J)
        CP = linear(J_data, CP_data, J)
        eta = linear(J_data, eta_data, J)

        CQ = CP/(2*pi)

        T = CT * p_rho .* n.^2 * p_D^4
        Q = CQ * p_rho .* n.^2 * p_D^5
        return Q, T, eta
    end

    function operatingpoint(m_vbatt, m_Kv, m_i0, m_R, p_data, p_D, p_rho, throttle, V)
        # solve for when torque's are equal
        function difference(Omega_g)
            _, Qm, _ = motor(m_vbatt, m_Kv, m_i0, m_R, Omega_g, throttle)
            Qp, _, _ = prop(p_data, p_D, p_rho, Omega_g, V)
            residual = Qm-Qp
            return residual
        end

        # range
        Jmin = p_data[1, 1]
        Jmax = last(p_data[:,1])
        Omegamax = 2*pi*V/(p_D*Jmin)-1e-3
        Omegamin = 2*pi*V/(p_D*Jmax)+1e-3
        if difference(Omegamin) * difference(Omegamax) < 0  # solution exists
            Omega = find_zero(difference, (Omegamin+1e-3,Omegamax-1e-3)) #find root
            # evaluate models at this Omega
            i, _, etam = motor(m_vbatt, m_Kv, m_i0, m_R, Omega, throttle)
            Q, T, etap = prop(p_data, p_D, p_rho, Omega, V)
            T = T * p_num

        else  # no solution
            i = NaN
            Q = NaN
            T = NaN
            Omega = NaN
            etam = NaN
            etap = NaN
        end
        return Q, T, i, Omega, etam, etap
    end


    function drag(ac_mass, rho, b, S, ac_einv, ac_CDp, V)

        # compute drag (and thus required thrust)
        g = 9.81
        W = ac_mass*g
        q = 0.5*rho*V^2  # dynamic pressure
        AR = b^2/S  # aspect ratio
        e = 1.0 / (1.0/ac_einv + 0.38*ac_CDp*pi*AR)  # Oswald efficiency factor
        D = ac_CDp*q*S + W^2/(q*pi*b^2*e)  # total drag
        return D
    end


    # --- run a sweep versus airspeed -----
    nsweep = 100


    # limits
    Vmin = 2
    Vmax = 1.5*Vdesign

    # initialize vectors
    Vvec = LinRange(Vmin, Vmax, nsweep)
    Qvec = zeros(nsweep, 1)
    Tvec = zeros(nsweep, 1)
    Dvec = zeros(nsweep, 1)
    ivec = zeros(nsweep, 1)
    Omegavec = zeros(nsweep, 1)
    etamvec = zeros(nsweep, 1)
    etapvec = zeros(nsweep, 1)

    for i = 1:nsweep
        Qvec[i], Tvec[i], ivec[i], Omegavec[i], etamvec[i], etapvec[i] = operatingpoint(m_vbatt, m_Kv, m_i0, m_R, p_data, p_D, p_rho, throttle, Vvec[i])
        Dvec[i] = drag(ac_mass, rho, b, S, ac_einv, ac_CDp, Vvec[i])
    end


    # ---- sweep performance vs Omega -----
    nomega = 150

    # design rotation speed
    _, _, _, Omega_d, _, _ = operatingpoint(m_vbatt, m_Kv, m_i0, m_R, p_data, p_D, p_rho, throttle, Vdesign)

    # limits of rotation speed for prop
    Jmin = p_data[1, 1]
    Jmax = last(p_data[:,1])
    Omegamax = 2*pi*Vdesign/(p_D*Jmin)-1e-3
    Omegamin = 2*pi*Vdesign/(p_D*Jmax)+1e-3
    Omegavecprop = LinRange(Omegamin, Omegamax, nomega)

    # limits of rotation speed for motor
    Omegamax = 0.9*m_vbatt*m_Kv
    Omegavecmotor = LinRange(0, Omegamax, nomega)

    # initialize
    etam = zeros(nomega, 1)
    etap = zeros(nomega, 1)

    for i = 1:nomega
        _, _, etam[i] = motor(m_vbatt, m_Kv, m_i0, m_R, Omegavecmotor[i], throttle)
        _, _, etap[i] = prop(p_data, p_D, p_rho, Omegavecprop[i], Vdesign)
    end


    return Omegavecmotor, etam, etap, Omegavecprop, Omega_d, Vvec, Tvec, Dvec, etamvec, etapvec, Omegavec, ivec, Vdesign
end

function plotmotorprop(a::prelim, Kvr, i0, R, v, d, num)
    Omegavecmotor, etam, etap, Omegavecprop, Omega_d, Vvec, Tvec, Dvec, etamvec, etapvec, Omegavec, ivec, Vdesign = motorprop(a, Kvr, i0, R, v, d, num)
    p1 = plot(Omegavecmotor*30/pi, etam, xlabel = "Omega (RPM)", ylabel = "efficiency", label = "motor", ylims =  (0, 1), xlims = (0, max(last(Omegavecmotor), last(Omegavecprop)*30/pi)))
    plot!(Omegavecprop*30/pi, etap, label = "prop")
    vline!([Omega_d*30/pi], label = "design point")

    p2 = plot(Vvec, Tvec, label = "Thrust", xlabel = "airspeed (m/s)", ylabel = "force (N)", ylims = (0, maximum(Tvec)))
    plot!(Vvec, Dvec, label = "Drag")
    vline!([Vdesign], label = "design speed")

    p3 = plot(Vvec, etamvec, label = "motor", xlabel = "airspeed (m/s)", ylabel = "efficiency", ylims = (0,1))
    plot!(Vvec, etapvec, label = "prop")
    plot!(Vvec, etamvec.*etapvec, label = "total")
    vline!([Vdesign], label = "design speed")

    p4 = plot(Vvec, Omegavec*30/pi, label = "Ω", xlabel = "airspeed (m/s)", ylabel = "rotation rate (RPM)")
    vline!([Vdesign], label = "design speed")

    p5 = plot(Vvec, ivec, label = "V", xlabel = "airspeed (m/s)", ylabel = "current (amps)")
    vline!([Vdesign], label = "design speed")

    #plot(p1,p2,p3,p4, or p5)
    plot(p1)

end

"""
flightenvelope(a::prelim)

Inputs
A struct of type 'prelim'
`ρ::Float64` : air density, assumed to be 1.225

Outputs the altitude and speed limits of the aircraft

"""

function flightenvelope(a::prelim)
    #aircraft properties
    b = a.span
    S = b*(a.rootchord + a.tipchord)/2
    tc = a.tipchord/a.rootchord
    AR = b^2/S
    fr = a.finenessratio
    Λ = a.sweep
    w = a.weight
    Clmax = a.Clmax
    md = a.machdive
    qm = a.qmax

    ceiling = 4000 #legal limit for RC aircraft or ceilingcalc(s, tc, AR, fr, Λ, ρ) for larger aircraft

    alts = range(1, stop = 40000, length = 5000)
    vs = range(1, stop = 50, length = 5000)

    #atmosphere properties
    ρ = zeros(length(alts))
    a = zeros(length(alts))
    for (i,_) in enumerate(alts)
        ρ[i], _, a[i] = atmospherefit(alts[i])
    end

    #Vstall
    function Vstall(w, S, Clmax, ρ)
        Vstall = zeros(length(ρ))
        for (i,_) in enumerate(ρ)
            Vstall[i] = sqrt(2*w/(ρ[i]*S*Clmax))
        end
        return Vstall
    end


    """
    ceilingcalc() calculates the true ceiling of an aircraft.  In this case we assume the ceiling is the legal altitude limit of RC aircraft
    """
    function ceilingcalc(a::prelim, vs, tc, AR, fr, Λ, ρ)
        conceptual = concept(a.area, a.weight, "conventional", "front singular")

        #velocity sweep
        d = zeros(length(ρ))
        for (i,_) in enumerate(ρ)
            v = Vstall(a,ρ)
            d[i] = dragcalculator(conceptual, v, tc, AR, fr, Λ, ρ[i])
            #more calculations necessary to determine ceiling
        end

    end

    function lowlimit(a, ρ, md, ceiling)
        lowlim = zeros(length(a))
        eas = sqrt(2*qm/1.225)
        for (j,_) in enumerate(a)
            i = length(alts) - j + 1
            lowlim[i] = a[i]*md
            q = .5*ρ[i]*lowlim[i]^2
            if q > qm
                lowlim[i] = eas/sqrt(ρ[i]/1.225)
            end
        end
        return lowlim
    end

    Vstall = Vstall(w,S,Clmax,ρ)
    ceilingIndex = findfirst(>=(ceiling), alts)
    plot(Vstall[1:ceilingIndex],alts[1:ceilingIndex], ylims = (0, ceiling[1]+ceiling[1]/2), label = "Vstall", xlabel = "speed", ylabel = "altitude")


    lowlim = lowlimit(a, ρ, md, ceiling)
    plot!(lowlim[1:ceilingIndex], alts[1:ceilingIndex], label = "speed limit")


    endCeiling = findfirst(>=(lowlim[ceilingIndex]), Vstall)
    ceilings = repeat([ceiling], length(vs))
    plot!(Vstall[ceilingIndex:endCeiling], ceilings[ceilingIndex:endCeiling], label = "ceiling")



end

"""
Vndiagram(a::prelim)

Inputs
A struct of type 'prelim'

Outputs the load factor and speed limits of the aircraft

"""


function Vndiagram(a::prelim, ρ = 1.225)

    #aircraft properties
    w = a.weight
    b = a.span
    S = b*(a.rootchord + a.tipchord)/2
    qm = a.qmax
    vd = sqrt(2*qm/1.225)
    vc = vd/1.25
    Clmax = a.Clmax

    #stall
    vs = range(1, stop = vd + vd*1.5, length = 1000)
    nstall = zeros(length(vs))
    for (i,v) in enumerate(vs)
        nstall[i] = (.5*Clmax*ρ*S*v^2)/w
    end
    lim1 = findfirst(>=(a.maxloadfactor), nstall)
    println("Sherlock\n\tlim1:$lim1")
    plot(vs[1:lim1], nstall[1:lim1], label = "stall", xlabel = "speed", ylabel = "load factor", ylims = (-2, a.maxloadfactor + a.maxloadfactor*1.5))


    #max maneuver load factor
    nmax = repeat([a.maxloadfactor], length(vs))
    lim2 = findfirst(>=(vd), vs)
    plot!(vs[lim1:lim2], nmax[lim1:lim2], label = "nmmax")


    #speed limit
    ns = range(0,stop = nmax[1], length = length(vs))
    border = repeat([vd],length(ns))
    plot!(border, ns, label = "vd")

    #negative stall
    negnstall = zeros(length(vs))
    for (i,v) in enumerate(vs)
        negnstall[i] = -(.5*Clmax*ρ*S*v^2)/w
    end
    neglim1 = findfirst(<=(-1), negnstall)
    plot!(vs[1:neglim1], negnstall[1:neglim1], label = "stall")

    #minimum load factor
    nmin = repeat([-1], length(vs))
    neglim2 = findfirst(>=(vc), vs)
    plot!(vs[neglim1:neglim2], nmin[neglim1:neglim2], label = "nmin")

    #linear
    linear  = zeros(length(vs))
    for (i,v) in enumerate(vs)
        linear[i] = -(v-vd)/(vc-vd)
    end
    plot!(vs[neglim2:lim2], linear[neglim2:lim2])


end

############################
#####    Structures    #####
############################
