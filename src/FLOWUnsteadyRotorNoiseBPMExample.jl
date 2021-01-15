module FLOWUnsteadyRotorNoiseBPMExample

import FLOWUnsteady

import DelimitedFiles
import Makie
import AbstractPlotting.MakieLayout
import AbstractPlotting

function doit()
    # Aliases
    uns = FLOWUnsteady
    vlm = FLOWUnsteady.vlm
    noise = FLOWUnsteady.noise
    gt = FLOWUnsteady.gt

    # Create temps folder
    if ispath("data/temps")==false; mkpath("data/temps"); end;


    save_path = "data/temps/dji9443_single_bemnoise_00/"  # Where to save the simulation


    # ------------ PARAMETERS --------------------------------------------------
    # Rotor geometry
    rotor_file      = "DJI9443.csv"        # Rotor geometry
    data_path       = uns.def_data_path    # Path to rotor database
    pitch           = 0.0                  # (deg) collective pitch of blades
    n               = 50                   # Number of blade elements
    CW              = true                 # Clock-wise rotation
    xfoil           = true                 # Whether to run XFOIL


    # Read radius of this rotor and number of blades
    R, B            = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

    # Simulation parameters
    RPM             = 5400                 # RPM
    J               = 0.0001               # Advance ratio Vinf/(nD)
    rho             = 1.071778             # (kg/m^3) air density
    mu              = 1.85508e-5           # (kg/ms) air dynamic viscosity
    speedofsound    = 342.35               # (m/s) speed of sound
    ReD             = 2*pi*RPM/60*R * rho/mu * 2*R   # Diameter-based Reynolds number

    magVinf         = J*RPM/60*(2*R)
    Vinf(X,t)       = magVinf*[1.0, 0, 0]  # (m/s) freestream velocity

    # Aerodynamic solver
    # VehicleType   = uns.UVLMVehicle      # Unsteady solver
    VehicleType     = uns.QVLMVehicle      # Quasi-steady solver

    # Solver parameters
    nrevs           = 10                   # Number of revolutions in simulation
    nsteps_per_rev  = 120                  # Time steps per revolution
    ttot            = nrevs/(RPM/60)       # (s) total simulation time
    nsteps          = nrevs*nsteps_per_rev # Number of time steps

    # (Unsteady solver parameters---not used in this example)
    p_per_step      = 2                    # Sheds per time step
    lambda          = 2.125                # Core overlap
    overwrite_sigma = lambda * 2*pi*R/(nsteps_per_rev*p_per_step) # Smoothing core size
    surf_sigma      = R/10                 # Smoothing radius of lifting surface
    vlm_sigma       = surf_sigma           # Smoothing radius of VLM
    shed_unsteady   = true                 # Shed particles from unsteady loading
                                           # Max particles for memory pre-allocation
    max_particles   = ((2*n+1)*B)*nrevs*nsteps_per_rev*p_per_step

    # OUTPUT OPTIONS
    run_name        = "singlerotor"
    nsteps_save     = 1                    # Save vtks every this many steps
    save_wopwopin   = true                 # Generate inputs for PSU-WOPWOP
    prompt          = false                # Whether to promp the user
    verbose         = true
    plot_disc       = true                 # Plot blade discretization for debugging




    # ------------ SIMULATION SETUP --------------------------------------------
    # Generate rotor
    rotor = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n, CW=CW, ReD=ReD, altReD=[RPM, J, mu/rho],
                                            verbose=verbose, xfoil=xfoil,
                                            data_path=data_path,
                                            plot_disc=plot_disc)

     # ----- VEHICLE DEFINITION
     # System of all FLOWVLM objects
     system = vlm.WingSystem()
     vlm.addwing(system, "Rotor1", rotor)

     # Systems of rotors
     rotors = vlm.Rotor[rotor]   # Defining this rotor as its own system
     rotor_systems = (rotors,)

     # Wake-shedding system (don't include the rotor if quasi-steady vehicle)
     wake_system = vlm.WingSystem()

     if VehicleType != uns.QVLMVehicle
         vlm.addwing(wake_system, "Rotor1", rotor)
     else
         # Mute colinear warnings. This is needed since the quasi-steady solver
         #   will probe induced velocities at the lifting line of the blade
         uns.vlm.VLMSolver._mute_warning(true)
     end

     # Define vehicle object
     vehicle = VehicleType(   system;
                                 rotor_systems=rotor_systems,
                                 wake_system=wake_system
                              )


     # ----- MANEUVER DEFINITION
     RPM_fun(t) = 1.0                # RPM (normalized by reference RPM) as a
                                     # function of normalized time
     angle = ()                      # Angle of each tilting system (none in this case)
     sysRPM = (RPM_fun, )            # RPM of each rotor system
     Vvehicle(t) = zeros(3)          # Translational velocity of vehicle over Vcruise
     anglevehicle(t) = zeros(3)      # (deg) angle of the vehicle

     # Define Maneuver object
     maneuver = uns.KinematicManeuver(angle, sysRPM, Vvehicle, anglevehicle)

     # Plot maneuver path and controls
     uns.plot_maneuver(maneuver; vis_nsteps=nsteps)


     # ----- SIMULATION DEFINITION
     RPMref = RPM
     Vref = 0.0

     tinit = 0.0                                  # (s) initial time
     Vinit = Vref*maneuver.Vvehicle(tinit/ttot)   # (m/s) initial vehicle velocity
     angle1 = maneuver.anglevehicle(tinit/ttot)   # (rad/s) initial vehicle angular velocity
     angle2 = maneuver.anglevehicle(tinit/ttot + 1e-12)
     Winit = pi/180 * (angle2-angle1)/(ttot*1e-12)

     simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                 Vinit=Vinit, Winit=Winit, t=tinit)


     # ----- MONITOR DEFINITION
     monitor = uns.generate_monitor_rotors( rotors, J, rho, RPM, nsteps;
                                         t_scale=RPM/60,        # Time scaling factor for plots
                                         t_lbl="Revolutions",   # x-axis label
                                         # OUTPUT OPTIONS
                                         save_path=save_path,
                                         run_name="rotor",
                                         figname="monitor_rotor",
                                         disp_conv=verbose,)


     # ------------ RUN SIMULATION ----------------------------------------------
     pfield = uns.run_simulation(simulation, nsteps;
                                       # SIMULATION OPTIONS
                                       Vinf=Vinf,
                                       rho=rho,
                                       mu=mu,
                                       sound_spd=speedofsound,
                                       # SOLVERS OPTIONS
                                       p_per_step=p_per_step,
                                       overwrite_sigma=overwrite_sigma,
                                       vlm_sigma=vlm_sigma,
                                       surf_sigma=surf_sigma,
                                       max_particles=max_particles,
                                       shed_unsteady=shed_unsteady,
                                       extra_runtime_function=monitor,
                                       # OUTPUT OPTIONS
                                       save_path=save_path,
                                       run_name=run_name,
                                       nsteps_save=nsteps_save,
                                       save_wopwopin=save_wopwopin,
                                       prompt=prompt,
                                       verbose=verbose)
    
    # ------------ PARAMETERS --------------------------------------------------
    # NOTE: Make sure that these parameters match what was used in the 
    #       aerodynamic solution.

    # Rotor geometry
    rotor_file      = "DJI9443.csv"        # Rotor geometry
    data_path       = uns.def_data_path    # Path to rotor database
    pitch           = 0.0                  # (deg) collective pitch of blades
    n               = 50                   # Number of blade elements
    CW              = true                 # Clock-wise rotation


    # Read radius of this rotor and number of blades
    R, B            = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

    # Simulation parameters
    RPM             = 5400                 # RPM
    J               = 0.0001               # Advance ratio Vinf/(nD)
    rho             = 1.071778             # (kg/m^3) air density
    mu              = 1.85508e-5           # (kg/ms) air dynamic viscosity
    speedofsound    = 342.35               # (m/s) speed of sound

    magVinf         = J*RPM/60*(2*R)
    Vinf(X,t)       = magVinf*[1.0, 0, 0]  # (m/s) freestream velocity

    # BPM parameters
    noise_correction= 0.65                 # Calibration parameter
    TE_thickness    = 16.0                 # (deg) trailing edge thickness
    freq_bins       = uns.BPM.default_f    # Frequency bins (default is one-third octave band)

    # Observer definition: Circular array of microphones
    sph_R           = 1.905                # (m) radial distance from rotor hub
    sph_nR          = 0
    sph_nphi        = 0
    sph_ntht        = 360                  # Number of microphones
    sph_thtmin      = 0                    # (deg) first microphone's angle
    sph_thtmax      = 360                  # (deg) last microphone's angle
    sph_phimax      = 180
    sph_rotation    = [90, 0, 0]           # Rotation of grid of microphones

    # Observer definition: Single microphone
    Rmic = 1.905                           # (m) radial distance from rotor hub
    anglemic = 90*pi/180                   # (rad) microphone angle from plane of rotation (- below, + above)
                                           # 0deg is at the plane of rotation, 90deg is upstream
    microphoneX = nothing                  # Comment and uncomment this to switch from array to single microphone
    # microphoneX = Rmic*[-sin(anglemic), cos(anglemic), 0]

    # OUTPUT OPTIONS
    prompt          = true                 # Whether to promp the user
    verbose         = false
    plot_disc       = false                # Plot blade discretization for debugging

    # ------------ GENERATE GEOMETRY -------------------------------------------
    # Generate rotor
    rotor = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n, CW=CW, ReD=0.0,
                                            verbose=verbose, xfoil=false,
                                            data_path=data_path,
                                            plot_disc=plot_disc)
    rotors = vlm.Rotor[rotor]

    # ------------ RUN BPM -----------------------------------------------------
    uns.run_noise_bpm(rotors, RPM, Vinf, rho, mu, speedofsound,
                                    save_path;
                                    # ---------- OBSERVERS -------------------------
                                    sph_R=sph_R,
                                    sph_nR=sph_nR, sph_ntht=sph_ntht,
                                    sph_nphi=sph_nphi, sph_phimax=sph_phimax,
                                    sph_rotation=sph_rotation,
                                    sph_thtmin=sph_thtmin, sph_thtmax=sph_thtmax,
                                    microphoneX=microphoneX,
                                    # ---------- BPM OPTIONS -----------------------
                                    noise_correction=noise_correction,
                                    TE_thickness=TE_thickness,
                                    freq_bins=freq_bins,
                                    # ---------- OUTPUT OPTIONS --------------------
                                    prompt=prompt
                                    );

    # Read the data from the files FLOWUnsteady writes out.
    freq, spl = get_bpm_spl(save_path)

    # Blade passing frequency in Hertz.
    bpf = RPM/60*B

    # Get the θ=-45° spl value.
    freq = freq[:, 1]
    spl = spl[:, -(-45) + 181]

    # Now we'll plot the BPM stuff.
    scene, layout = AbstractPlotting.MakieLayout.layoutscene()

    ax1 = layout[1, 1] = AbstractPlotting.MakieLayout.LAxis(scene, xlabel="freq/bpf", ylabel="SPL")
    line = AbstractPlotting.lines!(ax1, log10.(freq./bpf), spl)
    AbstractPlotting.MakieLayout.xlims!(ax1, log10(6*10^-1), log10(2*10^2))
    ax1.xticks = [0, 1, 2]
    AbstractPlotting.MakieLayout.ylims!(ax1, 0.0, 60.0)
    ax1.yticks = 0:10:60

    AbstractPlotting.save("spl_bpm.png", scene)

    return rotor
end

function get_bpm_spl(save_path)
    # Need to process the files a bit before I can read it with DelimitedFiles.
    # I'd like to 
    #   * delete/skip the first 4 lines
    #   * delete any leading and trailing spaces
    #   * replace any remaining run of spaces with a single comma
    s = open("$(save_path)/frequencies.tec", "r") do f
        read(f, String)
    end
    s = replace(s, r"^  *"=>"")  # Delete any leading spaces.
    s = replace(s, r"  *\n"=>"\n")  # Delete any spaces before a newline.
    s = replace(s, r"\n  *"=>"\n")  # Delete any spaces after a newline.
    s = replace(s, r"  *$"=>"")  # Delete any trailing spaces.
    s = replace(s, r"  *"=>",")  # Replace any remaining run of spaces with a single comma.
    io = IOBuffer(s)
    freq = DelimitedFiles.readdlm(io, ',', Float64, '\n'; skipstart=4, skipblanks=true)

    s = open("$(save_path)/spl_spectrum.tec", "r") do f
        read(f, String)
    end
    s = replace(s, r"^  *"=>"")  # Delete any leading spaces.
    s = replace(s, r"  *\n"=>"\n")  # Delete any spaces before a newline.
    s = replace(s, r"\n  *"=>"\n")  # Delete any spaces after a newline.
    s = replace(s, r"  *$"=>"")  # Delete any trailing spaces.
    s = replace(s, r"  *"=>",")  # Replace any remaining run of spaces with a single comma.
    io = IOBuffer(s)
    spl = DelimitedFiles.readdlm(io, ',', Float64, '\n'; skipstart=4, skipblanks=true)

    return freq, spl
end

end # module
