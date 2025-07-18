import math
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# The path to the location of Basilisk
# Used to get the location of supporting data.
from Basilisk import __path__
from Basilisk.simulation import spacecraft
from Basilisk.utilities import SimulationBaseClass
from Basilisk.utilities import macros
from Basilisk.utilities import orbitalMotion
from Basilisk.utilities import simIncludeGravBody
from Basilisk.utilities import simIncludeThruster
from Basilisk.simulation import thrusterStateEffector
from Basilisk.utilities import unitTestSupport  # general support file with common unit test functions
# attempt to import vizard
from Basilisk.utilities import vizSupport
from Basilisk.simulation import thrusterDynamicEffector

CHASER_deltaV_PATH = "../ephemeris_data/chaser_delta_v.csv"
CHASER_deltaV2_PATH = "../ephemeris_data/chaser_delta_v2.csv"
LOG_OUTPUT_PATH = "../logs/mwpt2_results.csv"

firstDeltaV_data = pd.read_csv(CHASER_deltaV_PATH)
secondDeltaV_data = pd.read_csv(CHASER_deltaV2_PATH)

# first_DeltaV = np.array([firstDeltaV_data['DeltaVx'][0], firstDeltaV_data['DeltaVy'][0], firstDeltaV_data['DeltaVz'][0]])
# second_DeltaV = np.array([secondDeltaV_data['DeltaVx'][0], secondDeltaV_data['DeltaVy'][0], secondDeltaV_data['DeltaVz'][0]])

first_DeltaV = np.array([-0.158742232009,0.532,-1.25])
second_DeltaV = np.array([0.251754096160,0,0])
# Variables to track burn status
burn1_done = False
burn2_done = False

bskPath = __path__[0]
fileName = os.path.basename(os.path.splitext(__file__)[0])

def check_and_apply_deltaV(scSim, velRef, posRef, r_trigger, deltaV, burn_done_flag, label):

    TOLERANCE = 100  # [m] tolerance for triggering event

    r_vec = unitTestSupport.EigenVector3d2np(posRef.getState())
    # print(f"here is the initail postion: {posRef.getState()} and ")
    r_mag = np.linalg.norm(r_vec)
    r_mag_km = r_mag/1000
    r_trigger_km = r_trigger/1000
    if not burn_done_flag and abs((r_mag_km - r_trigger_km)/1000) < TOLERANCE:

        v_vec = unitTestSupport.EigenVector3d2np(velRef.getState())
        print(f" current distance : {r_mag_km} \n and {label} : {r_trigger_km}")
        print(f"difference between burning position and {label}:  {r_mag_km - r_trigger_km}")
    
        print(v_vec)
        v_vec += (deltaV)*1000  # apply delta-V
        print(v_vec)
        velRef.setState(v_vec)
        print(f"{label} burn applied at r = {r_mag/1000:.3f} km")
        return True
    return burn_done_flag


def run(show_plots, maneuverCase):
    # Create simulation variable names
    simTaskName = "simTask"
    simProcessName = "simProcess"

    #  Create a sim module as an empty container
    scSim = SimulationBaseClass.SimBaseClass()

    #
    #  create the simulation process
    #
    dynProcess = scSim.CreateNewProcess(simProcessName)

    # create the dynamics task and specify the integration update time
    simulationTimeStep = macros.sec2nano(10.)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))

    #
    #   setup the simulation tasks/objects
    #

    # initialize spacecraft object and set properties
    scObject = spacecraft.Spacecraft()
    scObject.ModelTag = "target"
    
    scObject2 = spacecraft.Spacecraft()
    scObject2.ModelTag = "chaser"

    # add spacecraft object to the simulation process
    scSim.AddModelToTask(simTaskName, scObject)
    scSim.AddModelToTask(simTaskName, scObject2)

    # setup Gravity Body and SPICE definitions
    gravFactory = simIncludeGravBody.gravBodyFactory()
    earth = gravFactory.createEarth()
    earth.isCentralBody = True  # ensure this is the central gravitational body

    # attach gravity model to spacecraft
    gravFactory.addBodiesTo(scObject)
    gravFactory.addBodiesTo(scObject2)

    # setup spice library for Earth ephemeris
    timeInitString = "2024 September 21, 21:00:00.0 TDB"
    spiceObject = gravFactory.createSpiceInterface(time=timeInitString, epochInMsg=True)
    spiceObject.zeroBase = 'Earth'

    # need spice to run before spacecraft module as it provides the spacecraft translational states
    scSim.AddModelToTask(simTaskName, spiceObject)

    #
    #   setup orbit and simulation time
    #
    # setup the orbit using classical orbit elements
    oe = orbitalMotion.ClassicElements()
    rLEO = 9000000  # meters
    rGEO = 12000000 
    oe.a = rLEO
    oe.e = 0.01
    r_appogee = (oe.a)*(1+(oe.e))
    r_perigee = (oe.a)*(1-(oe.e))
    oe.i = 25.0 * macros.D2R
    oe.Omega = 0 * macros.D2R
    altitude = (oe.a)*(1 - (oe.e)) - 6378000
    print(f"chaser's altitude is { altitude/1000} km\n",)

    oe.omega = 0 * macros.D2R
    oe.f = 0 * macros.D2R
    rN, vN = orbitalMotion.elem2rv(earth.mu, oe)
    scObject.hub.r_CN_NInit = rN  # m - r_CN_N
    scObject.hub.v_CN_NInit = vN  # m - v_CN_N

    oe2 = orbitalMotion.ClassicElements()
    rLEO2 = 7000000  # meters
    # rGEO1 = math.pow(earth.mu / math.pow((2. * np.pi) / (24. * 3600.), 2), 1. / 3.)
    oe2.a = rLEO2
    oe2.e = 0.01
    altitude2 = (oe2.a)*(1 - (oe2.e)) - 6378000
    r_appogee2 = (oe2.a)*(1+(oe2.e))
    r_perigee2=(oe2.a)*(1-(oe2.e))
    print(f"target's altitude is {altitude2/1000} km \n")
    print(altitude2)
    oe2.i = 25 * macros.D2R
    oe2.Omega = 0* macros.D2R
    oe2.omega = 0 * macros.D2R
    oe2.f =0 * macros.D2R
    rN2, vN2 = orbitalMotion.elem2rv(earth.mu, oe2)
    scObject2.hub.r_CN_NInit = rN2  # m - r_CN_N
    scObject2.hub.v_CN_NInit = vN2  # m - v_CN_N

    # set the simulation time
    n = np.sqrt(earth.mu / oe.a / oe.a / oe.a)
    P = 2. * np.pi / n
    simulationTime = macros.sec2nano(10* P)



    #
    #   Setup data logging before the simulation is initialized
    #
    numDataPoints = 100
    samplingTime = unitTestSupport.samplingTime(simulationTime, simulationTimeStep, numDataPoints)
    dataRec = scObject.scStateOutMsg.recorder(samplingTime)
    dataRec1 = scObject2.scStateOutMsg.recorder(samplingTime)
    scSim.AddModelToTask(simTaskName, dataRec)
    scSim.AddModelToTask(simTaskName, dataRec1)
    if vizSupport.vizFound:
        # if this scenario is to interface with the BSK Viz, uncomment the following lines
        viz = vizSupport.enableUnityVisualization(scSim, simTaskName, [scObject, scObject2],
                                                    oscOrbitColorList=[vizSupport.toRGBA255("red"),vizSupport.toRGBA255("yellow")],
                                                    trueOrbitColorList=[vizSupport.toRGBA255("green"),vizSupport.toRGBA255("blue")],
                                                    saveFile="target_ahead")
        viz.settings.mainCameraTarget = "earth"
        viz.settings.showCelestialBodyLabels = 1
        viz.settings.showSpacecraftLabels = 1
        viz.settings.truePathRelativeBody = "earth"
        viz.settings.trueTrajectoryLinesOn = 3  # relative to celestial body inertial frame
    scSim.InitializeSimulation()
    posRef = scObject.dynManager.getStateObject(scObject.hub.nameOfHubPosition)
    velRef = scObject.dynManager.getStateObject(scObject.hub.nameOfHubVelocity)
    print(velRef.getState())
# Variables to track burn status
    burn1_done = False
    burn2_done = False

    n_check_steps = 300  # number of total checks in simulation
    check_time_step = macros.sec2nano(10.)  # Check every 10 sec
    
    # Calculate period in seconds
    P1 = 2 * np.pi * np.sqrt(oe2.a**3 / earth.mu)

# Convert to nanoseconds if using Basilisk time units
    sim_time_ns = macros.sec2nano(P1)
    scSim.ConfigureStopTime(sim_time_ns)
    scSim.ExecuteSimulation()





# Total simulation time should be long enough to reach apogee 
    for i in range(n_check_steps):
    # Run sim for the next time step (eventually, this will reach perigee then apogee)
        scSim.ConfigureStopTime(scSim.TotalSim.CurrentNanos + check_time_step)
        scSim.ExecuteSimulation()
        burn1_done = check_and_apply_deltaV(scSim, velRef, posRef, r_appogee, first_DeltaV, burn1_done, "Apogee")
        
    
    # Check apogee
    print(burn1_done)
    
    print(burn1_done)
    # Check perigee
    if burn1_done and not burn2_done:
        print(burn2_done)
        burn2_done = check_and_apply_deltaV(scSim, velRef, posRef, r_perigee2, second_DeltaV, burn2_done, "perigee")
        print(burn2_done)


    T2 = macros.sec2nano(P * 0.25)
    scSim.ConfigureStopTime(simulationTime + T2)
    scSim.ExecuteSimulation()
    print(burn2_done)
    print(f"Position : {posRef.getState()} \n Velocity : {velRef.getState()}")

    gravFactory.unloadSpiceKernels()

    posData = dataRec.r_BN_N
    plt.figure()
    for idx in range(3):
        plt.plot(dataRec.times() * macros.NANO2HOUR, posData[:, idx] / 1000., label=f"Position {idx}")
    plt.legend()
    plt.xlabel('Time [h]')
    plt.ylabel('Inertial Position [km]')
    if show_plots:
        plt.show()
    plt.close("all")

    # Plot velocity magnitude vs time
    velData = dataRec.v_BN_N  # Velocity components in m/s
    velMag = np.linalg.norm(velData, axis=1)  # Magnitude at each time step
    timeData = dataRec.times() * macros.NANO2SEC  # Convert nanoseconds to seconds

    plt.figure()
    plt.plot(timeData,velData, color='orange', label="Velocity Magnitude")
    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity Magnitude [m/s]')
    plt.title('Chaser Velocity vs Time')
    plt.grid(True)

    if show_plots:
        plt.show()
    plt.close("all")


if __name__ == "__main__":
    run(True, 0)
