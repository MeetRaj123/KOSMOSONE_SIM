import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from Basilisk.simulation import spacecraft
from Basilisk.utilities import SimulationBaseClass, macros
from Basilisk.utilities import orbitalMotion, simIncludeGravBody, unitTestSupport, vizSupport



# === Configuration ===
# Paths to your GMAT-exported ephemeris files
CHASER_CSV_PATH = "../ephemeris_data/Chaser_ephemeris.csv"
RECEIVER_CSV_PATH = "../ephemeris_data/Receiver_ephemeris.csv"
LOG_OUTPUT_PATH = "../logs/mwpt_results.csv"

# MWPT parameters (example values)
TRANSMIT_POWER = 1000  # Watts
GAIN_TX = 10           # Linear scale (not dB)
GAIN_RX = 10
WAVELENGTH = 0.03      # Meters (for ~10 GHz microwave)

# === Load Data ===
chaser_data = pd.read_csv(CHASER_CSV_PATH)
receiver_data = pd.read_csv(RECEIVER_CSV_PATH)

results = []

# === MWPT Power Calculation Function ===
def mwpt_power(transmit_power, gain_tx, gain_rx, wavelength, distance):
    path_loss = (4 * np.pi * distance / wavelength) ** 2
    received_power = transmit_power * gain_tx * gain_rx / path_loss
    return received_power

# === Process Data ===
for i in range(len(chaser_data)):
    chaser_pos = np.array([chaser_data['MWPT.EarthMJ2000Eq.X'][i], chaser_data['MWPT.EarthMJ2000Eq.Y'][i], chaser_data['MWPT.EarthMJ2000Eq.Z'][i]])
    receiver_pos = np.array([receiver_data['target.EarthMJ2000Eq.X'][i], receiver_data['target.EarthMJ2000Eq.Y'][i], receiver_data['target.EarthMJ2000Eq.Z'][i]])
    distance = np.linalg.norm(chaser_pos - receiver_pos)
    received_power = mwpt_power(TRANSMIT_POWER, GAIN_TX, GAIN_RX, WAVELENGTH, distance)

    results.append([i, distance, received_power])

# === Save Results ===
with open(LOG_OUTPUT_PATH, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Time Step", "Distance (m)", "Received Power (W)"])
    writer.writerows(results)

print(f"MWPT results saved to {LOG_OUTPUT_PATH}")

# === Plot Results ===
results_df = pd.DataFrame(results, columns=["Time Step", "Distance (m)", "Received Power (W)"])

plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(results_df['Time Step'], results_df['Distance (m)'])
plt.title("MWPT Simulation Results")
plt.ylabel("Distance (m)")
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(results_df['Time Step'], results_df['Received Power (W)'])
plt.xlabel("Time Step")
plt.ylabel("Received Power (W)")
plt.grid(True)

plt.tight_layout()
plt.show()


def run():
    simTaskName = "simTask"
    simProcessName = "simProcess"

    scSim = SimulationBaseClass.SimBaseClass()
    simulationTime = macros.min2nano(500.)  # 5 minutes simulated time converted to nanoseconds

    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(0.5)  # 0.5 second step
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))

    # First spacecraft
    scObject = spacecraft.Spacecraft()
    scObject.ModelTag = "Satellite1"
    I = [900., 0., 0.,
         0., 800., 0.,
         0., 0., 600.]
    scObject.hub.mHub = 750.0 #Sets the mass of the spacecraft’s main hub (center body) to 750 kg
    scObject.hub.IHubPntBc_B = unitTestSupport.np2EigenMatrix3d(I) #converts to eigen matrix

    # Second spacecraft
    scObject2 = spacecraft.Spacecraft()
    scObject2.ModelTag = "Satellite2"
    I2 = [700., 0., 0.,
          0., 700., 0.,
          0., 0., 500.]
    scObject2.hub.mHub = 500.0
    scObject2.hub.IHubPntBc_B = unitTestSupport.np2EigenMatrix3d(I2)

    scSim.AddModelToTask(simTaskName, scObject)
    scSim.AddModelToTask(simTaskName, scObject2)

    gravFactory = simIncludeGravBody.gravBodyFactory()#adding gravity factory object to manage the gravity model lllike earth and moon
    earth = gravFactory.createEarth() # create earth with gravity
    earth.isCentralBody = True #making earth the central body
    mu = earth.mu #value of GM

    gravFactory.addBodiesTo(scObject)#attach the satellite
    gravFactory.addBodiesTo(scObject2)#attach the satellite

    # Orbit elements for Satellite1
    #orbitalMotion python module to work with satellite
    oe = orbitalMotion.ClassicElements() # create an orbital element with 6 keplerian parameters
    oe.a = 7000000.0
    oe.e = 0.0
    oe.i = 33.3 * macros.D2R
    oe.Omega = 48.2 * macros.D2R
    oe.omega = 90.0 * macros.D2R
    oe.f = 0.0 * macros.D2R
    rN, vN = orbitalMotion.elem2rv(mu, oe)#will give position vector and velocity vector
    scObject.hub.r_CN_NInit = rN#postion of first satellite
    scObject.hub.v_CN_NInit = vN #velocity of first satellite
    scObject.hub.sigma_BNInit = [[0.0], [0.0], [0.0]] #Initial attitude (orientation) of the spacecraft body frame (B) relative to the inertial frame (N)
    scObject.hub.omega_BN_BInit = [[0.0], [0.0], [0.0]] #Initial angular velocity of the spacecraft body relative to the inertial frame.

    # Orbit elements for Satellite2 (slightly different)
    oe2 = orbitalMotion.ClassicElements()
    oe2.a = 9500000.0
    oe2.e = 0.0001
    oe2.i = 33.3 * macros.D2R
    oe2.Omega = 50.0 * macros.D2R
    oe2.omega = 92.0 * macros.D2R
    oe2.f = 5.0 * macros.D2R
    rN2, vN2 = orbitalMotion.elem2rv(mu, oe2)
    scObject2.hub.r_CN_NInit = rN2
    scObject2.hub.v_CN_NInit = vN2
    scObject2.hub.sigma_BNInit = [[0.0], [0.0], [0.0]]
    scObject2.hub.omega_BN_BInit = [[0.0], [0.0], [0.0]]

    if vizSupport.vizFound:
        viz = vizSupport.enableUnityVisualization(scSim, simTaskName, [scObject, scObject2], saveFile="two_satellites")
        viz.settings.orbitLinesOn = 1
        viz.settings.trueTrajectoryLinesOn = 1

    scSim.InitializeSimulation()
    scSim.ConfigureStopTime(simulationTime)
    print("Starting two satellites simulation...")
    scSim.ExecuteSimulation()
    print("Simulation completed.")

if __name__ == "__main__":
    run()
