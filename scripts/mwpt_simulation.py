import numpy as np
from Basilisk.simulation import spacecraft
from Basilisk.simulation import gravityEffector
from Basilisk.utilities import SimulationBaseClass, macros
from Basilisk.architecture import bskLogging

# # === Setup Simulation ===
scSim = SimulationBaseClass.SimBaseClass()
processName = "simProcess"
simTaskName = "simTask"
scSim.CreateNewProcess(processName)
scSim.CreateNewTask(simTaskName, macros.sec2nano(200))

# # === Create Spacecraft ===
chaser = spacecraft.Spacecraft()
receiver = spacecraft.Spacecraft()

chaser.hub.r_CN_NInit = [[1000.0], [0.0], [0.0]]
receiver.hub.r_CN_NInit = [[-1000.0], [0.0], [0.0]]

scSim.AddModelToTask(simTaskName, chaser)
scSim.AddModelToTask(simTaskName, receiver)

# # === Gravity Effector ===
gravity = gravityEffector.GravityEffector()
scSim.AddModelToTask(simTaskName, gravity)
chaser.gravField.gravBodies = gravity.gravBodies
receiver.gravField.gravBodies = gravity.gravBodies

# # === Recorders ===
chaserRec = chaser.scStateOutMsg.recorder()
receiverRec = receiver.scStateOutMsg.recorder()
scSim.AddModelToTask(simTaskName, chaserRec)
scSim.AddModelToTask(simTaskName, receiverRec)
# # === Set Stop Time Before Initializing ===
scSim.ConfigureStopTime(macros.sec2nano(0.001))  # Set as required
scSim.InitializeSimulation()
scSim.ExecuteSimulation()

# # === Write Logs for Vizard ===
# logger = bskLogging.BSKLogger()
# logger.addMessage(chaser.scStateOutMsg, chaserRec)
# logger.addMessage(receiver.scStateOutMsg, receiverRec)
# logger.save("logs/mwpt_vizard_log.bsk")

# print("Simulation complete. Logs saved for Vizard.")

from Basilisk.utilities import SimulationBaseClass
from Basilisk.utilities import macros


def run():
    """
    Illustration of Basilisk process and task creation
    """

    #  Create a sim module as an empty container
    scSim = SimulationBaseClass.SimBaseClass()

    #  create the simulation process
    dynProcess = scSim.CreateNewProcess("dynamicsProcess")
    fswProcess = scSim.CreateNewProcess("fswProcess")

    # create the dynamics task and specify the integration update time
    dynProcess.addTask(scSim.CreateNewTask("dynamicsTask", macros.sec2nano(5.)))
    dynProcess.addTask(scSim.CreateNewTask("sensorTask", macros.sec2nano(10.)))
    fswProcess.addTask(scSim.CreateNewTask("fswTask", macros.sec2nano(10.)))

    #  initialize Simulation:
    scSim.InitializeSimulation()

    #   configure a simulation stop time and execute the simulation run
    scSim.ConfigureStopTime(macros.sec2nano(20.0))
    scSim.ExecuteSimulation()

    return


if __name__ == "__main__":
    run()

