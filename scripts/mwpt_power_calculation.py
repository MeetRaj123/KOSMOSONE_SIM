import numpy as np
import matplotlib.pyplot as plt
from Basilisk.architecture import bskDataStream

chaserStream = bskDataStream.BSKDataStream()
chaserStream.open("logs/chaser_log.bsk")

# Example: read one message
data = chaserStream.readNextRecord()

print(data)


# === Load Logs ===
chaserRec = bskLogging.readLog("logs/chaser_log.bsk")
receiverRec = bskLogging.readLog("logs/receiver_log.bsk")

TRANSMIT_POWER = 1000
GAIN_TX = 10
GAIN_RX = 10
WAVELENGTH = 0.03

def mwpt_power(transmit_power, gain_tx, gain_rx, wavelength, distance):
    path_loss = (4 * np.pi * distance / wavelength) ** 2
    received_power = transmit_power * gain_tx * gain_rx / path_loss
    return received_power

results = []

for i in range(len(chaserRec["r_BN_N"])):
    chaser_pos = np.array(chaserRec["r_BN_N"][i])
    receiver_pos = np.array(receiverRec["r_BN_N"][i])
    distance = np.linalg.norm(chaser_pos - receiver_pos)
    power = mwpt_power(TRANSMIT_POWER, GAIN_TX, GAIN_RX, WAVELENGTH, distance)
    results.append([i, distance, power])

# === Plot ===
results = np.array(results)

plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(results[:, 0], results[:, 1])
plt.title("MWPT Distance and Power (From Recorder Logs)")
plt.ylabel("Distance (m)")
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(results[:, 0], results[:, 2])
plt.xlabel("Time Step")
plt.ylabel("Received Power (W)")
plt.grid(True)

plt.tight_layout()
plt.show()
