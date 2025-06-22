###
# plot the data
###

import matplotlib.pyplot as plt
import numpy as np

def plot_energy_scan(bond_lengths: list[float], energies: list[float],
                     title: str = "Energy vs. Bond Length Scan",
                     xlabel: str = "Bond Length (Angstroms)",
                     ylabel: str = "Energy (Hartree)") -> str:
    """
    Generates an energy vs. bond length plot and returns it as a base64 encoded PNG string.

    Args:
        bond_lengths: List of bond lengths (x-axis data).
        energies: List of energies (y-axis data).
        title: Title of the plot.
        xlabel: Label for the x-axis.
        ylabel: Label for the y-axis.

    Returns:
        A base64 encoded PNG image string.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(bond_lengths, energies, marker='o', linestyle='-', color='blue')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.tight_layout()

    # Save plot to an in-memory binary stream
    buffer = io.BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0) # Rewind the buffer to the beginning

    # Encode to Base64
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    plt.close() # Close the plot to free memory

    return img_base64


