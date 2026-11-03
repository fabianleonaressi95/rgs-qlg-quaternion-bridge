import matplotlib.pyplot as plt

# Data extracted from simulation runs
lattices = ['IQC Truss (Aperiodic)', 'Cubic Truss (Periodic)']
energy = [21495.04, 9.38]

plt.figure(figsize=(8, 6))
bars = plt.bar(lattices, energy, color=['#4C72B0', '#55A868'])
plt.ylabel('Final System Energy (J)')
plt.title('Energy Dissipation Comparison: IQC vs. Periodic Lattice')
plt.yscale('log') # Use log scale due to the large magnitude difference

# Add labels on bars
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.2f} J', va='bottom', ha='center')

plt.tight_layout()
plt.savefig('/workspace/energy_dissipation_plot.png')
print("Plot saved to /workspace/energy_dissipation_plot.png")
