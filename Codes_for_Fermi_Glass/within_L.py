from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

def poisson_fit_func(s, Δ):
    """Poisson-like fitting function."""
    return (np.exp(-s / Δ)) / Δ

def load_and_plot():
    # Load data
    energy_spacing = np.load("S=3.npy")

    # Create histogram
    bins_count = int(np.sqrt(len(energy_spacing))) + 1
    densities, bin_edges, _ = plt.hist(energy_spacing, bins=bins_count, density=True, color="b", alpha=0.7)

    # Use the centers of the bins for curve fitting
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    parameters, _ = curve_fit(poisson_fit_func, bin_centers, densities)
    Δ_fit = parameters[0]

    # Plotting the fitted curve
    x_fit = np.linspace(0, max(energy_spacing), 100)
    y_fit = poisson_fit_func(x_fit, Δ_fit)

    # Setup plot
    plt.plot(x_fit, y_fit, 'r-', label=f"Poisson fit: Δ={Δ_fit:.2f}")
    plt.title(f'Grain size=3, Δ={Δ_fit:.2f}')
    plt.xlabel(f"Energy level spacing / Δ")
    plt.ylabel("Normalized frequency")
    plt.legend()

    # Custom tick formatter
    def custom_x_formatter(x, _):
        return f"{x / Δ_fit:.2f}"

    def custom_y_formatter(y, _):
        return f"{y * Δ_fit:.2f}"

    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FuncFormatter(custom_x_formatter))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(custom_y_formatter))
    
    # Show plot
    plt.show()

load_and_plot()
