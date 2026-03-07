import matplotlib.pyplot as plt
import numpy as np



def plotear_casos():
    fig, axs = plt.subplots(2,1, sharex=True)

    casos = [
    "./tests/output_axi_test/"
    ]

    for caso in casos:
        csv_dir = caso + "results.csv"
        out = np.loadtxt(csv_dir, delimiter=',', skiprows=1)
        t, p, v, wplus, wminus = out[:, 0], out[:, 1], out[:, 2], out[:, 3], out[:, 4]
        axs[0].plot(t, p, label=caso)
        axs[1].plot(t, wplus - wminus, label=caso)

    for ax in axs.flatten():
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend()

    axs[0].set_ylabel("Presión [Pa]")
    axs[1].set_ylabel("Apertura [m]")
    axs[1].set_xlabel("Tiempo [s]")

    tteo = np.linspace(1e-7, 1e-1, 100)
    axs[0].plot(t, 3e3*tteo**(-1/3), "k--", label="KGD")
    axs[0].plot(t, 1.3e1*tteo**(-2/3), "r--", label="Radial")

    axs[0].legend()

    plt.show()

def plot_press():
    fig, axs = plt.subplots(1,3, sharex=False, figsize=(16,9))

    casos = [
    "./output_transition_3_final/",
    "./output_transition_4_final/",
    ]

    for caso in casos:
        csv_dir = caso + "results.csv"
        out = np.loadtxt(csv_dir, delimiter=',', skiprows=1)
        t, p, v, wplus, wminus = out[:, 0], out[:, 1], out[:, 2], out[:, 3], out[:, 4]
        axs[0].plot(t, p, label=caso)
        steps = np.arange(len(t))
        total_steps = steps.shape[0]
        axs[1].plot(steps[2:], t[1:total_steps-1] - t[:total_steps-2], label=caso)
        axs[2].plot(t, wplus - wminus, label=caso)

    axs[0].set_xscale("log")
    axs[0].set_yscale("log")

    axs[1].set_yscale("log")
    
    axs[2].set_xscale("log")
    axs[2].set_yscale("log")
    axs[2].legend()

    axs[0].set_ylabel("Presión [Pa]")
    axs[1].set_ylabel("Delta time [m^3]")
    axs[2].set_ylabel("Apertura [m]")
    axs[2].set_xlabel("Tiempo [s]")

    tteo = np.linspace(1e-7, 1e-1, 100)
    axs[0].plot(tteo, 3e3*tteo**(-1/3), "k--", label="KGD")
    axs[0].plot(tteo, 1.02e4*tteo**(-1/5), "r--", label="Radial")

    axs[0].legend()

    axs[1].set_xlim(0, 700)

    plt.show()



def limpiar_csv(csv_dir):
    # Leer el archivo con los errores
    with open(csv_dir, 'r') as f:
        contenido = f.read()

    # Reemplazar el texto "\n" por un salto de línea real
    contenido_corregido = contenido.replace('\\n', '\n')

    # Guardar el resultado
    with open(csv_dir, 'w') as f:
        f.write(contenido_corregido)

    print("¡Archivo corregido guardado como 'archivo_corregido.csv'!")


plot_press()
