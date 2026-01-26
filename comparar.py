import matplotlib.pyplot as plt
import numpy as np



def plotear_casos():
    fig, axs = plt.subplots(2,1, sharex=True)

    casos = [
    "./output_legible_dt2/",
    "./output_legible_dt5/",
    "./output_adaptativo/",
    "./output_adaptativo_2/",
    "./output_adaptativo_3/"
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


plotear_casos()