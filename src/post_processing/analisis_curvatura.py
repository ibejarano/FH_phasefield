import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
plt.rcParams.update({'font.size': 14})

def analisis_whewell_suavizada(case, window_size=500):
    plt.figure(figsize=(9, 6))

    df = pd.read_csv(case + 'output_coords.csv')

    # Sort the points by the 'x' coordinate to ensure a continuous path
    df = df.sort_values(by='x').reset_index(drop=True)

    # Calculate the fracture length
    lengths = np.sqrt(np.diff(df['x'])**2 + np.diff(df['y'])**2)
    cumulative_length = np.cumsum(lengths)

    # Calculate the angle of each segment
    angles = np.arctan2(np.diff(df['y']), np.diff(df['x'])) * 180 / np.pi

    # Create a DataFrame for angles
    angles_df = pd.DataFrame({
        'distance': cumulative_length,
        'angle': angles
    })

    # Smooth the angle data using a moving average
    angles_df['smoothed_angle'] = angles_df['angle'].rolling(window=window_size).mean()

    # Plot the original and smoothed angle data
    plt.plot(angles_df['distance'], angles_df['angle'], label='Curva original', alpha=0.5)
    plt.plot(angles_df['distance'], angles_df['smoothed_angle'], label=f'Curva suavizada (np={window_size})', color='red', linewidth=2)
    plt.xlabel('Distancia a lo largo de la fractura')
    plt.ylabel('Angulo del segmento (grados)')
    plt.legend()
    plt.grid(True)
    plt.show()

def analisis_whewell(cases=[], labels=[], window_size=500):
    plt.figure(figsize=(9, 6))

    for case_dir, label in zip(cases, labels):
        df = pd.read_csv(case_dir + 'output_coords.csv')

        # Sort the points by the 'x' coordinate to ensure a continuous path
        #df = df.sort_values(by='x').reset_index(drop=True)

        # Calculate the fracture length
        lengths = np.sqrt(np.diff(df['x'])**2 + np.diff(df['y'])**2)
        cumulative_length = np.cumsum(lengths)

        # Calculate the angle of each segment
        angles = np.arctan2(np.diff(df['y']), np.diff(df['x'])) * 180 / np.pi

        # Create a DataFrame for angles
        angles_df = pd.DataFrame({
            'distance': cumulative_length,
            'angle': angles
        })

        # Smooth the angle data using a moving average
        angles_df['smoothed_angle'] = angles_df['angle'].rolling(window=window_size).mean()

        # Plot the original and smoothed angle data
        plt.plot(angles_df['distance'], angles_df['smoothed_angle'], label=label, linewidth=2)

    plt.xlabel('Distancia a lo largo de la fractura')
    plt.ylabel('Angulo del segmento (grados)')
    plt.legend()
    plt.grid(True)
    plt.show()

def analisis_cesaro(cases=[], labels=[], window_size=500):
    plt.figure(figsize=(9, 6))
    # Ecuacion de Cesaro

    for case_dir, label in zip(cases, labels):
        df = pd.read_csv(case_dir + 'output_coords.csv')

        # Sort the points by the 'x' coordinate to ensure a continuous path
        #df = df.sort_values(by='x').reset_index(drop=True)

        # Calculate the tangential angle (phi) in degrees
        lengths = np.sqrt(np.diff(df['x'])**2 + np.diff(df['y'])**2)

        # Create a DataFrame for easier handling
        s = np.cumsum(lengths)
        phi_deg = np.arctan2(np.diff(df['y']), np.diff(df['x'])) * 180 / np.pi
        data = pd.DataFrame({'s': s, 'phi_deg': phi_deg})

        # Smooth the angle data using a moving average
        data['phi_deg_smoothed'] = data['phi_deg'].rolling(window=window_size, center=True).mean()
        data['phi_rad_smoothed'] = np.deg2rad(data['phi_deg_smoothed'])
        valid_data = data.dropna()
        curvature = np.gradient(valid_data['phi_rad_smoothed'], valid_data['s'])
        valid_data['curvature'] = curvature

        # Plot the curvature
        plt.plot(valid_data['s'], valid_data['curvature'], label=label)
    plt.xlabel('Distancia a lo largo de la fractura (s)')
    plt.ylabel('Curvatura ($\kappa = d\phi/ds$)')
    plt.grid(True)
    plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
    #plt.savefig('fracture_curvature.png')
    plt.show()

def comparar_suavizado_whewell(case, window_sizes=[100, 500, 1000]):
    """
    Compares the effect of different smoothing window sizes on the Whewell plot for a single case.

    Args:
        case (str): Path to the case directory.
        window_sizes (list): A list of integer window sizes to compare for the moving average.
    """
    plt.figure(figsize=(10, 7))

    df = pd.read_csv(case + 'output_coords.csv')
    df = df.sort_values(by='x').reset_index(drop=True)

    # Calculate the fracture length
    lengths = np.sqrt(np.diff(df['x'])**2 + np.diff(df['y'])**2)
    cumulative_length = np.cumsum(lengths)

    # Calculate the angle of each segment
    angles = np.arctan2(np.diff(df['y']), np.diff(df['x'])) * 180 / np.pi

    # Create a DataFrame for angles
    angles_df = pd.DataFrame({
        'distance': cumulative_length,
        'angle': angles
    })

    # Plot the original, unsmoothed data for reference
    plt.plot(angles_df['distance'], angles_df['angle'], label='Curva original', alpha=0.3, color='gray', linestyle='--')

    # Plot smoothed data for each window size
    for window_size in window_sizes:
        smoothed_angle = angles_df['angle'].rolling(window=window_size, center=True).mean()
        plt.plot(angles_df['distance'], smoothed_angle, label=f'Suavizado (ventana={window_size})', linewidth=2.5)

    plt.xlabel('Distancia a lo largo de la fractura')
    plt.ylabel('Angulo del segmento (grados)')
    case_name = os.path.basename(os.path.dirname(case))
    plt.title(f'Comparación de Suavizado para: {case_name}')
    plt.legend()
    plt.grid(True)
    plt.show()

def plotear_coordenadas_suavizadas(case, window_sizes=[100, 500, 1000]):
    """
    Plotea las coordenadas x, y originales y las compara con versiones suavizadas.

    Args:
        case (str): Ruta al directorio del caso.
        window_sizes (list): Lista de tamaños de ventana para el suavizado.
    """
    plt.figure(figsize=(10, 7))

    df = pd.read_csv(case + 'output_coords.csv')
    df = df.sort_values(by='x').reset_index(drop=True)

    # Plotear los puntos originales como una nube de puntos
    plt.scatter(df['x'], df['y'], label='Puntos Originales', s=5, alpha=0.4, color='gray')

    # Plotear las coordenadas suavizadas para cada tamaño de ventana
    for window_size in window_sizes:
        # Aplicar el promedio móvil a las coordenadas X e Y
        x_suavizado = df['x'].rolling(window=window_size, center=True).mean()
        y_suavizado = df['y'].rolling(window=window_size, center=True).mean()
        
        plt.plot(x_suavizado, y_suavizado, label=f'Promediado (np={window_size})', linewidth=2.5)

    plt.xlabel('Coordenada X')
    plt.ylabel('Coordenada Y')
    case_name = os.path.basename(os.path.dirname(case))
    plt.title(f'Comparación de Suavizado de Coordenadas para: {case_name}')
    plt.legend()
    plt.grid(True)
    plt.axis('equal') # Asegura que la escala de los ejes sea la misma para no distorsionar la forma
    plt.show()


job_30 = [147, 148, 149]
job_40 = [152, 151, 150]

algoId = [2, 4, 6]
cases_30 = [f"./results/h_30_algo{algoNum}_job_{job}/" for job, algoNum in zip(job_30, algoId)]
cases_40 = [f"./results/h_40_algo{algoNum}_job_{job}/" for job, algoNum in zip(job_40, algoId)]

labels= ["Mallado 2", "Mallado 4", "Mallado 6"]
analisis_whewell(cases_40, labels, window_size=500)
analisis_cesaro(cases_40, labels, window_size=500)


# comparar_suavizado_whewell(cases[0])

# --- Ejemplo de uso para la nueva función de coordenadas suavizadas ---
# plotear_coordenadas_suavizadas(cases_30[0], window_sizes=[50, 200, 400])