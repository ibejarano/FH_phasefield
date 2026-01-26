# Importar módulos necesarios de ParaView
from paraview.simple import *

# 1. Crear una nueva sesión de ParaView (automático al ejecutar con pvpython/pvbatch)
# Connect() # Solo si se conecta a un servidor remoto, no es necesario para local

# 2. Abrir el archivo de datos (Reemplace 'your_data.h5' con su archivo real)
# Supongamos que su archivo es un VTK XML Structured Grid (.vts) o similar
# Si es un archivo .h5, necesitará el lector de ParaView que lo entienda.
# Si el .h5 lee .xml, ParaView debería manejarlo si los lectores están disponibles.
# Ajuste el lector según su tipo de archivo.

def gemini(filename):
    # Ejemplo con un archivo VTK genérico

    #reader = OpenDataFile('your_data.h5') # Asegúrese de que ParaView tiene un lector para su .h5
    reader = Xdmf3ReaderS(registrationName='output.xdmf', FileName=['/home/ignacio/repos/FH-phasefield/results/slurm_job_76/output.xdmf'])

    # Si su archivo tiene pasos de tiempo, asegúrese de que están cargados
    # En ParaView, los timesteps se cargan automáticamente si el lector los soporta.

    # 3. Crear el filtro "Plot Over Line"
    # Define los puntos inicial y final de tu línea (x1, y1, z1) a (x2, y2, z2)
    plotOverLine = PlotOverLine(Input=reader)
    plotOverLine.Point1 = [-1.0, 0.0, 0.0] # Coordenadas del primer punto de la línea
    plotOverLine.Point2 = [1.0, 0.0, 0.0]  # Coordenadas del segundo punto de la línea
    plotOverLine.Resolution = 100 # Número de puntos a lo largo de la línea

    # 4. Guardar la salida del filtro "Plot Over Line" a un archivo CSV
    # Especifica el nombre del archivo de salida
    output_csv_file = 'plot_over_line_data.csv'

    # Para exportar los datos del filtro PlotOverLine
    SaveData(output_csv_file, proxy=plotOverLine)

    print(f"Datos exportados a: {output_csv_file}")

    # Para animaciones o múltiples pasos de tiempo, necesitaría un bucle
    # Ejemplo para guardar cada paso de tiempo (si su archivo tiene timesteps)
    timesteps = GetAnimationScene().TimeSteps
    print("pasos de tiempo", timesteps)

    # if timesteps:
    #     for i, t in enumerate(timesteps):
    #         print(f"Procesando paso de tiempo: {t}")
    #         GetAnimationScene().AnimationTime = t
    #         output_csv_file_timestep = f'plot_over_line_data_timestep_{i:04d}.csv'
    #         SaveData(output_csv_file_timestep, proxy=plotOverLine)
    #     print("Exportación de todos los pasos de tiempo completada.")
    # else:
    #     print("No se encontraron pasos de tiempo en los datos.")

    # Limpiar (opcional, útil en scripts más complejos)
    # Delete(plotOverLine)
    # Delete(reader)
    # del plotOverLine
    # del reader

def tracer_pv(filename):

    # create a new 'Xdmf3 Reader S'
    reader = Xdmf3ReaderS(registrationName='output.xdmf', FileName=['/home/ignacio/repos/FH-phasefield/results/slurm_job_76/output.xdmf'])

    UpdatePipeline(time=0.009999999999999929, proxy=reader)
    timesteps = GetAnimationScene().TimeSteps
    print("pasos de tiempo", timesteps)
    # create a new 'Plot Over Line'
    plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=reader)

    # Properties modified on plotOverLine1
    plotOverLine1.Point1 = [0.0, 0.0, 0.0]
    plotOverLine1.Point2 = [1.0, 0.0, 0.0]

    UpdatePipeline(time=0.009999999999999929, proxy=plotOverLine1)

    UpdatePipeline(time=0.009999999999999929, proxy=plotOverLine1)

gemini("asd")
tracer_pv("asd")