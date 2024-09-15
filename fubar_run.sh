#!/bin/bash

# Directorio base que contiene las carpetas cluster_1, cluster_2, ..., cluster_1244
BASE_DIR="data/phylogeny"

# Número de núcleos de CPU a usar para HYPHY
NUM_CPUS=4

# export HYPHYLIB="/home/vivianbass/Documents/hyphy/res/TemplateBatchFiles"

# Ruta al ejecutable de HYPHY (ajusta si es necesario)
HYPHY_EXEC="/home/vivianbass/Documents/hyphy"

# Iterar sobre cada carpeta
for i in $(seq 180 200); do
    # Construir el nombre de la carpeta y el archivo .bf
    FAM_DIR="${BASE_DIR}/cluster_${i}"
    FUBAR_BF="${FAM_DIR}/cluster_${i}_hyphy_fubar_input.bf"

    # Verificar si el archivo .bf existe
    if [ -f "$FUBAR_BF" ]; then
        echo "Procesando ${FUBAR_BF}..."

        # Ejecutar HYPHYMP para el archivo .bf
        hyphy "${FUBAR_BF}"

        echo "Procesamiento de ${FUBAR_BF} completado."
    else
        echo "Archivo ${FUBAR_BF} no encontrado. Saltando..."
    fi
done

echo "Todos los trabajos de FUBAR han sido procesados."
