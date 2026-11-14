# KAN preservada

Esta carpeta conserva el trabajo propio de la arquitectura KAN y del emparejamiento espacial. Los datos crudos y el flujo EDA/FDA de Jacobo permanecen en la raíz del proyecto.

## Integración con Jacobo

La integración ya utiliza directamente las salidas de Jacobo:

1. leer los pares elegibles desde `../eda_fase0/pares_estacion_anio_validos_preliminar.csv`;
2. reconstruir la columna de partición por año: entrenamiento 1981–1995, validación 1996–2000 y prueba 2001–2003;
3. leer las curvas funcionales desde `../fda_fase1/`;
4. escribir las matrices y resultados del piloto dentro de `KAN/resultados/`, sin modificar las salidas reproducibles de Jacobo;
5. usar los catálogos espaciales centralizados en `../Datos/`.

No se debe seleccionar arquitectura, transformación, número de bases ni hiperparámetros con validación o prueba. La prueba se abre una sola vez al final.

La auditoría actual valida 4.485 curvas terrestres de 248 estaciones. El cruce espacial completo —coordenadas, elevación y pertenencia a la grilla— produce 2.051 pares de 105 estaciones: 1.349 de entrenamiento, 445 de validación y 257 de prueba. Los detalles reproducibles están en `resultados/auditoria_entradas_jacobo.json`.
