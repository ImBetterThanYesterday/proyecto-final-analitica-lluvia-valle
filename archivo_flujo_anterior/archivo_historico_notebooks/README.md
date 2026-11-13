# Archivo histórico: preparación y prototipos Python

Estos ocho notebooks fueron inspeccionados antes de moverlos y se conservan
íntegros (celdas, salidas y metadatos). No son el flujo funcional vigente.
El flujo vigente es `../proyecto_funcional.R`; no depende de ejecutar notebooks.
Las rutas relativas originales se interpretan desde la raíz del proyecto,
no desde esta carpeta. No ejecutar automáticamente: contienen descargas y
afirmaciones previas de procedencia/sesgo que la guía posterior limita.

| Notebook | Avance preservado / dependencia |
|---|---|
| `01_EDA_calidad_datos.ipynb` | Depuración preliminar, banderas, selección de pares y cobertura. Umbrales heredados, no validados como receta de clase. |
| `01b_valores_unicos.ipynb` | Inventario por variable, estación y año; no imputa. |
| `02_datos_externos.ipynb` | ONI, catálogos, máscara de tierra CHIRPS, contraste de productos y NASA POWER; preservar para futura red. No establece identidad del producto entregado. |
| `03_construccion_funcional.ipynb` | Preparación de matrices observadas, calendario sin 29-feb, tasas pentadales y coeficientes históricos. Las matrices `.csv.gz` existentes son las entradas del flujo R. |
| `04_eda_funcional.ipynb` | Selección de tierra por claves y prototipo EDA; sus resultados incrustados son históricos. |
| `EDA_completo_faltantes_precipitacion.ipynb` | Auditoría de nulos explícitos/implícitos y fechas, sin eliminación automática. |
| `EDA_precipitacion_valle_cauca.ipynb` | Exploración inicial sin corrección automática. |
| `EDA_precipitacion_valle_cauca copy.ipynb` | Variante con avances adicionales; no se descartó como duplicado. |

No había notebook de entrenamiento de red en la raíz. No se creó ni entrenó
ninguna red. Se conservaron `eda_fase0/`, `datos_web/`, `fda_fase1/` y
`fda_fase2/` sin sobrescribir sus datos; las salidas nuevas van a `resultados_R/`.
Los verificadores Python en la raíz se conservan para auditar el historial;
se ajustó su localizador para encontrar aquí los notebooks.
