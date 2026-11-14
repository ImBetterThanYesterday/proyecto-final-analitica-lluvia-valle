# EDA inicial del proyecto de precipitación

Este inventario describe estructura y calidad básica. No elimina ni corrige datos.
Los valores atípicos son alertas que requieren revisión con metadatos.

## Resumen de fuentes

| source | file | rows_raw | rows_valid | stations | year_min | year_max | value_min | value_median | value_mean | value_p99 | value_p999 | value_max | zero_count | zero_percent | negative_count | over_500_count | duplicate_key_rows | invalid_calendar_dates | missing_required_values | data_tiff_count | years_with_72_tiffs | years_with_other_tiff_count | tiff_count_2013 | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| estaciones_1980_1990 | P_1980_1990.csv | 997791.0 | 997791.0 | 302.0 | 1980 | 1990 | 0.0 | 0.0 | 6.611451496355449 | 60.0 | 126.0 | 859375.0 | 578166.0 | 57.94459962056182 | 0.0 | 11.0 | 0.0 | 0.0 | 0.0 | nan | nan | nan | nan | nan |
| estaciones_1991_2003 | P_1991_2003.csv | 1010677.0 | 1010677.0 | 284.0 | 1991 | 2003 | 0.0 | 0.0 | 5.469389626953023 | 60.0 | 129.0 | 600.0 | 597017.0 | 59.07099894427201 | 0.0 | 3.0 | 0.0 | 0.0 | 0.0 | nan | nan | nan | nan | nan |
| estaciones_1930_1979 | P_1930_1979.xlsx | 984226.0 | 984226.0 | 264.0 | 1930 | 1979 | 0.0 | 0.0 | 5.427457758685505 | 55.0 | 104.0 | 736.0 | 559657.0 | 56.86265146419623 | 0.0 | 14.0 | 0.0 | 0.0 | 0.0 | nan | nan | nan | nan | nan |
| satelite_pentadal | ImagenesSatelitePentadal | nan | nan | nan | 1981 | 2013 | nan | nan | nan | nan | nan | nan | nan | nan | nan | nan | nan | nan | nan | 2328.0 | 32.0 | 1.0 | 24.0 | El script cuenta archivos; no interpreta valores raster ni sus unidades. |

## Archivos generados

- `resumen_fuentes.csv`: métricas generales por archivo.
- `cobertura_estaciones_todos_los_archivos.csv`: cobertura por estación y año.
- `conteo_tiff_por_anio.csv`: cantidad de imágenes satelitales por año.
- `alertas_mayores_500_*.csv`: observaciones superiores a 500, para inspección.

Nota: el número de días esperados es calendario; no significa que una estación deba tener todos esos registros.