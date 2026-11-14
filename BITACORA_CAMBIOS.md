# Bitácora de cambios

Qué cambió en cada sesión, en qué commit y por qué. El detalle de cada decisión está en `DECISIONES.md` (D-xx) y el estado frente a la rúbrica en `CONTROL_RUBRICA.md`.

## Sesión del 14-nov-2026 (rama `gus-cierre-funcional`)

Punto de partida: commit `fe34b7d` (parte funcional cerrada en un solo archivo, R4 calculada).

### 1. Correcciones rápidas — commit `46e3852`

| Archivo | Cambio | Por qué |
|---|---|---|
| `lluvia_valle_funcional.R` | "375 curvas" → **370**; título "1996-1999" → **1996-2000** | Eran errores de texto; los conteos salen de la tabla de pares |
| `lluvia_valle_funcional.R` | Línea del "mismo K": ahora dice que el diagnóstico global no cambia, pero el veredicto día a día coincide solo en **281 de 365 días** | Antes exageraba ("las diferencias no salen de usar K distinto") |
| `lluvia_valle_funcional.R` | **Prueba 9 reescrita:** prueba de no fuga real. Se reemplazan por ruido validación y prueba, se reajusta y la media, el FPCA, las escalas y B salen idénticos; con la prueba alterada se vuelven a elegir M y λ iguales | La prueba vieja solo revisaba las máscaras de partición (era tautológica) |
| Los dos `.R` | Fijan **UTF-8** solos al arrancar | En la Mac del grupo LANG está vacío y el script se detenía al comparar "validación" |
| `lluvia_valle_funcional.R` | Orden fijo de los grupos de altitud en las tablas F3 y F4 | El orden dependía del idioma de la máquina |
| `.md` | P-4 resuelta en README y CONTEXTO_EDA; CONTROL_RUBRICA y AUDITORIA sin referencias viejas; P-20 registrada | Estaban desactualizados |

### 2. Datos del curso localizados

`datos_curso/` nunca estuvo en git (está en `.gitignore`). Los archivos están en la rama de Federico, `origin/Arquitectura-KAN-Base`, carpeta `Datos/`. Se sacan sin cambiar de rama:

```bash
mkdir -p datos_curso
git archive origin/Arquitectura-KAN-Base Datos | tar -x --strip-components=1 -C datos_curso
find datos_curso -name '._*' -delete
```

⚠️ Esa rama sube los datos crudos al repositorio y quitó `datos_curso/` del `.gitignore`: **no fusionarla así**.
También se alineó el `main` local, que era un commit suelto: `git branch -f main origin/main`.

### 3. R0 y lectura del satélite — commits `23db96f` y `b1bd344`

| Cambio | Por qué | Resultado |
|---|---|---|
| **Regla R0** en `lluvia_valle_datos.R`: valor < 0 o > 1 825 mm → faltante (D-50) | La rúbrica 4.1 pide una regla explícita para valores imposibles; 1 825 mm es el récord mundial de lluvia en 24 h de la OMM | En el Valle no quita nada (máximo 600 mm); en el crudo completo detecta 859 375 mm de una estación de fuera. Fila "1b. R0" en la trazabilidad y prueba nueva: **27 pruebas** |
| El satélite se lee **por lotes de 200** archivos | Abrir los 2 328 GeoTIFF a la vez superaba el límite de archivos abiertos de macOS | El script corre completo en la Mac |
| **P-20 cerrada: K = 8** (D-51) | Con K = 16 la fecha del pico parecía mejorar, pero el K que gana cambia con la referencia (medias móviles de 15, 31 o 61 días) | Tablas `F2_K_picos_vs_referencia.csv` y `F2_K_que_gana_por_referencia.csv`; el script se detiene si algún día una sola K ganara con todas |

### 4. Red neuronal funcional y ablación — `lluvia_valle_red.py` (D-52)

Archivo nuevo en Python. Lee solo `salidas_R_funcional/para_python/`.

| Paso | Qué se hizo | Cómo se verificó |
|---|---|---|
| 1. Insumos | Lee pares, coeficientes en √P, W, R, FPCA y la especificación de R | W igual a la de R (8e-14); FPCA igual (1e-14); reproduce las métricas de prueba de la regresión de R |
| 2. Red | Codificador z = A W c (8 neuronas) + altitud + ONI; decodificador FPCA de R congelado (M = 6); cabezas lineal, KAN [10,4,6] y MLP [10,37,6] | — |
| 3. Lineal = regresión | Demostración escrita en el código | La red con los coeficientes de R da sus predicciones (3,6e-15); la identidad con β(s,t) a 1,5e-14 (criterio 1e-4) |
| 4. KAN y MLP | KAN probada con una función sintética; MLP con los mismos parámetros | RMSE 2e-3 (criterio 1e-2); 640 vs 635 parámetros (0,8 %) |
| 5. Ablación | Mismo presupuesto: 5 pruebas de λ_w por cabeza y 3 semillas | 14 de 14 pruebas pasan |

**Tres cambios frente al script viejo (`archivo_flujo_anterior/05_red_funcional.py`):**

1. **Rugosidad relativa** ∫w''²/∫w²: con z estandarizado, la penalización vieja se podía bajar encogiendo las w_k sin suavizarlas.
2. **Optimizador:** con lr = 3e-3 y paciencia 300, la cabeza lineal paraba antes de converger (pérdida de validación 142, peor que mínimos cuadrados). Ahora llega a 126, igual que la regresión de R. Sin eso, la ventaja de KAN y MLP habría sido falsa.
3. **Rejilla de λ_w** {0, 1e6, …, 1e9}: la vieja (≤ 0,1) no penalizaba nada.

**Resultado en prueba (2001–2003):** L² relativo satélite 0,513 → lineal 0,440, KAN 0,434, MLP 0,420. El MLP supera a la lineal por poco; la KAN no se distingue de la lineal en L² ni KGE, pero reduce el sesgo. La no linealidad se ve más en validación (que tiene el Niño de 1997 y las Niñas de 1999 y 2000) que en prueba.

### Qué sigue

1. 4.5: graficar las w_k(t) y las activaciones de la KAN, ajuste simbólico y peores casos.
2. 4.4: métricas de las tres cabezas por fase ENSO.
3. Decisiones abiertas: P-16 (ONI anual), P-19 (1930–1979) y confirmar con el profesor el CSV de estaciones reservadas.
4. Unir los dos R en uno para la entrega; informe PDF; anexo de IA con los prompts.

### Cómo correr todo

```bash
Rscript lluvia_valle_datos.R        # ~5 min, 27 pruebas (necesita datos_curso/)
Rscript lluvia_valle_funcional.R    # ~1 min, 12 pruebas
python lluvia_valle_red.py          # ~15 min, 14 pruebas (pip install torch numpy pandas scipy matplotlib)
```
