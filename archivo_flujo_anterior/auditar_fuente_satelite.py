"""Inventario de TIFF originales y contraste con CHIRPS ya descargado.

No reconstruye curvas ni sustituye productos. Solo lee archivos y genera
evidencia_fuente_satelite.json en esta carpeta.
Uso: python auditar_fuente_satelite.py --raw-dir /ruta/ImagenesSatelitePentadal
Dependencias: numpy y tifffile.
"""

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import numpy as np
import tifffile

ROOT = Path(__file__).resolve().parent


def auditar(raw):
    archivos=sorted(p for p in raw.glob('*.tif') if re.fullmatch(r'ppt\d{6}\.tif',p.name))
    if not archivos:
        raise ValueError(f'No se encontraron TIFF pptAAAAPP.tif en {raw}')
    anios=defaultdict(list); metadata=Counter(); ejemplos=[]; fallos=[]
    maximo=-np.inf; minimo=np.inf; iguales_200=0; sin_dato=0
    pares={}
    for p in archivos:
        anio,pentada=int(p.stem[3:7]),int(p.stem[7:9]);anios[anio].append(pentada)
        with tifffile.TiffFile(p) as t:
            a=t.asarray();tags={tag.name:tag.value for tag in t.pages[0].tags.values()}
            nodata=tags.get('GDAL_NODATA')
            valido=np.isfinite(a)
            if nodata is not None:valido &= a != float(nodata)
            sin_dato += int((~valido).sum())
            if valido.any():
                maximo=max(maximo,float(a[valido].max()));minimo=min(minimo,float(a[valido].min()))
                iguales_200 += int(((a==200)&valido).sum())
            for key in ('ImageDescription','Software','DateTime','GDAL_METADATA','GDAL_NODATA'):
                if key in tags:metadata[(key,str(tags[key]))]+=1
            if p in (archivos[0],archivos[-1]) or p.name=='ppt201301.tif':
                scale=tags.get('ModelPixelScaleTag');tie=tags.get('ModelTiepointTag')
                ejemplos.append({'archivo':p.name,'forma':list(a.shape),'dtype':str(a.dtype),
                    'tags':{k:str(v) for k,v in tags.items()},
                    'bounds':([tie[3],tie[4]-a.shape[0]*scale[1],tie[3]+a.shape[1]*scale[0],tie[4]] if scale and tie else None)})
            if anio==2013:pares[pentada]=np.where(valido,a.astype(float),np.nan)
    claves=(ROOT/'datos_web/chirps_v2_valle_1981_2013_claves.txt').read_text().split()
    ch=np.fromfile(ROOT/'datos_web/chirps_v2_valle_1981_2013.bin',dtype='<f4').reshape(len(claves),39,37)
    total2013=0; dif=[]; cap=[]
    for key,a in zip(claves,ch):
        year,month,p=map(int,key.split('.'))
        if year!=2013:continue
        total2013+=1;pentada=(month-1)*6+p
        if pentada in pares:
            m=np.isfinite(a)&(a>-100)&np.isfinite(pares[pentada])
            dif.extend((pares[pentada][m]-a[m]).tolist())
            cap.extend(a[m].tolist())
    resultado={
        'fecha_utc':datetime.now(timezone.utc).isoformat(),'directorio_original_inspeccionado':str(raw.resolve()),
        'total_tiff':len(archivos),'primero':archivos[0].name,'ultimo':archivos[-1].name,
        'por_anio':{str(y):{'n_archivos':len(ps),'indices':sorted(ps),'ausentes_respecto_1_72':sorted(set(range(1,73))-set(ps))} for y,ps in sorted(anios.items())},
        'ejemplos_tiff':ejemplos,
        'metadatos_descriptivos':[{'tag':k,'valor':v,'n_archivos':n} for (k,v),n in metadata.items()],
        'valores_validos':{'minimo':minimo,'maximo':maximo,'conteo_igual_200':iguales_200,'sin_dato':sin_dato},
        'contraste_CHIRPS_local':{'clave_origen':'Notebook 02 §3: descarga externa CHIRPS v2; no implica identidad con los TIFF entregados',
            'pentadas_2013':total2013,'pentadas_2013_comparadas':len(pares),'diferencia_absoluta_media':float(np.mean(np.abs(dif))),
            'diferencia_absoluta_maxima':float(np.max(np.abs(dif))), 'mismos_valores':bool(np.all(np.abs(dif)<1e-6))},
        'alcance':['La identidad del producto original no se demuestra por compartir malla o correlacion.',
            'No se confirma la unidad ni la version original si los tags no la declaran.',
            '72 pasos por año es el patron de la entrega y del calendario mensual CHIRPS; no todos los pasos duran cinco dias.',
            '24 pasos de 2013 no forman un año completo en ese calendario. No preguntar al profesor como ajustar: solicitar los archivos que faltan.',
            'La descarga externa CHIRPS contiene 2013 completo pero no se mezcla silenciosamente con el producto original.'],
    }
    destino=ROOT/'evidencia_fuente_satelite.json';destino.write_text(json.dumps(resultado,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:resultado[k] for k in ('total_tiff','primero','ultimo','metadatos_descriptivos','valores_validos','contraste_CHIRPS_local')},indent=2,ensure_ascii=False))
    print('Años con cobertura distinta de 72:',{y:d for y,d in resultado['por_anio'].items() if d['n_archivos']!=72})
    print('Evidencia guardada:',destino)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',type=Path,required=True)
    auditar(parser.parse_args().raw_dir)
