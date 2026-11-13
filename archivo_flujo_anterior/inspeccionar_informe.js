// macOS: PDFKit mediante JXA, sin controlar aplicaciones ni instalar paquetes.
// osascript -l JavaScript inspeccionar_informe.js archivo.pdf carpeta_existente
ObjC.import('PDFKit');
ObjC.import('AppKit');
ObjC.import('Foundation');
function run(args) {
    if (args.length !== 2) throw Error('Indique PDF y carpeta de salida existente');
    const pdf = $.PDFDocument.alloc.initWithURL($.NSURL.fileURLWithPath(args[0]));
    const n = Number(pdf.pageCount);
    const paginas = [];
    console.log('Páginas PDF: ' + n);
    for (let i = 0; i < n; i++) {
        const page = pdf.pageAtIndex(i);
        const texto = ObjC.unwrap(page.string) || '';
        const image = page.thumbnailOfSizeForBox($.NSMakeSize(1000, 1414), 0);
        const bitmap = $.NSBitmapImageRep.imageRepWithData(image.TIFFRepresentation);
        const data = bitmap.representationUsingTypeProperties(4, $.NSDictionary.dictionary);
        const nombre = 'pagina-' + ('0' + (i + 1)).slice(-2) + '.png';
        data.writeToFileAtomically(args[1] + '/' + nombre, true);
        paginas.push({pagina: i + 1, caracteres: texto.length, texto: texto});
        console.log('Página ' + (i + 1) + ': ' + texto.length + ' caracteres — ' + texto.split('\n').slice(0, 2).join(' / '));
    }
    $.NSString.alloc.initWithUTF8String(JSON.stringify(paginas, null, 2))
        .writeToFileAtomicallyEncodingError(args[1] + '/inspeccion_paginas.json', true, $.NSUTF8StringEncoding, null);
    return 'PDF renderizado a PNG y texto por página guardado.';
}
