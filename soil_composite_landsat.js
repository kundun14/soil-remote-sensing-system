// POLYGON 
var BUFFER = ee.FeatureCollection('projects/ee-sysiflorida/assets/buffer');

Map.addLayer(BUFFER, {color: 'FF0000'}, 'Buffer Shapefile');
Map.centerObject(BUFFER, 10); 



// #TILES
// # path 4 y row 69
// # 3 69
// # 4 70 
// # 3 70 

// FILTROS

var colFilter = ee.Filter.and( 
  // ee.Filter.eq('WRS_PATH', 4), 
  // ee.Filter.eq('WRS_ROW', 69),
    ee.Filter.lt('CLOUD_COVER', 30), // cobertura de nubes menor a 20%
    ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10), // GEOMETRIC_RMSE_MODEL menor a 10 , validacion respecto a GCP
    ee.Filter.or(
        ee.Filter.eq('IMAGE_QUALITY', 9), // 9 = Best
        ee.Filter.eq('IMAGE_QUALITY_OLI', 9)) // 9 = Best
    ); 


// Scale y offset (se aplica para todas las imagenes de la coleccion 2)


function applyScaleFactors(image) {
  var opticalBands = image.select("SR_B.").multiply(0.0000275).add(-0.2); // solo bandas opticas;
  return image.addBands(opticalBands, null, true); //adjuntar las bandas scaladas a la lista original
}

// Renombrar las bandas de OLI (L8)

function renameOli(img) {
    return img.select(
    ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL"],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', "QA_PIXEL"]);
}

// Renombrar las bandas del ETM+.(L4,5,7)

function renameEtm(img) {
    return img.select(
    ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7", "QA_PIXEL"],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', "QA_PIXEL"]);
}


// TRANSFORMACION ESPACIO ESPECTRAL ETM -> OLI
// Coeficientes de Roy et al. (2016)

var coefficients = {
    itcps: ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]),
    slopes: ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
};

function etmToOli(img) {
    return img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
        .multiply(coefficients.slopes)
        .add(coefficients.itcps)
        .addBands(img.select("QA_PIXEL"));
}

// MASCARA CON QA_PIXEL

// L8 tiene diferentes valores de bit en el QA_PIXEL que L457?


// ver documentacion para valores del QA_PIXEL
// https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/LSDS-1619_Landsat8-C2-L2-ScienceProductGuide-v3.pdf

function bitwiseExtract(value, fromBit, toBit) {
    if (toBit === undefined) toBit = fromBit;
        var maskSize = ee.Number(1).add(toBit).subtract(fromBit);
        var mask = ee.Number(1).leftShift(maskSize).subtract(1);
    return value.rightShift(fromBit).bitwiseAnd(mask);
}

function fmask(img) {
    var qa = img.select("QA_PIXEL");
    var cloudShadowBitMask = 1 << 4;
    var dcloudsBitMask = 1 << 1;
    var cirrusBitMask = 1 << 2;
    var cloudsBitMask = 1 << 3;
    var waterBitMask = 1 << 7;
    var snowBitMask = 1 << 5;
    var dilaBitMask = 0 << 6;
    var cloudConBitMask = bitwiseExtract(qa, 8, 9); //Cloud Confidence
    var cloudShadowConBitMask = bitwiseExtract(qa, 10, 11); //Cloud Shadow Confidence
    var snowConBitMask = bitwiseExtract(qa, 12, 13); // Snow/Ice Confidence
    var mask = qa.bitwiseAnd(cloudShadowBitMask)
            .eq(0)
            .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
            .and(qa.bitwiseAnd(waterBitMask).eq(0))
            .and(qa.bitwiseAnd(snowBitMask).eq(0))
            .and(qa.bitwiseAnd(cirrusBitMask).eq(0))
            .and(qa.bitwiseAnd(dcloudsBitMask).eq(0))
            .and(qa.bitwiseAnd(dilaBitMask).eq(0));
                // .and(qa.bitwiseAnd(cloudConBitMask).eq(3)) // si el Cloud Confidence es alto el valor del bit es 3
                // .and(qa.bitwiseAnd(cloudShadowConBitMask).eq(3)) // si el Cloud Shadow Confidence es alto el valor del bit es 3
                // .and(qa.bitwiseAnd(snowConBitMask).eq(3)); // si el Snow/Ice Confidencee es alto el valor del bit es 3
    return img.updateMask(mask);
}


// Pre OLI

function prepOli(img) {
    var orig = img;
    img = applyScaleFactors(img);
    img = renameOli(img);
    // img = filtro_reflec_valid(img); //  filtra valores de pixel en el rango valido de reflectancia
    return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

//  Pre ETM+ .

function prepEtm(img) {
    var orig = img;
    img = applyScaleFactors(img);
    img = renameEtm(img);
    // img = filtro_reflec_valid(img); //  filtra valores de pixel en el rango valido de reflectancia
    img = etmToOli(img);
    return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

// RANGO VALIDO REFLECTANCIA


// INDICES PAPERS

// ADD INDICES

// function add_vegind(image) {
// var ndvi = image.normalizedDifference(["NIR", "Red"]).rename("ndvi");
// var nbr2 = image.normalizedDifference(["SWIR1","SWIR2" ]).rename("nbr2");
// return image.addBands(ndvi).addBands(nbr2);
// }


// ADD INDICES

var add_indices = function(image) {
  // Calculate NDVI
    var ndvi = ee.Image(0).expression(
    "(NIR - Red) / (NIR + Red)", {
        "NIR": image.select("NIR"),
        "Red": image.select("Red")
    }).rename("NDVI");

  // Calculate NBR2
    var nbr2 = ee.Image(0).expression(
    "(SWIR1 - SWIR2) / (SWIR1 + SWIR2)", {
        "SWIR1": image.select("SWIR1"),
        "SWIR2": image.select("SWIR2")
    }).rename("NBR2");

  // Calculate BSI
    var bsi = ee.Image(0).expression(
    "((SWIR2 + Red) - (NIR + Blue)) / ((SWIR2 + Red) + (NIR + Blue))", {
        "SWIR2": image.select("SWIR2"),
        "Red": image.select("Red"),
        "NIR": image.select("NIR"),
        "Blue": image.select("Blue")
    }).rename("BSI");

  // Calculate BSI1
    ar bsi1 = ee.Image(0).expression(
    "((SWIR1 + Red) - (NIR + Blue)) / ((SWIR1 + Red) + (NIR + Blue))", {
        "SWIR1": image.select("SWIR1"),
        "Red": image.select("Red"),
        "NIR": image.select("NIR"),
        "Blue": image.select("Blue")
    }).rename("BSI1");

  // Calculate BSI2
    var bsi2 = ee.Image(0).expression(
    "100 * sqrt((SWIR2 - Green) / (SWIR2 + Green))", {
        "SWIR2": image.select("SWIR2"),
        "Green": image.select("Green")
    }).rename("BSI2");

  // Calculate BSI3
    var bsi3 = ee.Image(0).expression(
    "(((SWIR2 + Red) - (NIR + Blue)) / ((SWIR2 + Red) + (NIR + Blue))) * 100 + 100", {
        "SWIR2": image.select("SWIR2"),
        "Red": image.select("Red"),
        "NIR": image.select("NIR"),
        "Blue": image.select("Blue")
    }).rename("BSI3");

  // Calculate NDSI1
    var ndsi1 = ee.Image(0).expression(
    "(SWIR1 - NIR) / (SWIR1 + NIR)", {
        "SWIR1": image.select("SWIR1"),
        "NIR": image.select("NIR")
    }).rename("NDSI1");

  // Calculate NDSI2
    var ndsi2 = ee.Image(0).expression(
    "(SWIR2 - Green) / (SWIR2 + Green)", {
        "SWIR2": image.select("SWIR2"),
        "Green": image.select("Green")
    }).rename("NDSI2");

  // Calculate BI
    var bi = ee.Image(0).expression(
    "(SWIR1 - SWIR2) / (SWIR1 + SWIR2)", {
        "SWIR1": image.select("SWIR1"),
        "SWIR2": image.select("SWIR2")
    }).rename("BI");

  // Calculate MBI
    var mbi = ee.Image(0).expression(
    "((SWIR1 - SWIR2 - NIR) / (SWIR1 + SWIR2 + NIR)) + 0.5", {
        "SWIR1": image.select("SWIR1"),
        "SWIR2": image.select("SWIR2"),
        "NIR": image.select("NIR")
    }).rename("MBI");

    return image
    .addBands([ndvi, nbr2, bsi, bsi1, bsi2, bsi3, ndsi1, ndsi2, bi, mbi]);
};




// DATA

var DN_L4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2");
var DN_L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2");
var DN_L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2");
var DN_L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");



// DATE FILTER

// PARAMETERS

var MONTHSTART = 5;
var MONTHEND = 9; 


var dateRangeBoundFilter = ee.Filter.and(
    ee.Filter.bounds(BUFFER), ee.Filter.calendarRange(MONTHSTART, MONTHEND, 'month'));

var DN_L4 = DN_L4.filter(dateRangeBoundFilter);
var DN_L5 = DN_L5.filter(dateRangeBoundFilter);
var DN_L7 = DN_L7.filter(dateRangeBoundFilter);
var DN_L8 = DN_L8.filter(dateRangeBoundFilter);


// print(DN_L5);


var sr_L4 = DN_L4.filter(colFilter).map(fmask).map(prepEtm).map(add_indices);
var sr_L5 = DN_L5.filter(colFilter).map(fmask).map(prepEtm).map(add_indices);
var sr_L7 = DN_L7.filter(colFilter).map(fmask).map(prepEtm).map(add_indices);
var sr_L8 = DN_L8.filter(colFilter).map(fmask).map(prepOli).map(add_indices);


// print(sr_L5);

//TESS

// Unir las coleccion L4-5-7-8 

var tess = sr_L4.merge(sr_L5).merge(sr_L7).merge(sr_L8); // 

// PRINT TESS

// print("TESS (L4+L5+L7+L8)", tess);  // propiedades del TESS

// TESS stack para exportar

var selectedBands = ['Red','Green','Blue', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'NBR2', 'BSI', 'BSI1', 'BSI2', 'BSI3', 'NDSI1', 'NDSI2', 'BI', 'MBI'];

var tess = tess.select(selectedBands);

// homogeneizar el tipo de dato de cada banda como float

var tess = tess.map(function(img){
var bandas = img.select(selectedBands);
return img.addBands(bandas, null,  true).toFloat();
});

// // SAVE ALL COLECTION

print(tess);

// EPSG:32718


var SCALE = 30 ;
var batch = require('users/fitoprincipe/geetools:batch');

batch.Download.ImageCollection.toDrive(tess, 'ls-florida', {
    name: 'ls_FiltFlorida{system:index}',
    type: 'float',
    scale: SCALE,
    region: BUFFER
});

