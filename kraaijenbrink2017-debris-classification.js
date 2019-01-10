////////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                            ////
////   Earth Enginge script to classify the surface of all glaciers in High Mountain Asia to    ////
////   distinguish between debris-covered and debris-free ice, and supraglacial ponds.          ////
////                                                                                            ////
////   Script also generates a map, data inpection panel, and a quick access panel.             ////
////                                                                                            ////
////                                                                                            ////
////   Author:                                                                                  ////
////   Philip Kraaijenbrink (p.d.a.kraaijenbrink@uu.nl)                                         ////
////                                                                                            ////
////   Reference:                                                                               ////
////   Kraaijenbrink PDA, Bierkens MFP, Lutz AF and Immerzeel WW (2017).                        ////
////   Impact of a global temperature rise of 1.5 degrees Celsius on Asia’s glaciers. Nature.   ////
////   http://doi.org/10.1038/nature23878                                                       ////
////                                                                                            ////
////////////////////////////////////////////////////////////////////////////////////////////////////


// Data import =====================================================================================

var srtm30     = ee.Image("USGS/SRTMGL1_003");
var ls8rfl     = ee.ImageCollection("LANDSAT/LC8_L1T_TOA");
var rgihull    = ee.FeatureCollection("users/philipkraaijenbrink/hma-15-degrees/rgi50_asia_10kmBuff_hull");
var rgimask    = ee.Image("users/philipkraaijenbrink/hma-15-degrees/rgi50_arealim025_05arcsec_rgiid");
var rgiregions = ee.FeatureCollection("users/philipkraaijenbrink/hma-15-degrees/rgi50_O2Regions_asia");
var astem      = ee.Image("NASA/ASTER_GED/AG100_003");

    

// Data preparation ================================================================================

// rgi extent
var rgiext = rgihull.geometry().bounds();

// enlarge each 1 arcsec RGI mask with 2 px to ensure full rgi glacier coverage later on
var rgimask = rgimask.reduceNeighborhood({reducer: ee.Reducer.max(),kernel: ee.Kernel.square(2)});

// filter the landsat 8 collection
var ls8rfl = ls8rfl
        .filterBounds(rgihull)                         // only use scenes within the RGI asia coverage
        .filter(ee.Filter.lte('CLOUD_COVER',20));      // only use scenes that have a relatively low cloud cover
        
// make a single composite image for which each pixel is ordered by the warmest in the whole collection
var ordered = ls8rfl.qualityMosaic('B10');

// mask the ordered image using RGI raster mask
var ordered = ordered.updateMask(rgimask);



// TIR prep =====================================================================================================

// get and correct landsat 8 thermal data
var tbrightK   = ordered.select('B10');                    // LS8 brightness temp in celsius

// get emissivity map for the landsat B10 emmisivity range
var emissivity  = astem.select('emissivity_band13','emissivity_band14').reduce(ee.Reducer.mean()).divide(1000).unmask()
var emissivity  = emissivity.where(emissivity.eq(0), 0.96).updateMask(rgimask) // fill masked aster data with fixed e value

// correct to surface temperature by TS=TB/[1+(?*TB/c2)*ln(e)] (Weng et al., 2004)
var planck     = ee.Image(6.626e-34);
var boltzmann  = ee.Image(1.38e-23);
var lightvel   = ee.Image(2.998e8);
var c2         = planck.multiply(lightvel).divide(boltzmann).multiply(1e6);
var labdaB10   = ee.Image(10.8);
var tsurfK     = tbrightK.divide(ee.Image(1).add((labdaB10.multiply(tbrightK).divide(c2)).multiply(emissivity.log())))
var tbrightC   = tbrightK.subtract(273.15);
var tsurfC     = tsurfK.subtract(273.15);



// Calculate indices and derivatives =================================================================

var ndsi      = ordered.normalizedDifference(['B5','B6']);                  // snow/ice index
var ndvi      = ordered.normalizedDifference(['B5','B4']);                  // veg index
var ndbi      = ordered.normalizedDifference(['B4','B1']);                  // blue index
var brnsvis   = ordered.select('B2','B3','B4')                              // brightness in visual light
                       .reduce(ee.Reducer.mean()).rename('vis_brightness');
var srtmslope = ee.Terrain.slope(srtm30);                                   // slope of the 1arcsec SRTM
var hsv       = ordered.select(['B4', 'B3', 'B5']).rgbToHsv();              // hue,sat,val transform
var wsi       = hsv.normalizedDifference(['value','hue']).rename('WSI')     // water resistant snow index (Sharma et al, 2016)
var wateridx = ee.Image(-1).multiply(ndvi);                                 // final water index used



// Image classification ================================================================================

// classification settings
var maxpondslope = 20;    // max slope for supraglacial pond class
var maxdebslope  = 24;    // max slope for debris class
var ndsithresh   = 0.25;
var watidxthresh = 0.1;
var tempthreshdeb = 10;
var tempthreshpond = 5;

// debris classifation
var debris = ndsi.lt(ndsithresh)

// ice classification
var ice = debris.eq(0);

// get small patches of ice
var icesize = ice.updateMask(ice).connectedPixelCount(50).reproject({crs: ndsi.projection(), scale: 30});
var smallpatches = icesize.lt(50).unmask().updateMask(rgimask);

// water classification
// var water = wsi.lt(-0.6).and(debris);
var water = wateridx.gt(watidxthresh)
                    .and(debris)
                    .and(tsurfC.gt(tempthreshpond))
                    .and(srtmslope.lte(maxpondslope));

// warm small ice patches below certain slope are actually also water
var water = water.where(smallpatches
                        .and(tsurfC.gt(tempthreshpond))
                        .and(srtmslope.lte(maxpondslope)),1)

// rock outcrop classification
var rock  = debris.and(srtmslope.gt(maxdebslope));
var rock  = debris.multiply(0);

// final classification merge
var cls = debris.multiply(0);
var cls = cls.where(ice,0)
             .where(debris,1)
             .where(water,2)
             .where(rock,3)
             .add(1)
             .reproject({crs: ndsi.projection(), scale: 30});

// get steep debris
var steepdeb = ee.Terrain.slope(srtm30).gt(24).and(cls.eq(2));
var steepdeb = steepdeb.updateMask(steepdeb);





// make map =======================================================================================

var plotcols = ['lightblue','gray','darkblue', 'brown'];
var plotnames = ['Debris-free ice','Debris-covered ice','Supraglacial pond', 'Debris on steep slope'];

// apply tight rgi mask again for display
Map.addLayer(ee.Image(),{},'Placeholder', 0);
Map.addLayer(ordered.updateMask(rgimask), {bands: ['B6','B5','B4'], min: ['0.0','0.0','0.0'], max: ['0.8','0.8','0.8']}, 'LS8 Imagery (RGB=SWIR-NIR-RED)', 0);
Map.addLayer(tsurfC.updateMask(rgimask),   {min: '-2', max: '30', palette: ['darkblue','purple','green','yellow','orange','red']}, 'LS8 Uncal. Surface Temperature', 0);
Map.addLayer(ndsi.updateMask(rgimask), {min: '-0.4', max: '0.8', palette: ['darkred','yellow','green','blue','darkblue']}, 'Snow Index', 0);
Map.addLayer(wateridx.updateMask(rgimask),{},'Water index',0);
Map.addLayer(cls.updateMask(rgimask), {min: '1', max: '3', palette: plotcols.slice(0,3)},'Debris classification', 1, 1);
Map.addLayer(steepdeb.updateMask(rgimask), {palette: plotcols[3]}, 'Debris on steep slope', 1, 1)

// get number of layers in map
var numlayers = Map.layers().length();



// make data inspection panel =====================================================================

// Initiate widget panel
var inspectpanel = ui.Panel();
inspectpanel.style().set('width', '350px');

// Create an intro panel with labels.
var intro = ui.Panel([
  ui.Label({
    value: 'Click map to inspect',
    style: {fontSize: '15px', fontWeight: 'bold'}
  }),
]);
inspectpanel.add(intro);

// Function to run on map click
Map.onClick(function(coords) {
  
  // Add a red dot for the point clicked on.
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var dot = ui.Map.Layer(point, {color: 'FF0000'}, 'Inspection point');
  Map.layers().set(numlayers, dot);
  
  // get raster values to add to panels
  var elev_samp = srtm30.reduceRegion(ee.Reducer.mode(),point).get('elevation');
  var rgi_samp  = rgimask.reduceRegion(ee.Reducer.mode(),point).get('b1_max');
  var tir_samp  = ee.Number(tsurfC.reduceRegion(ee.Reducer.mean(),point).get('B10'));
  var ndsi_samp = ee.Number(ndsi.reduceRegion(ee.Reducer.mean(),point,30).get('nd'));
  
  // Create panels that display lon/lat/elev values.
  var lon = ui.Label('lon: ' + coords.lon.toFixed(3),{fontSize: '12px'});
  var lat = ui.Label('lat: ' + coords.lat.toFixed(3),{fontSize: '12px'});
  var elev = ui.Label(undefined,{fontSize: '12px'});
  elev_samp.evaluate(function(x){elev.setValue('z: ' + x)});
  var coordpan = ui.Panel([lon, lat, elev], ui.Panel.Layout.flow('horizontal'));
  inspectpanel.widgets().set(1, coordpan);

  // Create panels that displays RGI id values.
  var rgi = ui.Label(undefined,{fontSize: '12px'});
  rgi_samp.evaluate(function(x){rgi.setValue('Randolph Glacier Inventory ID: '+x)});
  var rgipan = ui.Panel([rgi],ui.Panel.Layout.flow('vertical'));
  inspectpanel.widgets().set(2, rgipan);

  // Create spectral profile
  var chartimg = ordered.select(['B[1-7]']);
  var bandChart = ui.Chart.image.regions({
    image: chartimg,
    regions: point,
    scale: 30,
    seriesProperty: 'label'
  });
  bandChart.setChartType('LineChart');
  bandChart.setOptions({
    title: 'Spectral profile',
    hAxis: {title: 'Landsat 8 band'},
    vAxis: {title: 'TOA Reflectance'},
    lineWidth: 1,
    pointSize: 2,
    series: {0: {color: 'darkblue'}}
  });
  inspectpanel.widgets().set(3, bandChart);
  
  // Create panels that displays temperature
  var bt = ui.Label(undefined,{fontSize: '12px'});
  tir_samp.evaluate(function(x){bt.setValue('Tsurf: '+x.toFixed(1)+' °C')});
  var btpan = ui.Panel([bt],ui.Panel.Layout.flow('horizontal'));
  inspectpanel.widgets().set(4, btpan);
  
  // Create panels that displays NDSI
  var ndsilab = ui.Label(undefined,{fontSize: '12px'});
  ndsi_samp.evaluate(function(x){ndsilab.setValue('NDSI: '+x.toFixed(2))});
  var ndsipan = ui.Panel([ndsilab],ui.Panel.Layout.flow('horizontal'));
  inspectpanel.widgets().set(5, ndsipan);
});

Map.style().set('cursor', 'crosshair');

// Add the panel to the ui.root.
ui.root.insert(0, inspectpanel);



// panel for quick access to specific locations ========================================================

var glaclist = {// name lat lon zoom
               Baltoro: [76.360,35.750,10],
               BaraShigri: [77.663,32.183,11],
               Chongce: [81.151,35.311,11],
               Engilchek: [79.883,42.163,11],
               Fedchenko: [72.279,38.767,10],
               Hispar: [75.317,36.085,10],
               Khumbu: [86.838,27.982,12],
               Langtang: [85.712,28.294,12],
               Ngozumpa: [86.694,27.991,12],
               Siachen: [77.136,35.390,10]
};

var reglist    = rgiregions.toList(100);
var namelist   = reglist.map(function(x){return ee.Feature(x).get('Label')});
var regiondict = ee.Dictionary.fromLists(namelist,reglist).getInfo();

var selectpanel = ui.Panel({
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {position: 'top-center', padding: '0', margin: '0'}
});
var selectregion = ui.Select({
  items: Object.keys(regiondict),
  onChange: function(key) {
    Map.centerObject(ee.Feature(regiondict[key]));
    var regmap = ui.Map.Layer(ee.Feature(regiondict[key]), {color: '#FFFF55'}, key+' outline', 1, 0.67);
    Map.layers().set(0, regmap)
  },
  style: {width: '140px'},
  placeholder: '[Navigate to region...]',
});
var selectglac = ui.Select({
  items: Object.keys(glaclist),
  onChange: function(key) {
    Map.setCenter(glaclist[key][0], glaclist[key][1],glaclist[key][2]);
  },
  style: {width: '140px'},
  placeholder: '[Navigate to glacier...]',
});

// Add panel to map
selectpanel.add(selectregion).add(selectglac);
Map.add(selectpanel);



// add a grayscale basemap  =============================================================

var styles = {'Custom': [
  {featureType: 'all',stylers: [{saturation: -75},{lightness:50}]},
  {featureType: 'all', elementType: "labels", stylers: [{visibility: 'off' }]},
  {featureType: 'road', stylers: [{visibility: 'off' }]},
  {featureType: 'transit', stylers: [{visibility: 'off' }]},
  {featureType: 'poi', stylers: [{visibility: 'off' }]}
  ]};
Map.setOptions(null, styles);

// center map on full region after creating all panels
Map.centerObject(rgihull);



// add legend    =============================================================

var legend = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px'
  }
});
var legendTitle = ui.Label({
  value: 'Glacier classification',
  style: {
    fontWeight: 'bold',
    fontSize: '16px',
    margin: '0 0 10px 0',
    padding: '0 0 0 0'
  }
});
legend.add(legendTitle);
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
  // Create the label that is actually the colored box.
  var colorBox = ui.Label({
    style: {
      backgroundColor: color,
      // Use padding to give the box height and width.
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });
  // Create the label filled with the description text.
  var description = ui.Label({
    value: name,
    style: {margin: '0 0 4px 6px'}
  });
  return ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
};
for (var i = 0; i < plotcols.length; i++) {
    legend.add(makeRow(plotcols[i], plotnames[i]));
  }
Map.add(legend);

