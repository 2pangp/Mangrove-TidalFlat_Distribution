#模型构建
library(tidyverse)
library(terra)
library(raster)
library(randomForest)
library(ranger)
library(gbm)
library(caret)
library(dplyr)
library(ggplot2)

tidal = rast("E:/tidalflat/zhuanyi/tidal1996.tif")
tidal = tidal*0

PA = rast("E:/tidalflat/tidal/PA.tif")
PA1 = crop(PA, tidal, mask = T)
pa = merge(PA1,tidal,first = T)

writeRaster(pa, "E:/tidalflat/inputdata/PA.tif", overwrite=T)

plot(pa)
PA = rast("E:/Mangrove/DataMap/PAs.tif")
hf = rast("E:/tidalflat/result/human17_19.tif")
city = rast("E:/tidalflat/result/city.tif")
road = rast("G:/data/road/road22.tif")
river = rast("E:/NUS/new data/river1.tif")
coastline = rast("E:/NUS/new data/coastline2.tif")
loss_prev = rast("E:/tidalflat/tidal1km/tidal2019_lossprev.tif")



PA = resample(PA, tidal, method = 'average')
road  = project(road, crs(tidal))
road = resample(road, tidal, method = 'average')
river = resample(river, tidal, method = 'average')
coastline = resample(coastline, tidal, method = 'average')
plot(hf)

PA = crop(PA,tidal,mask = T)
hf = crop(hf,tidal,mask = T)
city = crop(city,tidal,mask = T)
road = crop(road,tidal,mask = T)
river = crop(river,tidal,mask = T)
coastline = crop(coastline,tidal,mask = T)
loss_prev = crop(loss_prev,tidal,mask = T)

PA = merge(PA,tidal,first = T)
hf = merge(hf,tidal,first = T)
city = merge(city,tidal,first = T)
road = merge(road,tidal,first= T)
river = merge(river,tidal,first = T)
coastline = merge(coastline,tidal,first = T)
loss_prev = merge(loss_prev,tidal,first= T)

writeRaster(PA, "E:/tidalflat/inputdata/PA.tif", overwrite=T)
writeRaster(hf, "E:/tidalflat/inputdata/hf.tif", overwrite=T)
writeRaster(city, "E:/tidalflat/inputdata/city.tif", overwrite=T)
writeRaster(road, "E:/tidalflat/inputdata/road.tif", overwrite=T)
writeRaster(river, "E:/tidalflat/inputdata/river.tif", overwrite=T)
writeRaster(coastline, "E:/tidalflat/inputdata/coastline.tif", overwrite=T)
writeRaster(loss_prev, "E:/tidalflat/inputdata/loss_prev.tif", overwrite=T)

FileName = list.files("E:/tidalflat/agriculture/future3")
for (f in FileName)
{
  rastTemp = rast(paste0("E:/tidalflat/agriculture/future3/", f))
  rastTemp= lapp(rastTemp, fun = function(x){ifelse((x < 0), 0, x)})
  writeRaster(rastTemp, paste0("E:/tidalflat/agriculture/future4/", f), overwrite = T)
  print(f)
}



loss_prev = rast("E:/tidalflat/inputdata/loss_prev.tif")
hf = rast("E:/tidalflat/inputdata/hf.tif")
PA = rast("E:/tidalflat/inputdata/PA.tif")
city = rast("E:/tidalflat/inputdata/city.tif")
road = rast("E:/tidalflat/inputdata/road.tif")
river = rast("E:/tidalflat/inputdata/river.tif")
coastline = rast("E:/tidalflat/inputdata/coastline.tif")
gdp = rast("E:/tidalflat/agriculture/ssp3/GDP2100_ssp5.tif")
pop = rast("E:/tidalflat/agriculture/ssp3/SSP5_2100.tif")
coconut = rast("E:/tidalflat/agriculture/suitability3/coconut_rcp85_2100.tif")
coffee = rast("E:/tidalflat/agriculture/suitability3/coffee_rcp85_2100.tif")
olp = rast("E:/tidalflat/agriculture/suitability3/olp_rcp85_2100.tif")
rice = rast("E:/tidalflat/agriculture/suitability3/rice_rcp85_2100.tif")
rubber = rast("E:/tidalflat/agriculture/suitability3/rubber_rcp85_2100.tif")
sugarcane = rast("E:/tidalflat/agriculture/suitability3/sugarcane_rcp85_2100.tif")
barley = rast("E:/tidalflat/agriculture/future4/2100_Barley_RCP8.5.tif")
maize = rast("E:/tidalflat/agriculture/future4/2100_Maize_RCP8.5.tif")
rapeseed = rast("E:/tidalflat/agriculture/future4/2100_Rapeseed_RCP8.5.tif")
sorghum = rast("E:/tidalflat/agriculture/future4/2100_Sorghum_RCP8.5.tif")
soybean = rast("E:/tidalflat/agriculture/future4/2100_Soybean_RCP8.5.tif")
sunflower = rast("E:/tidalflat/agriculture/future4/2100_Sunflower_RCP8.5.tif")
wheat = rast("E:/tidalflat/agriculture/future4/2100_Wheat_RCP8.5.tif")


pred_list = stack(raster(gdp),raster(hf), raster(pop), raster(city), raster(road), raster(PA), 
                  raster(coconut), raster(coffee), raster(olp), raster(rice), raster(rubber), 
                  raster(sugarcane),raster(loss_prev),raster(river), raster(coastline),
                  raster(barley),raster(maize),raster(rapeseed),raster(sorghum),raster(soybean),
                  raster(sunflower),raster(wheat))

pred_list = c(gdp, hf, pop, city, road, PA, coconut, coffee, olp, rice, rubber, sugarcane,loss_prev,river,
              coastline,barley,maize,rapeseed,sorghum,soybean,sunflower,wheat)
names(pred_list) = c("gdp", "hf", "pop", "city", "road", "PA", "coconut", "coffee", "olp", "rice", 
                     "rubber", "sugarcane","loss_prev","river","coastline","barley","maize",
                     "rapeseed","sorghum","soybean","sunflower","wheat")


writeRaster(pred_list, "E:/tidalflat/inputdata/2100ssp5pred_list.tif",overwrite =T)
plot(pred_list)

coastline = rast("E:/NUS/new data/coastline2.tif")
tidal = rast("E:/tidalflat/zhuanyi/tidal1996.tif")
coastline = resample(coastline,tidal)


RFmodel = readRDS("E:/tidalflat/data/累积/RF1.rds")
pred_set = rast("E:/tidalflat/inputdata/pred/2030ssp1pred_list.tif")
PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "E:/tidalflat/data/累积/2030ssp1pred.tif",overwrite =T)

pred_set = rast("E:/tidalflat/inputdata/pred/2030ssp5pred_list.tif")
PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "E:/tidalflat/data/累积/2030ssp5pred.tif",overwrite =T)

pred_set = rast("E:/tidalflat/inputdata/pred/2050ssp1pred_list.tif")
PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "E:/tidalflat/data/累积/2050ssp1pred_list.tif",overwrite =T)

pred_set = rast("E:/tidalflat/inputdata/pred/2050ssp5pred_list.tif")
PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "E:/tidalflat/data/累积/2050ssp5pred_list.tif",overwrite =T)

pred_set = rast("E:/tidalflat/inputdata/pred/2100ssp1pred_list.tif")
PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "E:/tidalflat/data/累积/2100ssp1pred_list.tif",overwrite =T)

pred_set = rast("E:/tidalflat/inputdata/pred/2100ssp5pred_list.tif")
PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "E:/tidalflat/data/累积/2100ssp5pred_list.tif",overwrite =T)
plot(PredRas)


rasTemp = rast(ncol = , nrow = , extent = c(-180, 180, -90, 90), crs = crs(Country))

pa3 = vect("E:/tidalflat/tidalmian/PA_3GE1.shp")
pa012 = vect("E:/tidalflat/tidalmian/PA012.shp")

PA = rbind(pa012,pa3)
writeVector(PA,"E:/tidalflat/tidalmian/PA.shp")
pa_iucn = PA[,'IUCN_CAT']

unique(pa_iucn$IUCN_CAT)

pa_cat = function(iucn_cat)
{
  if (is.na(iucn_cat)){cat_val = 8}
  else if (iucn_cat == 'Ia'){cat_val = 1}
  else if (iucn_cat == 'Ib'){cat_val = 2}
  else if (iucn_cat == 'II'){cat_val = 3}
  else if (iucn_cat == 'III'){cat_val = 4}
  else if (iucn_cat == 'IV'){cat_val = 5}
  else if (iucn_cat == 'V'){cat_val = 6}
  else if (iucn_cat == 'VI'){cat_val = 7}
  else {cat_val = 8}
  
  return (cat_val)
}

pa_cat(NA)
is.na(NA)

test = data.frame(pa_iucn$IUCN_CAT)
is.na(test$pa_iucn.IUCN_CAT)

cat_value = apply(test, MARGIN = 1, FUN = pa_cat)
pa_iucn$cat_value = cat_value

R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(tidal))
pa_ras = rasterize(pa_iucn, R, field = "cat_value")
plot(pa_ras)
writeRaster(pa_ras, "E:/tidalflat/tidal/PA.tif",overwrite=TRUE)



