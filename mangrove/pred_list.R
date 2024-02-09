library(terra)

Country = vect("E:/NUS/WorldCountries/WorldCountriesWithISO.shp")
Country$val0 = 0
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Country))
Country_rast = rasterize(Country, R, field = "val0")

plot(Country_rast)

hf17 = rast("G:/NUS/humanfootprint/hfp2017.tif")
hf18 = rast("G:/NUS/humanfootprint/hfp2018.tif")
hf19 = rast("G:/NUS/humanfootprint/hfp2019.tif")
hf17_19 = (hf17+hf18+hf19)/3
writeRaster(hf17_19,"E:/Mangrove/crop/data3/P1719/hf17_19.tif")

Loss2007 = merge(Loss07, Country_rast, first = T)
Loss2010 = merge(Loss10, Country_rast, first = T)
Loss2013 = merge(Loss13, Country_rast, first = T)

Loss04 = rast("E:/NUS/data/Loss2004merge.tif")
Loss07 = rast("E:/NUS/data/Loss2007merge.tif")
Loss10 = rast("E:/NUS/data/Loss2010merge.tif")
Loss13 = rast("E:/NUS/data/Loss2013merge.tif")
Loss16 = rast("E:/NUS/mian/2016loss22.tif")
Loss19 = rast("E:/NUS/mian/2019loss11.tif")

Loss2004 = resample(Loss04, PA, method = 'average')
Loss2007 = resample(Loss07, PA, method = 'average')
Loss2010 = resample(Loss10, PA, method = 'average')
Loss2013 = resample(Loss13, PA, method = 'average')
plot(Loss2013)
writeRaster(Loss2004,"E:/NUS/data/Loss2004_10km.tif")
writeRaster(Loss2007,"E:/NUS/data/Loss2007_10km.tif")
writeRaster(Loss2010,"E:/NUS/data/Loss2010_10km.tif")
writeRaster(Loss2013,"E:/NUS/data/Loss2013_10km.tif")



Loss04 = vect("E:/NUS/mian/2004mian3.shp")
Loss07 = vect("E:/NUS/mian/2007mian3.shp")
Loss10 = vect("E:/NUS/mian/2010mian3.shp")
Loss13 = vect("E:/NUS/mian/2013mian3.shp")
Loss16 = vect("E:/NUS/mian/2016mian3.shp")
Loss19 = vect("E:/NUS/mian/2019mian3.shp")

R = rast(ext = c(-180, 180, -90, 90), ncol = 43200, nrow = 21600)
R = rast(ext(Country_rast), ncol = ncol(Country_rast), nrow = nrow(Country_rast))
Loss2004 = rasterize(Loss04, R, field = "gridcode", background = NA, touches = T)
plot(Loss2004)
writeRaster(Loss2004,"E:/NUS/mangroveloss/Loss2004.tif")


hf2017 = rast("G:/NUS/humanfootprint/hfp2017.tif")
hf2018 = rast("G:/NUS/humanfootprint/hfp2018.tif")
hf2019 = rast("G:/NUS/humanfootprint/hfp2019.tif")
hf = (hf2017+hf2018+hf2019)/3
writeRaster(hf,"G:/data/human/human17_19.tif")

PA = rast("E:/Mangrove/crop/PA.tif")
hf = rast("G:/data/human/human17_19.tif")
Mangrove1996 = vect("G:/NUS/GMW_001_GlobalMangroveWatch/01_Data/GMW_1996_v2.shp")
Country$val0 = 0
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Country))
Country_rast = rasterize(Country, R, field = "val0")

hf = project(hf, crs(PA))
hf = resample(hf, PA, method = 'average')
hf = crop(hf, Mangrove1996, mask = T)
writeRaster(hf,"E:/Mangrove/crop/hf.tif",overwrite = T)

Loss2013 = rast("E:/NUS/mian/2013loss1111pro.tif")
Loss20131 = resample(Loss2013, PA, method = 'average')
writeRaster(Loss20131 ,"E:/Mangrove/crop/Loss2013.tif",overwrite = T)
plot(PA)

ext(Country_rast) = ext(Loss20131)
test = merge(Loss20131, Country_rast)
plot(test)

Loss19 = project(Loss19, crs(PA))
Mangrove1996 = vect("G:/NUS/GMW_001_GlobalMangroveWatch/01_Data/GMW_1996_v2.shp")
Mangrove1996 $val0 = 0
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Mangrove1996))
Mangrove = rasterize(Mangrove1996, R, field = "val0")
writeRaster(Mangrove, "E:/Mangrove/crop/Mangrove.tif", overwrite=T)
plot(Mangrove)

loss = vect("E:/NUS/mian/99_19mian.shp")
loss $val0 = 0
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Mangrove1996))
lossmian = rasterize(loss, R, field = "val0")
writeRaster(lossmian, "E:/Mangrove/crop/lossmian.tif", overwrite=T)


GDP = rast("G:/data/gdp_Kummu/GDP_per_capita_PPP_1990_2015_v2.nc")
#14_15
GDP14_15 = GDP[[25:26]]
GDP14_15_mean = mean(GDP14_15)

GDP = rast("G:/data/GDP_real1/GDP17_19.tif")
PA = rast("E:/Mangrove/DataMap/PAs.tif")
hf = rast("G:/data/data17_19/human17_18.tif")
pop = rast("G:/data/data17_19/Pop17_19.tif")
city = rast("G:/data/Accessibility to cities/accessibility_to_cities_2015_v1.0/accessibility_to_cities_2015_v1.0.tif")
road = rast("G:/data/road/road22.tif")
coconut = rast("G:/data/suitability/Coconut2010.tif")
coffee = rast("G:/data/suitability/Coffee2010.tif")
olp = rast("G:/data/suitability/olp2010.tif")
rice = rast("G:/data/suitability/rice2010.tif")
rubber= rast("G:/data/suitability/Rubber2010.tif")
sugarcane= rast("G:/data/suitability/Sugarcane2010.tif")

PA = resample(PA, Loss2019, method = 'average')
hf = project(hf, crs(Loss2019))
hf = resample(hf, Loss2019, method = 'average')
pop = project(pop, crs(Loss2019))
pop = resample(pop, Loss2019, method = 'average')
road = project(road, crs(Loss2019))
road = resample(road, Loss2019, method = 'average')
city = project(city, crs(Loss2019))
city = resample(city, Loss2019, method = 'average')
PA1 = crop(PA,Loss2019,mask = T)
hf1=crop(hf,Loss2019,mask = T)
pop1=crop(pop,Loss2019,mask = T)
road1=crop(road,Loss2019,mask = T)
city1=crop(city,Loss2019,mask = T)

writeRaster(PA1, "E:/Mangrove/crop/data3/PA.tif", overwrite=T)
writeRaster(hf1, "E:/Mangrove/crop/data3/hf.tif", overwrite=T)
writeRaster(pop1, "E:/Mangrove/crop/data3/pop.tif", overwrite=T)
writeRaster(road1, "E:/Mangrove/crop/data3/road.tif", overwrite=T)
writeRaster(city1, "E:/Mangrove/crop/data3/city.tif", overwrite=T)

coconut = project(coconut, crs(Country_rast))
coconut = resample(coconut, Country_rast, method = 'bilinear')
coffee = project(coffee, crs(Country_rast))
coffee = resample(coffee, Country_rast, method = 'bilinear')
olp = project(olp, crs(Country_rast))
olp = resample(olp, Country_rast, method = 'bilinear')
rice = project(rice, crs(Country_rast))
rice = resample(rice, Country_rast, method = 'bilinear')
rubber = project(rubber, crs(Country_rast))
rubber = resample(rubber, Country_rast, method = 'bilinear')
sugarcane = project(sugarcane, crs(Country_rast))
sugarcane = resample(sugarcane, Country_rast, method = 'bilinear')
gdp = project(GDP, crs(Country_rast))
gdp = resample(GDP, Country_rast, method = 'average')
plot(hf)
Mangrove1996 = vect("G:/NUS/GMW_001_GlobalMangroveWatch/01_Data/GMW_1996_v2.shp")
loss_prev1 = crop(loss_prev, Mangrove1996, mask = T)
writeRaster(loss_prev1,"E:/Mangrove/crop/loss_prev.tif")



ssp1gdp2040 = rast("G:/data/GDPfuture/GDP_SSP1_1km/GDP_SSP1_1km/GDP2040_ssp1.tif")
ssp1gdp2050 = rast("G:/data/GDPfuture/GDP_SSP1_1km/GDP_SSP1_1km/GDP2050_ssp1.tif")
ssp1gdp20401 = crop(ssp1gdp2040, Mangrove1996, mask = T)
ssp1gdp20501 = crop(ssp1gdp2050, Mangrove1996, mask = T)
writeRaster(ssp1gdp20401,"E:/Mangrove/crop/ssp1gdp2040.tif")
writeRaster(ssp1gdp20501,"E:/Mangrove/crop/ssp1gdp2050.tif")


gdp = rast("E:/Mangrove/crop/popssp_data/SSP1GDP.tif")
pop = rast("E:/Mangrove/crop/popssp_data/SSP1Pop.tif")
hf = rast("E:/Mangrove/crop/hf.tif")
city = rast("E:/Mangrove/crop/city.tif")
road = rast("E:/Mangrove/crop/road.tif")
PA = rast("E:/Mangrove/crop/PA.tif")
coconut = rast("E:/Mangrove/crop/popssp_data11/coconut_rcp26_2030.tif")
coffee = rast("E:/Mangrove/crop/popssp_data11/coffee_rcp26_2030.tif")
olp = rast("E:/Mangrove/crop/popssp_data11/olp_rcp26_2030.tif")
rice = rast("E:/Mangrove/crop/popssp_data11/rice_rcp26_2030.tif")
rubber = rast("E:/Mangrove/crop/popssp_data11/rubber_rcp26_2030.tif")
sugarcane = rast("E:/Mangrove/crop/popssp_data11/sugarcane_rcp26_2030.tif")
loss_prev = rast("E:/Nus/data/Loss2019crop.tif")

coconut = resample(coconut,loss_prev,method = 'average')
coffee = resample(coffee,loss_prev,method = 'average')
olp = resample(olp,loss_prev,method = 'average')
rice  = resample(rice,loss_prev,method = 'average')
rubber = resample(rubber,loss_prev,method = 'average')
sugarcane= resample(sugarcane,loss_prev,method = 'average')

pred_list = stack(raster(gdp),  raster(hf), raster(pop), raster(city), raster(road), raster(PA), raster(coconut), raster(coffee), raster(olp), raster(rice), raster(rubber), raster(sugarcane),raster(loss_prev))
pred_list = c(gdp, hf, pop, city, road, PA, coconut, coffee, olp, rice, rubber, sugarcane,loss_prev)
names(pred_list) = c("gdp", "hf", "pop", "city", "road", "PA", "coconut", "coffee", "olp", "rice", "rubber", "sugarcane","loss_prev")

writeRaster(pred_list, "E:/Mangrove/data3/crop_data/2030ssp1pred_list.tif")

plot(pred_list)



FileName = list.files("E:/Mangrove/crop/popssp_data/")
for (f in FileName)
{
  rastTemp = rast(paste0("E:/Mangrove/crop/popssp_data/", f))
  rastTemp = crop(rastTemp, Loss2019, mask = T)
  writeRaster(rastTemp, paste0("E:/Mangrove/crop/data3/", f), overwrite = T)
  print(f)
}

yanmo = rast("E:/Mangrove/crop/data2/Loss2019.tif")

library(terra)
Loss04 = rast("E:/NUS/data2/Loss2004.tif")
Loss07 = rast("E:/NUS/data2/Loss2007.tif")
Loss10 = rast("E:/NUS/data2/Loss2010.tif")
Loss13 = rast("E:/NUS/data2/Loss2013.tif")
Loss16 = rast("E:/NUS/data2/Loss2016.tif")
Loss19 = rast("E:/NUS/data2/Loss2019.tif")

loss_prev_04 = focal(Loss04, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_07 = focal(Loss07, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_10 = focal(Loss10, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_13 = focal(Loss13, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_16 = focal(Loss16, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_19 = focal(Loss19, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 

loss_prev2004=crop(loss_prev_04,Loss19,mask = T)
loss_prev2007=crop(loss_prev_07,Loss19,mask = T)
loss_prev2010=crop(loss_prev_10,Loss19,mask = T)
loss_prev2013=crop(loss_prev_13,Loss19,mask = T)
loss_prev2016=crop(loss_prev_16,Loss19,mask = T)
loss_prev2019=crop(loss_prev_19,Loss19,mask = T)

writeRaster(loss_prev2004, "E:/Mangrove/crop/loss_prev_04.tif", overwrite=T)
writeRaster(loss_prev2007, "E:/Mangrove/crop/loss_prev_07.tif", overwrite=T)
writeRaster(loss_prev2010, "E:/Mangrove/crop/loss_prev_10.tif", overwrite=T)
writeRaster(loss_prev2013, "E:/Mangrove/crop/loss_prev_13.tif", overwrite=T)
writeRaster(loss_prev2016, "E:/Mangrove/crop/loss_prev_16.tif", overwrite=T)
writeRaster(loss_prev2019, "E:/Mangrove/crop/loss_prev_19.tif", overwrite=T)

plot(loss_prev2004)

hf = rast("E:/Mangrove/crop/hf2.tif")
loss_prev = rast("E:/Mangrove/crop/loss_prev_2019.tif")
pre_set = rast("E:/Mangrove/data3/crop_data/2100ssp5pred_set.tif")
names(pre_set) = c("gdp", "hf", "pop", "city", "road", "PA", "coconut", "coffee", "olp", "rice", "rubber", "sugarcane","loss_prev")
pre_set$hf =hf 
pre_set$loss_prev = loss_prev 
writeRaster(pre_set, "E:/Mangrove/data3/crop_data/2100ssp5pred_set1.tif", overwrite=T)
plot(pre_set)

Loss19 = rast("E:/NUS/data2/Loss2019.tif")
PA = rast("./DataMap/PAs.tif")
plot(pa)
pa = resample(PA,Loss19,method = 'average')
pa = crop(pa,Loss19,mask = T)
writeRaster(pa, "E:/Mangrove/crop/data3/crop/pa.tif", overwrite=T)

pre_set = rast("E:/Mangrove/data3/pred_set.tif")
names(pre_set) = c("gdp", "hf", "pop", "city", "road", "PA", "coconut", "coffee", "olp", "rice", "rubber", "sugarcane","loss_prev")
loss_prev= rast("E:/Mangrove/crop/loss_prev_16.tif")
loss_prev1=crop(loss_prev,Mangrove1996,mask = T)
pre_set$loss_prev = loss_prev1 
writeRaster(pre_set, "E:/Mangrove/data3/pre_set.tif", overwrite=T)
plot(pre_set$loss_prev)

Loss2019 = rast("E:/NUS/data2/Loss2019.tif")
pop= rast("E:/Mangrove/crop/data3/P1719/Pop17_19.tif")

pop =crop(pop,Loss19,mask = T)
writeRaster(pop, "E:/Mangrove/crop/data3/xiu/pop.tif", overwrite=T)

river= rast("E:/NUS/new data/river.tif")
coastline= rast("E:/NUS/new data/coastline1.tif")
river1=crop(river,Loss2019,mask = T)
coastline1=crop(coastline,Loss2019,mask = T)

writeRaster(river1, "E:/Mangrove/crop/data3/river.tif", overwrite=T)
writeRaster(coastline1, "E:/Mangrove/crop/data3/coastline.tif", overwrite=T)



Loss19 = rast("E:/Mangrove/crop/loss_prev_19int.tif")
Loss19 = Loss19*0

FileName = list.files("E:/Mangrove/crop/data3/merge/")
for (f in FileName)
{
  rastTemp = rast(paste0("E:/Mangrove/crop/data3/merge/", f))
  rastTemp= crop(rastTemp, Loss19,mask=T)
  writeRaster(rastTemp, paste0("E:/Mangrove/crop/data3/crop/", f), overwrite = T)
  print(f)
}

Data = read.csv("./data3/p55/data_NEW3.csv")
Data[is.na(Data$PA), c('PA')] = 0
Data[Data$PA == 9, c('PA')] = 8





pa1 = vect("G:/WDPA_Feb2023_Public_shp_0/WDPA_Feb2023_Public_shp-polygons0.shp")
pa2 = vect("G:/WDPA_Feb2023_Public_shp_1/WDPA_Feb2023_Public_shp-polygons1.shp")
pa3 = vect("G:/WDPA_Feb2023_Public_shp_2/WDPA_Feb2023_Public_shp-polygons2.shp")
pa_chn = vect("G:/Additional data for protected area_Apr2023/China_CAS+WDPA_2019/PA_Union.shp")
pa_ind = vect("G:/Additional data for protected area_Apr2023/India_WDPA_Oct2018/India_PA.shp")
#pa_tur = vect("G:/Additional data for protected area_Apr2023/Turkey_OSM_Apr2023_nature_reserve/turkey_NR.shp")

samp07 = vect("E:/NUS/sample/p55/2007sample.shp")
samp10 = vect("E:/NUS/sample/p55/2010sample.shp")
samp13 = vect("E:/NUS/sample/p55/2013sample.shp")
samp16 = vect("E:/NUS/sample/p55/2016sample.shp")
samp19 = vect("E:/NUS/sample/p55/2019sample.shp") 
samp = rbind(samp07, samp10, samp13, samp16, samp19)

pa_wdpa = rbind(pa1, pa2, pa3)
pa_wdpa = pa_wdpa[,c('ISO3', 'IUCN_CAT', 'DESIG_TYPE')]
pa_chn = pa_chn[,c('ISO3', 'IUCN_CAT', 'DESIG_TYPE')]
pa_ind = pa_ind[,c('ISO3', 'IUCN_CAT', 'DESIG_TYPE')]

pa = rbind(pa_wdpa, pa_chn, pa_ind)
pa_iucn = pa[,'IUCN_CAT']

pa_iucn
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

R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Country))
pa_ras = rasterize(pa_iucn, R, field = "cat_value")
writeRaster(pa_ras,"E:/Mangrove/data3/p55/pa_iucn.tif")
plot(pa_ras)

Loss2019 = rast("E:/NUS/data2/Loss2019.tif")

pa_ras= rast("E:/Mangrove/data3/p55/pa_iucn.tif")
pa = merge(pa_ras,Loss19,first=T)
pa1 = crop(pa,Loss2019,mask=T)
writeRaster(pa1, "E:/Mangrove/crop/data3/crop/PA.tif", overwrite=T)

Data$PA = c(terra::extract(pa_ras, samp07)[,2],
            terra::extract(pa_ras, samp10)[,2],
            terra::extract(pa_ras, samp13)[,2],
            terra::extract(pa_ras, samp16)[,2],
            terra::extract(pa_ras, samp19)[,2])
write.csv(Data, "./data3/p55/data_NEW3.csv", row.names = F)
sum(is.na(Data$PA))

pa = rast("E:/Mangrove/crop/data3/crop/PA.tif")

plot(pa_ras)

gdpssp1= rast("E:/ssp/GDP2030_ssp1.tif")
popssp1= rast("E:/ssp/SSP1_2030.tif")

gdpssp1 = resample(gdpssp1,Loss2019)
popssp1 = resample(popssp1,Loss2019)
gdp = crop(gdpssp1,Loss2019,mask=T)
pop = crop(popssp1,Loss2019,mask=T)
writeRaster(gdp, "E:/Mangrove/crop/data3/crop/gdp.tif", overwrite=T)
writeRaster(pop, "E:/Mangrove/crop/data3/crop/pop.tif", overwrite=T)


FileName = list.files("E:/Mangrove/agriculture3/")
for (f in FileName)
{
  rastTemp = rast(paste0("E:/Mangrove/agriculture3/", f))
  rastTemp= scale(rastTemp)
  writeRaster(rastTemp, paste0("E:/Mangrove/ssp2/", f), overwrite = T)
  print(f)
}

