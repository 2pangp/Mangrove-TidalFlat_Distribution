# setwd("")

# Library
library(raster)
library(terra)

library(randomForest)
library(ranger)
library(caret)

samp04 = shapefile("D:/tpp/dian2/dian/2004sample.shp")
samp07 = shapefile("D:/tpp/dian2/dian/2007sample.shp")
samp10 = shapefile("D:/tpp/dian2/dian/2010sample.shp")
samp13 = shapefile("D:/tpp/dian2/dian/2013sample.shp")
samp16 = shapefile("D:/tpp/dian2/dian/2016sample.shp") 
samp19 = shapefile("D:/tpp/dian2/dian/2019sample.shp") 

GDPData = matrix(ncol = 1, nrow = 100000)
GDPData[,1] = raster::extract(GDP99_04_mean, samp16, buffer=2000, fun=mean, na.rm = T)

colnames(GDPData) = c("GDP99_04", "GDP05_07", "GDP08_10", "GDP11_13", "GDP14_15")
write.csv(GDPData, "./data2/GDP2016.csv", row.names = F)
sum(is.na(GDPData))

Country_rast = rast("D:/tpp/data/world300m.tif")
Country = vect("D:/tpp/WorldCountries/WorldCountriesWithISO.shp")
Country$val0 = 0
R = rast(ext = ext(Mangrove1996), ncol = ncol(Country_rast), nrow = nrow(Country_rast))
Country_rast = rasterize(Country, R, field = "val0")
Country_rast = project(Country_rast, crs(Country))
plot(Country_rast)

SSP1GDP = rast("D:/tpp/future/GDP2030_ssp1.tif")
SSP2GDP = rast("D:/tpp/future/GDP2030_ssp2.tif")
SSP3GDP = rast("D:/tpp/future/GDP2030_ssp3.tif")
SSP4GDP = rast("D:/tpp/future/GDP2030_ssp4.tif")
SSP5GDP = rast("D:/tpp/future/GDP2030_ssp5.tif")
Mangrove1996 = vect("D:/tpp/GMW_001_GlobalMangroveWatch/01_Data/GMW_1996_v2.shp")
Mangrove1996 $val0 = 0
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Mangrove1996))
Mangrove = rasterize(Mangrove1996, R, field = "val0")

Loss04= vect("D:/tpp/mian/2004mian3.shp")
Loss07= vect("D:/tpp/mian/2007mian3.shp")
Loss10= vect("D:/tpp/mian/2010mian3.shp")
Loss13= vect("D:/tpp/mian/2013mian3.shp")
Loss16= vect("D:/tpp/mian/2016mian3.shp")
Loss19= vect("D:/tpp/mian/2019mian33.shp")

loss2007 = merge(Loss04,Loss07)  
loss2010 = merge(Loss04,Loss07,Loss10)
loss2013 = merge(Loss04,Loss07,Loss10,Loss13)
loss2016 = merge(Loss04,Loss07,Loss10,Loss13,Loss16)
loss2019 = merge(Loss04,Loss07,Loss10,Loss13,Loss16,Loss19)

Mangrove1996 $val0 = 0
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Mangrove1996))
Mangrove = rasterize(Mangrove1996, R, field = "val0")


SSP1GDP2 = crop(SSP1GDP, Mangrove1996, mask = T)
SSP2GDP2 = crop(SSP2GDP, Mangrove1996, mask = T)
SSP3GDP2 = crop(SSP3GDP, Mangrove1996, mask = T)
SSP4GDP2 = crop(SSP4GDP, Mangrove1996, mask = T)
SSP5GDP2 = crop(SSP5GDP, Mangrove1996, mask = T)

writeRaster(SSP1GDP2, "D:/tpp/data2/SSP1GDP.tif")
writeRaster(SSP2GDP2, "D:/tpp/data2/SSP2GDP.tif")
writeRaster(SSP3GDP2, "D:/tpp/data2/SSP3GDP.tif")
writeRaster(SSP4GDP2, "D:/tpp/data2/SSP4GDP.tif")
writeRaster(SSP5GDP2, "D:/tpp/data2/SSP5GDP.tif")


SSP1Pop = rast("D:/tpp/future/SSP1_2030.tif")
SSP2Pop = rast("D:/tpp/future/SSP2_2030.tif")
SSP3Pop = rast("D:/tpp/future/SSP3_2030.tif")
SSP4Pop = rast("D:/tpp/future/SSP4_2030.tif")
SSP5Pop = rast("D:/tpp/future/SSP5_2030.tif")

SSP1Pop2 = crop(SSP1Pop, Mangrove1996, mask = T)
SSP2Pop2 = crop(SSP2Pop, Mangrove1996, mask = T)
SSP3Pop2 = crop(SSP3Pop, Mangrove1996, mask = T)
SSP4Pop2 = crop(SSP4Pop, Mangrove1996, mask = T)
SSP5Pop2 = crop(SSP5Pop, Mangrove1996, mask = T)

writeRaster(SSP1Pop2,"D:/tpp/data2/SSP1Pop.tif")
writeRaster(SSP2Pop2,"D:/tpp/data2/SSP2Pop.tif")
writeRaster(SSP3Pop2,"D:/tpp/data2/SSP3Pop.tif")
writeRaster(SSP4Pop2,"D:/tpp/data2/SSP4Pop.tif")
writeRaster(SSP5Pop2,"D:/tpp/data2/SSP5Pop.tif")


coconut = rast("D:/tpp/data2/coconut.tif")

Loss04 = rast("D:/tpp/loss3/2004reclass.tif")
Loss07 = rast("D:/tpp/loss3/2007reclass.tif")
Loss10 = rast("D:/tpp/loss3/2010reclass.tif")
Loss13 = rast("D:/tpp/loss3/2013reclass.tif")

samp07 = vect("D:/tpp/sample/2007sample.shp")
samp10 = vect("D:/tpp/sample/2010sample.shp")
samp13 = vect("D:/tpp/sample/2013sample.shp")
samp16 = vect("D:/tpp/sample/2016sample.shp") 

coconut = terra::project(coconut,crs(Loss04))
Loss04 = resample(Loss04,coconut,method='average')
Loss07 = resample(Loss07,coconut,method='average')
Loss10 = resample(Loss10,coconut,method='average')
Loss13 = resample(Loss13,coconut,method='average')

Loss04 = terra::project(Loss04,csr(coconut))
Loss07 = terra::project(Loss07,csr(coconut))
Loss10 = terra::project(Loss10,csr(coconut))
Loss13 = terra::project(Loss13,csr(coconut))

writeRaster(Loss04,"D:/tpp/sample/Loss04_10km.tif")
writeRaster(Loss07,"D:/tpp/sample/Loss07_10km.tif")
writeRaster(Loss10,"D:/tpp/sample/Loss10_10km.tif")
writeRaster(Loss13,"D:/tpp/sample/Loss13_10km.tif")

Data = as.data.frame(matrix(ncol=1,nrow = 200000))
colnames(Data) = c("loss_pre")
Data$loss_pre = c(terra::extract(Loss04,sample07),
                  terra::extract(Loss07,sample10),
                  terra::extract(Loss10,sample13),
                  terra::extract(Loss13,sample16))

coconut = rast("D:/tpp/data2/coconut.tif")
Loss2004 = rast("D:/tpp/mian/200403pro.tif")
Loss2007 = rast("D:/tpp/mian/200703pro_re1.tif")
Loss2010 = rast("D:/tpp/mian/201003pro_re1.tif")
Loss2013 = rast("D:/tpp/mian/201303pro.tif")

Loss04 = resample(Loss2004, coconut, method = 'average')
Loss07 = resample(Loss2007, coconut, method = 'average')
Loss10 = resample(Loss2010, coconut, method = 'average')
Loss13 = resample(Loss2013, coconut, method = 'average')
writeRaster(Loss07,"D:/tpp/data2//Loss2007.tif",overwrite = T)
writeRaster(Loss10,"D:/tpp/data2//Loss2010.tif",overwrite = T)
writeRaster(Loss13,"D:/tpp/data2//Loss2013.tif")

Country = rast("D:/tpp/data/world300m.tif")
Loss2004 = rast("D:/tpp/data/2004loss222pro.tif")
Loss2007 = rast("D:/tpp/data/2007loss1111pro.tif")
Loss2010 = rast("D:/tpp/data/2010loss1111pro.tif")
Loss2013 = rast("D:/tpp/data/2013loss1111pro.tif")
Loss2016 = rast("D:/tpp/data/2016loss22.tif")
Loss2019 = rast("D:/tpp/data/2019loss11.tif")
Country = vect("D:/tpp/WorldCountries/WorldCountriesWithISO.shp")
Country$val0 = 0
R=rast(ncol=129600, nrow=64800,ext=c(-180,180,-90,90),crs=crs(Country))
Country_rast = rasterize(Country, R, field = "val0")

plot(Country_rast)


Loss04 = resample(Loss2004,R)
Loss07 = resample(Loss2007,R)
Loss10 = resample(Loss2010,R)
Loss13 = resample(Loss2013,R)
Loss16 = resample(Loss2016,R)
Loss19 = resample(Loss2019,R)
writeRaster(Loss19,"D:/tpp/data/Loss2019.tif")

writeRaster(Loss07,"D:/tpp/data/Loss2007.tif")
writeRaster(Loss10,"D:/tpp/data/Loss2010.tif")
writeRaster(Loss13,"D:/tpp/data/Loss2013.tif")

Loss04 = rast("D:/tpp/data/Loss2004.tif")

plot(Country_rast)

Loss20191 = merge(Loss19, Country_rast, first = T)
writeRaster(Loss20191,"D:/tpp/data/Loss2019merge.tif")
coconut = rast("D:/tpp/data2/coconut.tif")
Loss2019 = resample(Loss20191, coconut, method = 'average')
writeRaster(Loss2019,"D:/tpp/data/Loss2019_1km.tif")
plot(Loss2019_crop)
Loss2019_crop = crop(Loss2019, Mangrove1996, mask = T)
writeRaster(Loss2019_crop,"D:/tpp/data/Loss2019crop.tif")

coconut = rast("D:/tpp/crop/coconut_rcp26_2030.tif")


FileName = list.files("D:/tpp/crop")
for (f in FileName)
{
  rastTemp = rast(paste0("D:/tpp/crop/", f))
  rastTemp = crop(rastTemp, Mangrove1996, mask = T)
  writeRaster(rastTemp, paste0("D:/tpp/popssp_data/", f), overwrite = T)
  print(f)
}

FileName = list.files("D:/tpp/data/gdp_ssp1")
for (f in FileName)
{
  rastTemp = rast(paste0("D:/tpp/data/gdp_ssp1/", f))
  rastTemp = crop(rastTemp, Mangrove1996, mask = T)
  writeRaster(rastTemp, paste0("D:/tpp/popssp_data/", f), overwrite = T)
  print(f)
}

rasTemp = rast(ncol = 43200, nrow = 21600 , extent = c(-180, 180, -90, 90), crs = crs(Loss04))

Loss04 = rast("D:/tpp/mangroveloss/2004mangrove.tif")
Loss07 = rast("D:/tpp/mangroveloss/2007mangrove.tif")
Loss10 = rast("D:/tpp/mangroveloss/2010mangrove.tif")
Loss13 = rast("D:/tpp/mangroveloss/2013mangrove.tif")
Loss16 = rast("D:/tpp/mangroveloss/2016mangrove.tif")
Loss19 = rast("D:/tpp/mangroveloss/2019mangrove.tif")
plot(Loss2019)

Loss2004 = resample(Loss04, rasTemp, method = 'average')
Loss2007 = resample(Loss07, rasTemp, method = 'average')
Loss2010 = resample(Loss10, rasTemp, method = 'average')
Loss2013 = resample(Loss13, rasTemp, method = 'average')
Loss2016 = resample(Loss16, rasTemp, method = 'average')
Loss2019 = resample(Loss19, rasTemp, method = 'average')

writeRaster(Loss2004,"D:/tpp/data2/Loss2004.tif")
writeRaster(Loss2007,"D:/tpp/data2/Loss2007.tif")
writeRaster(Loss2010,"D:/tpp/data2/Loss2010.tif")
writeRaster(Loss2013,"D:/tpp/data2/Loss2013.tif")
writeRaster(Loss2016,"D:/tpp/data2/Loss2016.tif")
writeRaster(Loss2019,"D:/tpp/data2/Loss2019.tif")


Mangrove = read.csv("D:/tpp/mangrove.csv")
head(Mangrove)
set.seed(10)
randidx = sample(nrow(Mangrove), 0.8*nrow(Mangrove))
train_Mangrove= Mangrove[randidx,]
test_Mangrove = Mangrove[-randidx,]
str(Mangrove)
sum(train_Mangrove$loss)
table(test_Mangrove$loss)

mod_f = as.formula(loss ~ gdp  + hf  + pop  + city + road + PA + coconut + coffee + olp + 
                     rice + rubber + sugarcane+loss_prev+coconut_prod + coffee_prod + 
                     olp_prod + rice_prod + rubber_prod+river+coastline)

RFmodel = randomForest(mod_f,data= train_Mangrove, ntrees = 300, mtry = 21/3, type = regression, importance = T)
predRF = predict(RFmodel, newdata = test_Mangrove)
RMSE(test_Mangrove$loss, RFPred$predictions)
cor(train_Mangrove$loss, predict(RFmodel, train_Mangrove, n.trees = 300)$predictions) #RSQ_Trn
cor(test_Mangrove$loss, predict(RFmodel, test_Mangrove, n.trees = 300)$predictions) #RSQ_Tst


coconut= vect("D:/tpp/production/coconut.shp")
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(coconut))


Loss19 = rast("D:/tpp/mangroveloss/Loss2019.tif")
Loss19 = Loss19*0


FileName = list.files("D:/tpp/crop/")
for (f in FileName)
{
  rastTemp = rast(paste0("D:/tpp/crop/", f))
  rastTemp = resample(rastTemp, Loss19)
  writeRaster(rastTemp, paste0("D:/tpp/agriculture1/", f), overwrite = T)
  print(f)
}


FileName = list.files("D:/tpp/agriculture1/")
for (f in FileName)
{
  rastTemp = rast(paste0("D:/tpp/agriculture1/", f))
  rastTemp = merge(rastTemp, Loss19,first = T)
  writeRaster(rastTemp, paste0("D:/tpp/agriculture2/", f), overwrite = T)
  print(f)
}


FileName = list.files("D:/tpp/agriculture2/")
for (f in FileName)
{
  rastTemp = rast(paste0("D:/tpp/agriculture2/", f))
  rastTemp = crop(rastTemp, Loss19, mask = T)
  writeRaster(rastTemp, paste0("D:/tpp/agriculture3/", f), overwrite = T)
  print(f)
}

Loss04 = rast("D:/tpp/tidal/2004tidal01.tif")
Loss07 = rast("D:/tpp/tidal/2007tidal01.tif")
Loss10 = rast("D:/tpp/tidal/2010tidal01.tif")
Loss13 = rast("D:/tpp/tidal/2013tidal01.tif")
Loss16 = rast("D:/tpp/tidal/2016tidal01.tif")
Loss19 = rast("D:/tpp/tidal/2019tidal01.tif")
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(Loss04))


loss_prev_04 = focal(Loss2004, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_07 = focal(loss20071, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_10 = focal(loss20101, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_13 = focal(loss20131, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_16 = focal(loss20161, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 
loss_prev_19 = focal(loss20191, w=matrix(1,11,11), fun= mean,na.rm=T, pad=T, padValue=NA) 

writeRaster(loss_prev_04, "D:/tpp/data2/lossprev/loss_prev_04.tif", overwrite=T)
writeRaster(loss_prev_07, "D:/tpp/data2/lossprev/loss_prev_07.tif", overwrite=T)
writeRaster(loss_prev_10, "D:/tpp/data2/lossprev/loss_prev_10.tif", overwrite=T)
writeRaster(loss_prev_13, "D:/tpp/data2/lossprev/loss_prev_13.tif", overwrite=T)
writeRaster(loss_prev_16, "D:/tpp/data2/lossprev/loss_prev_16.tif", overwrite=T)
writeRaster(loss_prev_19, "D:/tpp/data2/lossprev/loss_prev_19.tif", overwrite=T)

Loss04 = rast("D:/tpp/tidal/nlossotidal.tif")
Sample_Loss =spatSample(Loss04, 100000,exhaustive = T,as.points = T,na.rm = T)
colnames(Sample_Loss) = c("CID")

writeVector(Sample_Loss,"D:/tpp/tidal/noloss.shp")

plot(Sample_Loss)
