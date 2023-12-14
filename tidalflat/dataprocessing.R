# Library
library(terra)

setwd("E:/tidalflat")


samp07 = vect("E:/tidalflat/tidalpoint/2007sample.shp")
samp10 = vect("E:/tidalflat/tidalpoint/2010sample.shp")
samp13 = vect("E:/tidalflat/tidalpoint/2013sample.shp")
samp16 = vect("E:/tidalflat/tidalpoint/2016sample.shp")
samp19 = vect("E:/tidalflat/tidalpoint/2019sample.shp") 

samp07[,c('lon', 'lat')] = geom(samp07, df = T)[,c('x', 'y')]
samp10[,c('lon', 'lat')] = geom(samp10, df = T)[,c('x', 'y')]
samp13[,c('lon', 'lat')] = geom(samp13, df = T)[,c('x', 'y')]
samp16[,c('lon', 'lat')] = geom(samp16, df = T)[,c('x', 'y')]
samp19[,c('lon', 'lat')] = geom(samp19, df = T)[,c('x', 'y')]

samp = rbind(samp07, samp10, samp13, samp16, samp19)

Data = cbind(Data, data.frame(samp[,c('lon', 'lat')]))
write.csv(Data, "./data3/p55/Data.csv", row.names = F)

#数据提取
CityAccess = rast("E:/tidalflat/focalresult/city.tif")
RoadAccess = rast("G:/data/road/road22.tif")
RoadAccess = project(RoadAccess, crs(CityAccess))

Access = data.frame(matrix(ncol = 5, nrow = 200000))
Access[,1] = terra::extract(CityAccess, samp07, buffer=2000, fun=mean, na.rm = T)
Access[,2] = terra::extract(CityAccess, samp10, buffer=2000, fun=mean, na.rm = T)
Access[,3] = terra::extract(CityAccess, samp13, buffer=2000, fun=mean, na.rm = T)
Access[,4] = terra::extract(CityAccess, samp16, buffer=2000, fun=mean, na.rm = T)
Access[,5] = terra::extract(CityAccess, samp19, buffer=2000, fun=mean, na.rm = T)
colnames(Access) = c("07","10","13","16","19")
write.csv(Access, "E:/tidalflat/data/CityAccess.csv", row.names = F)

Access1 = data.frame(matrix(ncol = 5, nrow = 200000))
Access1[,1] = terra::extract(RoadAccess, samp07, buffer=2000, fun=mean, na.rm = T)
Access1[,2] = terra::extract(RoadAccess, samp10, buffer=2000, fun=mean, na.rm = T)
Access1[,3] = terra::extract(RoadAccess, samp13, buffer=2000, fun=mean, na.rm = T)
Access1[,4] = terra::extract(RoadAccess, samp16, buffer=2000, fun=mean, na.rm = T)
Access1[,5] = terra::extract(RoadAccess, samp19, buffer=2000, fun=mean, na.rm = T)
colnames(Access1) = c("07","10","13","16","19")
write.csv(Access1, "E:/tidalflat/data/RoadAccess.csv", row.names = F)

#HFP
R=rast(ncol = 43200, nrow = 21600,ext=c(-180,180,-90,90),crs=crs(gdp))

human07 = rast("G:/data/human/human05_07.tif")
human10 = rast("G:/data/human/human08_10.tif")
human13 = rast("G:/data/human/human11_13.tif")
human16 = rast("G:/data/human/human14_16.tif")
human19 = rast("G:/data/human/human17_19.tif")
human071 = project(human07, crs = crs(R))
human101 = project(human10, crs = crs(R))
human131 = project(human13, crs = crs(R))
human161 = project(human16, crs = crs(R))
human191 = project(human19, crs = crs(R))
human = data.frame(matrix(ncol = 5, nrow = 200000))
human[,1] = terra::extract(human071, samp07, buffer=2000, fun=mean, na.rm = T)
human[,2] = terra::extract(human101, samp10, buffer=2000, fun=mean, na.rm = T)
human[,3] = terra::extract(human131, samp13, buffer=2000, fun=mean, na.rm = T)
human[,4] = terra::extract(human161, samp16, buffer=2000, fun=mean, na.rm = T)
human[,5] = terra::extract(human191, samp19, buffer=2000, fun=mean, na.rm = T)

colnames(human) = c("07","10","13","16","19")
write.csv(human, "E:/tidalflat/data/human.csv", row.names = F)

Pop07 = rast("G:/data/Pop/Pop05_07.tif")
Pop10 = rast("G:/data/Pop/Pop08_10.tif")
Pop13 = rast("G:/data/Pop/Pop11_13.tif")
Pop16 = rast("G:/data/Pop/Pop14_16.tif")
Pop19 = rast("G:/data/Pop/Pop17_19.tif")
Pop = data.frame(matrix(ncol = 5, nrow = 200000))
Pop[,1] = terra::extract(Pop07, samp07, buffer=2000, fun=mean, na.rm = T)
Pop[,2] = terra::extract(Pop10, samp10, buffer=2000, fun=mean, na.rm = T)
Pop[,3] = terra::extract(Pop13, samp13, buffer=2000, fun=mean, na.rm = T)
Pop[,4] = terra::extract(Pop16, samp16, buffer=2000, fun=mean, na.rm = T)
Pop[,5] = terra::extract(Pop19, samp19, buffer=2000, fun=mean, na.rm = T)

colnames(Pop) = c("07","10","13","16","19")
write.csv(Pop, "E:/tidalflat/data/Pop.csv", row.names = F)

GDP07 = rast("G:/data/GDP_real1/GDP05_07.tif")
GDP10 = rast("G:/data/GDP_real1/GDP08_10.tif")
GDP13 = rast("G:/data/GDP_real1/GDP11_13.tif")
GDP16 = rast("G:/data/GDP_real1/GDP14_16.tif")
GDP19 = rast("G:/data/GDP_real1/GDP17_19.tif")

GDP07 = focal(GDP07, w=matrix(1,31,31), fun= sum,na.rm=T, pad=T, padValue=NA) 


GDP = data.frame(matrix(ncol = 5, nrow = 200000))
GDP[,1] = terra::extract(GDP07, samp07, buffer=2000, fun=mean, na.rm = T)
GDP[,2] = terra::extract(GDP10, samp10, buffer=2000, fun=mean, na.rm = T)
GDP[,3] = terra::extract(GDP13, samp13, buffer=2000, fun=mean, na.rm = T)
GDP[,4] = terra::extract(GDP16, samp16, buffer=2000, fun=mean, na.rm = T)
GDP[,5] = terra::extract(GDP19, samp19, buffer=2000, fun=mean, na.rm = T)

colnames(GDP) = c("07","10","13","16","19")
write.csv(GDP, "E:/tidalflat/data/GDP.csv", row.names = F)

#数据合并创建表格
CityAccess = read.csv("./data/CityAccess.csv")
RoadAccess = read.csv("./data/RoadAccess.csv")
human = read.csv("./data/human.csv")
GDP = read.csv("./data/GDP.csv")
Pop = read.csv("./data/Pop.csv")
Data = as.data.frame(matrix(ncol = 6, nrow = 1000000))
colnames(Data) = c("year", "gdp", "hf", "pop", "city", "road")
Data$year = c(rep(2007, 100000), rep(2010, 100000), rep(2013, 100000), rep(2016, 100000), rep(2019, 100000))
Data$gdp = c(GDP$X07,GDP$X10,GDP$X13,GDP$X16,GDP$X19)
Data$hf = c(human$X07,human$X10,human$X13,human$X16,human$X19)
Data$pop = c(Pop$X07,Pop$X10,Pop$X13,Pop$X16,Pop$X19)
Data$city = c(CityAccess$X07,CityAccess$X10,CityAccess$X13,CityAccess$X16,CityAccess$X19)
Data$road =  c(RoadAccess$X07,RoadAccess$X10,RoadAccess$X13,RoadAccess$X16,RoadAccess$X19)


PA = rast("E:/tidalflat/tidal/PA.tif")
plot(PA)
Data$PA= c(terra::extract(PA, samp07)[,2],
               terra::extract(PA, samp10)[,2],
               terra::extract(PA, samp13)[,2],
               terra::extract(PA, samp16)[,2],
               terra::extract(PA, samp19)[,2])

Data$PA[is.na(Data$PA)]=0


Coconut = rast("G:/data/suitability/Coconut2010.tif")
Coffee = rast("G:/data/suitability/Coffee2010.tif")
Olp = rast("G:/data/suitability/olp2010.tif")
Rice = rast("G:/data/suitability/rice2010.tif")
Rubber= rast("G:/data/suitability/Rubber2010.tif")
Sugarcane= rast("G:/data/suitability/Sugarcane2010.tif")

Coconut07 = terra::extract(Coconut, samp07)[,2]
Coconut10 = terra::extract(Coconut, samp10)[,2]
Coconut13 = terra::extract(Coconut, samp13)[,2]
Coconut16 = terra::extract(Coconut, samp16)[,2]
Coconut19 = terra::extract(Coconut, samp19)[,2]
Data$coconut = c(Coconut07,Coconut10,Coconut13,Coconut16,Coconut19)
Coffee07 = terra::extract(Coffee, samp07)[,2]
Coffee10 = terra::extract(Coffee, samp10)[,2]
Coffee13 = terra::extract(Coffee, samp13)[,2]
Coffee16 = terra::extract(Coffee, samp16)[,2]
Coffee19 = terra::extract(Coffee, samp19)[,2]
Data$coffee = c(Coffee07,Coffee10,Coffee13,Coffee16,Coffee19)
Olp07 = terra::extract(Olp, samp07)[,2]
Olp10 = terra::extract(Olp, samp10)[,2]
Olp13 = terra::extract(Olp, samp13)[,2]
Olp16 = terra::extract(Olp, samp16)[,2]
Olp19 = terra::extract(Olp, samp19)[,2]
Data$olp = c(Olp07,Olp10,Olp13,Olp16,Olp19)
Rice07 = terra::extract(Rice, samp07)[,2]
Rice10 = terra::extract(Rice, samp10)[,2]
Rice13 = terra::extract(Rice, samp13)[,2]
Rice16 = terra::extract(Rice, samp16)[,2]
Rice19 = terra::extract(Rice, samp19)[,2]
Data$rice = c(Rice07,Rice10,Rice13,Rice16,Rice19)
Rubber07 = terra::extract(Rubber, samp07)[,2]
Rubber10 = terra::extract(Rubber, samp10)[,2]
Rubber13 = terra::extract(Rubber, samp13)[,2]
Rubber16 = terra::extract(Rubber, samp16)[,2]
Rubber19 = terra::extract(Rubber, samp19)[,2]
Data$rubber = c(Rubber07,Rubber10,Rubber13,Rubber16,Rubber19)
Sugarcane07 = terra::extract(Sugarcane, samp07)[,2]
Sugarcane10 = terra::extract(Sugarcane, samp10)[,2]
Sugarcane13 = terra::extract(Sugarcane, samp13)[,2]
Sugarcane16 = terra::extract(Sugarcane, samp16)[,2]
Sugarcane19 = terra::extract(Sugarcane, samp19)[,2]
Data$sugarcane = c(Sugarcane07,Sugarcane10,Sugarcane13,Sugarcane16,Sugarcane19)

Barley = rast("E:/tidalflat/agriculture/Barley.tif")
Maize = rast("E:/tidalflat/agriculture/Maize.tif")
Rapeseed = rast("E:/tidalflat/agriculture/Rapeseed.tif")
Sorghum = rast("E:/tidalflat/agriculture/Sorghum.tif")
Soybean = rast("E:/tidalflat/agriculture/Soybean.tif")
Sunflower = rast("E:/tidalflat/agriculture/Sunflower.tif")
Wheat = rast("E:/tidalflat/agriculture/Wheat.tif")
Data$Barley= c(terra::extract(Barley, samp07)[,2],
               terra::extract(Barley, samp10)[,2],
               terra::extract(Barley, samp13)[,2],
               terra::extract(Barley, samp16)[,2],
               terra::extract(Barley, samp19)[,2])
Data$Maize= c(terra::extract(Maize, samp07)[,2],
               terra::extract(Maize, samp10)[,2],
               terra::extract(Maize, samp13)[,2],
               terra::extract(Maize, samp16)[,2],
               terra::extract(Maize, samp19)[,2])
Data$Rapeseed= c(terra::extract(Rapeseed, samp07)[,2],
              terra::extract(Rapeseed, samp10)[,2],
              terra::extract(Rapeseed, samp13)[,2],
              terra::extract(Rapeseed, samp16)[,2],
              terra::extract(Rapeseed, samp19)[,2])
Data$Sorghum= c(terra::extract(Sorghum, samp07)[,2],
                 terra::extract(Sorghum, samp10)[,2],
                 terra::extract(Sorghum, samp13)[,2],
                 terra::extract(Sorghum, samp16)[,2],
                 terra::extract(Sorghum, samp19)[,2])
Data$Soybean= c(terra::extract(Soybean, samp07)[,2],
                 terra::extract(Soybean, samp10)[,2],
                 terra::extract(Soybean, samp13)[,2],
                 terra::extract(Soybean, samp16)[,2],
                 terra::extract(Soybean, samp19)[,2])
Data$Sunflower= c(terra::extract(Sunflower, samp07)[,2],
                terra::extract(Sunflower, samp10)[,2],
                terra::extract(Sunflower, samp13)[,2],
                terra::extract(Sunflower, samp16)[,2],
                terra::extract(Sunflower, samp19)[,2])
Data$Wheat = c(terra::extract(Wheat, samp07)[,2],
                  terra::extract(Wheat, samp10)[,2],
                  terra::extract(Wheat, samp13)[,2],
                  terra::extract(Wheat, samp16)[,2],
                  terra::extract(Wheat, samp19)[,2])

Loss04 = rast("E:/tidalflat/tidal1km/tidal2004_lossprev.tif")
Loss07 = rast("E:/tidalflat/tidal1km/tidal2007_lossprev.tif")
Loss10 = rast("E:/tidalflat/tidal1km/tidal2010_lossprev.tif")
Loss13 = rast("E:/tidalflat/tidal1km/tidal2013_lossprev.tif")
Loss16 = rast("E:/tidalflat/tidal1km/tidal2016_lossprev.tif")
Data$loss_prev = c(terra::extract(Loss04, samp07)[,2],
               terra::extract(Loss07, samp10)[,2],
               terra::extract(Loss10, samp13)[,2],
               terra::extract(Loss13, samp16)[,2],
               terra::extract(Loss16, samp19)[,2])

samp2007 = vect("E:/tidalflat/tidalpoint/2007sample_Intersect1.shp")
samp2010 = vect("E:/tidalflat/tidalpoint/2010sample_Intersect1.shp")
samp2013 = vect("E:/tidalflat/tidalpoint/2013sample_Intersect1.shp")
samp2016 = vect("E:/tidalflat/tidalpoint/2016sample_Intersect1.shp")
samp2019 = vect("E:/tidalflat/tidalpoint/2019sample_Intersect1.shp")
Data$ISO3 = c(samp2007$ISO3,samp2010$ISO3,samp2013$ISO3,samp2016$ISO3,samp2019$ISO3)

river = rast("E:/NUS/new data/river1.tif")
coastline = rast("E:/NUS/new data/coastline2.tif")
Data$river =  c(terra::extract(river, samp07)[,2],
                terra::extract(river, samp10)[,2],
                terra::extract(river, samp13)[,2],
                terra::extract(river, samp16)[,2],
                terra::extract(river, samp19)[,2])
Data$coastline =  c(terra::extract(coastline, samp07)[,2],
                    terra::extract(coastline, samp10)[,2],
                    terra::extract(coastline, samp13)[,2],
                    terra::extract(coastline, samp16)[,2],
                    terra::extract(coastline, samp19)[,2])
sum(is.na(Data$coastline))

GDP07 = rast("E:/tidalflat/focalresult/GDP05_07.tif")
GDP10 = rast("E:/tidalflat/focalresult/GDP08_10.tif")
GDP13 = rast("E:/tidalflat/focalresult/GDP11_13.tif")
GDP16 = rast("E:/tidalflat/focalresult/GDP14_16.tif")
GDP19 = rast("E:/tidalflat/focalresult/GDP17_19.tif")
Pop07 = rast("E:/tidalflat/focalresult/Pop05_07.tif")
Pop10 = rast("E:/tidalflat/focalresult/Pop08_10.tif")
Pop13 = rast("E:/tidalflat/focalresult/Pop11_13.tif")
Pop16 = rast("E:/tidalflat/focalresult/Pop14_16.tif")
Pop19 = rast("E:/tidalflat/focalresult/Pop17_19.tif")
human07 = rast("E:/tidalflat/focalresult/human05_07.tif")
human10 = rast("E:/tidalflat/focalresult/human08_10.tif")
human13 = rast("E:/tidalflat/focalresult/human11_13.tif")
human16 = rast("E:/tidalflat/focalresult/human14_16.tif")
human19 = rast("E:/tidalflat/focalresult/human17_19.tif")
city = rast("E:/tidalflat/focalresult/city.tif")

Data$gdp = c(terra::extract(GDP07, samp07)[,2],
             terra::extract(GDP10, samp10)[,2],
             terra::extract(GDP13, samp13)[,2],
             terra::extract(GDP16, samp16)[,2],
             terra::extract(GDP19, samp19)[,2])
Data$pop = c(terra::extract(Pop07, samp07)[,2],
             terra::extract(Pop10, samp10)[,2],
             terra::extract(Pop13, samp13)[,2],
             terra::extract(Pop16, samp16)[,2],
             terra::extract(Pop19, samp19)[,2])
Data$hf = c(terra::extract(human07, samp07)[,2],
             terra::extract(human10, samp10)[,2],
             terra::extract(human13, samp13)[,2],
             terra::extract(human16, samp16)[,2],
             terra::extract(human19, samp19)[,2])
Data$city = c(terra::extract(city, samp07)[,2],
            terra::extract(city, samp10)[,2],
            terra::extract(city, samp13)[,2],
            terra::extract(city, samp16)[,2],
            terra::extract(city, samp19)[,2])
tidal2004 = rast("E:/tidalflat/tidal/tidal2004.tif")
tidal2007 = rast("E:/tidalflat/tidal/tidal2007.tif")
tidal2010 = rast("E:/tidalflat/tidal/tidal2010.tif")
tidal2013 = rast("E:/tidalflat/tidal/tidal2013.tif")
tidal2016 = rast("E:/tidalflat/tidal/tidal2016.tif")
tidal2019 = rast("E:/tidalflat/tidal/tidal2019.tif")
Data$loss = c(terra::extract(tidal2007, samp07)[,2],
              terra::extract(tidal2010, samp10)[,2],
              terra::extract(tidal2013, samp13)[,2],
              terra::extract(tidal2016, samp16)[,2],
              terra::extract(tidal2019, samp19)[,2])


tidal2007 = rast("E:/tidalflat/tidal1km/tidal20071km.tif")
tidal2010 = rast("E:/tidalflat/tidal1km/tidal20101km.tif")
tidal2013 = rast("E:/tidalflat/tidal1km/tidal20131km.tif")
tidal2016 = rast("E:/tidalflat/tidal1km/tidal20161km.tif")
tidal2019 = rast("E:/tidalflat/tidal1km/tidal20191km.tif")

Data$loss = c(terra::extract(tidal2007, samp07)[,2],
              terra::extract(tidal2010, samp10)[,2],
              terra::extract(tidal2013, samp13)[,2],
              terra::extract(tidal2016, samp16)[,2],
              terra::extract(tidal2019, samp19)[,2])


sum(is.na(Data$PA))
head(Data)
write.csv(Data, "./data/累积/data_raw.csv", row.names = F)


###
Data = read.csv("./data/data_raw.csv")
sum(is.na(Data$city))
Data$Barley[Data$Barley < '0']<- 0
Data$Maize[Data$Maize < '0']<- 0
Data$Rapeseed[Data$Rapeseed < '0']<- 0
Data$Sorghum[Data$Sorghum < '0']<- 0
Data$Soybean[Data$Soybean < '0']<- 0
Data$Sunflower[Data$Sunflower < '0']<- 0
Data$Wheat[Data$Wheat < '0']<- 0

colnames(Data)
Data = na.omit(Data)

IQR_gdp = IQR(Data$gdp, na.rm = T)
IQR_hf = IQR(Data$hf, na.rm = T)
IQR_pop = IQR(Data$pop, na.rm = T)

Q_gdp = quantile(Data$gdp, probs = c(0.25, 0.75), na.rm = T)
Q_hf = quantile(Data$hf, probs = c(0.25, 0.75), na.rm = T)
Q_pop = quantile(Data$pop, probs = c(0.25, 0.75), na.rm = T)

Data = Data[ - which((Data$gdp > Q_gdp[2] + 5 * IQR_gdp) | (Data$gdp < Q_gdp[1] - 5 * IQR_gdp) |
                    (Data$pop > Q_pop[2] + 5 * IQR_pop) | (Data$pop < Q_pop[1] - 5 * IQR_pop) ),]


#Data = Data[ - which((Data$gdp_var > Q_gdp_var[2] + 9 * IQR_gdp_var) | (Data$gdp_var < Q_gdp_var[1] - 9 * IQR_gdp_var)),]
#Data = Data[ - which((Data$hf > Q_hf[2] + 5 * IQR_hf) | (Data$hf < Q_hf[1] - 5 * IQR_hf)),]
#Data = Data[ - which((Data$hf_var > Q_hf_var[2] + 5 * IQR_hf_var) | (Data$hf_var < Q_hf_var[1] - 5 * IQR_hf_var)),]
#Data = Data[ - which((Data$pop > Q_pop[2] + 5 * IQR_pop) | (Data$pop < Q_pop[1] - 5 * IQR_pop)),]
#Data = Data[ - which((Data$pop_var > Q_pop_var[2] + 5 * IQR_pop_var) | (Data$pop_var < Q_pop_var[1] - 5 * IQR_pop_var)),]
#Data = Data[ - which((Data$city > Q_city[2] + 1.5 * IQR_city) | (Data$city < Q_city[1] - 1.5 * IQR_city)),]
#Data = Data[ - which((Data$road > Q_road[2] + 1.5 * IQR_road) | (Data$road < Q_road[1] - 1.5 * IQR_road)),]
write.csv(Data,"./data/累积/Data1.csv", row.names = F)

#模型构建
library(tidyverse)
library(terra)
library(raster)
library(randomForest)
library(ranger)
library(gbm)
library(caret)
tidalflat = read.csv("./data/累积/Data1.csv")
colnames(tidalflat)

tidalflat = tidalflat[c("loss","gdp", "hf", "pop", "city", "road","PA", "coconut", "coffee", "olp", 
                      "rice", "rubber", "sugarcane","loss_prev","river","coastline","Barley", 
                      "Maize","Rapeseed","Sorghum","Soybean","Sunflower","Wheat")]
colnames(tidalflat) = c("loss", "gdp", "hf", "pop", "city", "road", "PA", "coconut", "coffee", "olp", "rice", 
                       "rubber", "sugarcane","loss_prev","river","coastline","barley","maize",
                       "rapeseed","sorghum","soybean","sunflower","wheat")
#01 data split
head(tidalflat)
set.seed(10)
randidx = sample(nrow(tidalflat), 0.8*nrow(tidalflat))
train_tidalflat= tidalflat[randidx,]
test_tidalflat = tidalflat[-randidx,]
write.csv(test_tidalflat, "E:/tidalflat/inputdata/importance/test_tidalflat.csv", row.names = F)

mod_f = as.formula(loss ~ gdp  + hf  + pop  + city + road + PA + coconut + coffee + olp + 
                     rice + rubber + sugarcane+loss_prev+river+coastline+
                     barley+maize+rapeseed+sorghum +soybean+sunflower+wheat)

RFmodel = ranger(mod_f,data= train_tidalflat, num.trees = 300, mtry = 22/3,importance = 'impurity', verbose = T)
saveRDS(RFmodel, "./data/累积/RF1.rds")

# 获取特征的重要性排序
importance(RFmodel)
pdp = partial(RFmodel,
              train_tidalflat,
              pred.var = 'loss_prev')
line(pdp$loss_prev, pdp$lossratio, type = 'l')
plot(pdp)

###GrowPDSI
pdp_loss_prev <- partial(
  RFmodel,
  pred.var = 'wheat',
  progress = T,
  approx = T,
  grid.resolution = 30,
  plot = F
)
plot(pdp_loss_prev,type = 'l')
write.csv(pdp_loss_prev, "E:/tidalflat/inputdata/importance/wheat.csv", row.names = F)


###testing
RFmodel = readRDS("./data/累积/RF1.rds")
RFPred = predict(RFmodel, data = test_tidalflat)
RMSE(test_tidalflat$loss, RFPred$predictions)
cor(train_tidalflat$loss, predict(RFmodel, train_tidalflat, n.trees = 300)$predictions) #RSQ_Trn
cor(test_tidalflat$loss, predict(RFmodel, test_tidalflat, n.trees = 300)$predictions) #RSQ_Tst
write.csv(RFPred, "E:/tidalflat/inputdata/importance/RFPred.csv", row.names = F)

#yuce
RFmodel = readRDS("E:/tidalflat/data/累积/RF1.rds")
pred_set = rast("E:/tidalflat/inputdata/2100ssp5pred_list.tif")

PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "E:/tidalflat/inputdata/2100ssp5pred.tif",overwrite =T)

plot(PredRas)
summary(RFmodel)

brtMan = gbm(
  loss ~ gdp  + hf  + pop  + city + road + PA + coconut + coffee + olp + 
    rice + rubber + sugarcane+loss_prev+river+coastline,
  data = train_Mangrove,
  distribution = "bernoulli",
  n.trees = 500,
  interaction.depth = 4,
  cv.folds = 5, 
  bag.fraction = 0.5, 
  shrinkage = 0.01
)
summary(brtMan)
RFPred = predict(brtMan, data = test_Mangrove)
RMSE(test_Mangrove$loss, RFPred)
cor(train_Mangrove$loss, predict(brtMan, train_Mangrove, n.trees = 500)) #RSQ_Trn
cor(test_Mangrove$loss, predict(brtMan, test_Mangrove, n.trees = 500)) #RSQ_Tst
saveRDS(brtMan, "E:/Mangrove/data3/brtMan.rds")

plot(brtMan, i = "gdp")
plot(brtMan, i = "pop")
plot(brtMan, i = "hf")
plot(brtMan, i = "road")
plot(brtMan, i = "city")
plot(brtMan, i = "rice")
plot(brtMan, i = "sugarcane")
plot(brtMan, i = "coffee")
plot(brtMan, i = "coconut")
plot(brtMan, i = "rubber")
plot(brtMan, i = "olp")
plot(brtMan, i = "loss_prev")

varimp <- RFmodel$variable.importance
varimp <- data.frame(sort(varimp, decreasing = T))
varimp$Variables <- row.names(varimp)
row.names(varimp) <- NULL
colnames(varimp) <- c('Importance', 'Variables')
varimp$Importance_rel <- varimp$Importance / sum(varimp$Importance)
varimp <- varimp[,c('Variables', 'Importance', 'Importance_rel')]
write.csv(varimp, "E:/tidalflat/inputdata/importance/varimp.csv", row.names = F)
ggplot(data = varimp) +
  geom_bar(aes(x = Importance_rel, y = reorder(Variables, Importance_rel), fill = Importance_rel), stat = 'identity') +
  scale_fill_gradient(low = 'blue4', high = 'red4') +
  ggtitle("Variable importance (impurity)") +
  xlab("Reletive importance") +
  ylab("Variables") +
  labs(fill = NULL) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color = 'black', linewidth = 0.4, linetype = 1),
        panel.grid.major.x = element_line(color = 'black', linewidth = 0.2, linetype = 2),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 10))
ggsave("E:/Mangrove/data3/p55/varimp_impurity.tif", device = 'tiff', width = 3840 * 2, height = 2160 * 2, units = 'px', dpi = 600)

library(MASS)
gm1 = readRDS("C:/Users/tppji/Desktop/glmm1.rds")
gm2 = readRDS("C:/Users/tppji/Desktop/glmm2.rds")

summary(gm2)


bird = terra::vect("C:/Users/tppji/Desktop/BOTW_2022.2/bird.gdb")
