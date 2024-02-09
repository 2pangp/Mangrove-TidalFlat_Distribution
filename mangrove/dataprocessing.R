#library
library(rgdal)
library(sp)
library(sf)
library(terra)
library(raster)

setwd("E:/Mangrove")

samp04 = shapefile("E:/Nus/dian/2004sample.shp")
samp07 = shapefile("E:/Nus/dian/2007sample.shp")
samp10 = shapefile("E:/Nus/dian/2010sample.shp")
samp13 = shapefile("E:/Nus/dian/2013sample.shp")
samp16 = shapefile("E:/Nus/dian/2016sample.shp") 
samp19 = shapefile("E:/Nus/dian/2019sample.shp") 



#数据提取
CityAccess = raster("G:/data/Accessibility to cities/accessibility_to_cities_2015_v1.0/accessibility_to_cities_2015_v1.0.tif")
RoadAccess = raster("G:/data/road/road22.tif")
RoadAccess = projectRaster(RoadAccess, crs = crs(CityAccess))

CityAccess = calc(CityAccess, fun = function(x){ifelse((x < 0) & (x > -9999), 0, x)})
CityAccess = calc(CityAccess, fun = function(x){ifelse((x == -9999), NA, x)})
Access = data.frame(matrix(ncol = 6, nrow = 100000))
Access[,1] = raster::extract(CityAccess, samp04, buffer=2000, fun=mean, na.rm = T)
Access[,2] = raster::extract(CityAccess, samp07, buffer=2000, fun=mean, na.rm = T)
Access[,3] = raster::extract(CityAccess, samp10, buffer=2000, fun=mean, na.rm = T)
Access[,4] = raster::extract(CityAccess, samp13, buffer=2000, fun=mean, na.rm = T)
Access[,5] = raster::extract(CityAccess, samp16, buffer=2000, fun=mean, na.rm = T)
Access[,6] = raster::extract(CityAccess, samp16, buffer=2000, fun=mean, na.rm = T)

colnames(Access) = c("04","07","10","13","16")
write.csv(Access, "./data2/CityAccess.csv", row.names = F)

Access = data.frame(matrix(ncol = 5, nrow = 100000))
Access[,1] = raster::extract(RoadAccess, samp04, buffer=2000, fun=mean, na.rm = T)
Access[,2] = raster::extract(RoadAccess, samp07, buffer=2000, fun=mean, na.rm = T)
Access[,3] = raster::extract(RoadAccess, samp10, buffer=2000, fun=mean, na.rm = T)
Access[,4] = raster::extract(RoadAccess, samp13, buffer=2000, fun=mean, na.rm = T)
Access[,5] = raster::extract(RoadAccess, samp16, buffer=2000, fun=mean, na.rm = T)

colnames(Access) = c("04","07","10","13","16")
write.csv(Access, "./data2/RoadAccess.csv", row.names = F)

#HFP
humanFiles = list.files("G:/data/human")
Data = matrix(ncol = 5, nrow = 100000)
for (i in 1:length(humanFiles))
{
  dir = paste0("G:/data/human/", humanFiles[i])
  rastTemp = raster(dir)
  rastTemp = projectRaster(rastTemp, crs = crs(samp04))
  Data[,i] = raster::extract(rastTemp, samp07, buffer=2000, fun=mean, na.rm = T)
  print(i)
}
Data =as.data.frame(Data)
colnames(Data) = substr(humanFiles, 1, 8)
write.csv(Data, "./data2/human07.csv", row.names = F)
sum(is.na(Data))

humanchangeFiles = list.files("G:/data/humanchange")
HFP = matrix(ncol = 5, nrow = 100000)
for (i in 1:length(humanchangeFiles))
{
  dir = paste0("G:/data/humanchange/", humanchangeFiles[i])
  rastTemp = raster(dir)
  rastTemp = projectRaster(rastTemp, crs = crs(Man_samp))
  HFP[,i] = raster::extract(rastTemp, Man_samp, buffer=2000, fun=mean, na.rm = T)
  print(i)
}
HFP =as.data.frame(HFP)
colnames(HFP) = substr(humanchangeFiles, 1, 8)
write.csv(HFP, "./data1/humanchange.csv", row.names = F)
sum(is.na(HFP))


#POP
PopFiles = list.files("G:/data/Pop")
PopData = matrix(ncol = 5, nrow = 100000)
for (i in 1:length(PopFiles))
{
  dir = paste0("G:/data/Pop/", PopFiles[i])
  Pop = raster(dir)
  PopData[,i] = raster::extract(Pop, samp04, buffer=2000, fun=mean, na.rm = T)
  print(i)
}
PopData =as.data.frame(PopData)
colnames(PopData) = substr(PopFiles, 1, nchar(PopFiles)-4)
write.csv(PopData, "./data2/Pop2004.csv", row.names = F)
sum(is.na(PopData))



PopchangeFiles = list.files("G:/data/Popchange")
PopchangeData = matrix(ncol = 5, nrow = 100000)
for (i in 1:length(PopchangeFiles))
{
  dir = paste0("G:/data/Popchange/",PopchangeFiles[i])
  Popchange = raster(dir)
  PopchangeData[,i] = raster::extract(Popchange, Man_samp, buffer=2000, fun=mean, na.rm = T)
  print(i)
}
PopchangeData =as.data.frame(PopchangeData)
colnames(PopchangeData) = substr(PopchangeFiles, 1, nchar(PopchangeFiles)-4)
write.csv(PopchangeData, "./data1/Popchange.csv", row.names = F)
sum(is.na(PopchangeData))


library(ncdf4)
#GDP
GDP = list()
for (layer in 1:26)
{
  layerTemp = raster("G:/data/gdp_Kummu/GDP_per_capita_PPP_1990_2015_v2.nc", band = layer)
  GDP[[layer]] = layerTemp
  print(layer)
}
GDP =  stack(GDP)

#99_04
GDP99_04 = GDP[[10:15]]
GDP99_04_mean = mean(GDP99_04)
GDP99_04_change = (GDP99_04[[5]] - GDP99_04[[1]]) / GDP99_04[[1]]
#05_07
GDP05_07 = GDP[[16:18]]
GDP05_07_mean = mean(GDP05_07)
GDP05_07_change = (GDP05_07[[3]] - GDP05_07[[1]]) / GDP05_07[[1]]
#08_10
GDP08_10 = GDP[[19:21]]
GDP08_10_mean = mean(GDP08_10)
GDP08_10_change = (GDP08_10[[2]] - GDP08_10[[1]]) / GDP08_10[[1]]
#11_13
GDP11_13 = GDP[[22:24]]
GDP11_13_mean = mean(GDP11_13)
GDP11_13_change = (GDP11_13[[3]] - GDP11_13[[1]]) / GDP11_13[[1]]
#14_15
GDP14_15 = GDP[[25:26]]
GDP14_15_mean = mean(GDP14_15)
GDP14_15_change = (GDP14_15[[2]] - GDP14_15[[1]]) / GDP14_15[[1]]

GDPData = matrix(ncol = 5, nrow = 100000)
GDPData[,1] = raster::extract(GDP99_04_mean, samp16, buffer=2000, fun=mean, na.rm = T)
GDPData[,2] = raster::extract(GDP05_07_mean, samp16, buffer=2000, fun=mean, na.rm = T)
GDPData[,3] = raster::extract(GDP08_10_mean, samp16, buffer=2000, fun=mean, na.rm = T)
GDPData[,4] = raster::extract(GDP11_13_mean, samp16, buffer=2000, fun=mean, na.rm = T)
GDPData[,5] = raster::extract(GDP14_15_mean, samp16, buffer=2000, fun=mean, na.rm = T)
colnames(GDPData) = c("GDP99_04", "GDP05_07", "GDP08_10", "GDP11_13", "GDP14_15")
write.csv(GDPData, "./data2/GDP2016.csv", row.names = F)
sum(is.na(GDPData))

GDPData = matrix(ncol = 5, nrow = 100000)
GDPData[,1] = raster::extract(GDP99_04_change, Man_samp, buffer=2000, fun=mean, na.rm = T)
GDPData[,2] = raster::extract(GDP05_07_change, Man_samp, buffer=2000, fun=mean, na.rm = T)
GDPData[,3] = raster::extract(GDP08_10_change, Man_samp, buffer=2000, fun=mean, na.rm = T)
GDPData[,4] = raster::extract(GDP11_13_change, Man_samp, buffer=2000, fun=mean, na.rm = T)
GDPData[,5] = raster::extract(GDP14_15_change, Man_samp, buffer=2000, fun=mean, na.rm = T)
colnames(GDPData) = c("GDP99_04", "GDP05_07", "GDP08_10", "GDP11_13", "GDP14_15")
write.csv(GDPData, "./data1/GDPchange.csv", row.names = F)
sum(is.na(GDPData))


#数据合并创建表格
CityAccess = read.csv("./data2/CityAccess.csv")
RoadAccess = read.csv("./data2/RoadAccess.csv")
Access = read.csv("./data2/2019Access.csv")
human04 = read.csv("./data2/human04.csv")
human07 = read.csv("./data2/human07.csv")
human10 = read.csv("./data2/human10.csv")
human13 = read.csv("./data2/human13.csv")
human16 = read.csv("./data2/human16.csv")
human19 = read.csv("./data2/human19.csv")
GDP = read.csv("./data2/GDP.csv")
Pop04 = read.csv("./data2/Pop2004.csv")
Pop07 = read.csv("./data2/Pop2007.csv")
Pop10 = read.csv("./data2/Pop2010.csv")
Pop13 = read.csv("./data2/Pop2013.csv")
Pop16 = read.csv("./data2/Pop2016.csv")
Pop19 = read.csv("./data2/Pop2019.csv")

dim(Pop04)

Data = as.data.frame(matrix(ncol = 6, nrow = 600000))
colnames(Data) = c("year", "gdp", "hf", "pop", "city", "road")
Data$year = c(rep(2004, 100000), rep(2007, 100000), rep(2010, 100000), rep(2013, 100000), rep(2016, 100000), rep(2019, 100000))
Data$gdp = c(GDP$X1,GDP$X2,GDP$X3,GDP$X4,GDP$X5,GDP$X6)
Data$hf = c(human04$human99_, human07$human05_, human10$human08_,human13$human11_, human16$human14_,human19$X2019)
Data$pop = c(Pop04$Pop99_04, Pop07$Pop05_07, Pop10$Pop08_10, Pop13$Pop11_13, Pop16$Pop14_16, Pop19$X2019)
Data$city = c(CityAccess$X04,CityAccess$X07,CityAccess$X10,CityAccess$X13,CityAccess$X16,Access$city)
Data$road =  c(RoadAccess$X04,RoadAccess$X07,RoadAccess$X10,RoadAccess$X13,RoadAccess$X16,Access$road)

Man_samp_df = as.data.frame(samp04,samp07,samp10,samp13,samp16,samp19)
colnames(Man_samp_df) = c("y04","y07", "y11", "y13", "y16", "y19")
Data$loss = c(samp04$CID, samp07$CID,samp10$CID, samp13$CID,samp16$CID,samp19$CID)
table(Data$loss)


PA = raster("./DataMap/PAs.tif")
PA04 = raster::extract(PA, samp04, ID = F)
PA07 = raster::extract(PA, samp07, ID = F)
PA10 = raster::extract(PA, samp10, ID = F)
PA13 = raster::extract(PA, samp13, ID = F)
PA16 = raster::extract(PA, samp16, ID = F)
PA19 = raster::extract(PA, samp19, ID = F)
Data$PA = c(PA04,PA07,PA10,PA13,PA16,PA19)


Coconut = raster("G:/data/suitability/Coconut2010.tif")
Coffee = raster("G:/data/suitability/Coffee2010.tif")
Olp = raster("G:/data/suitability/olp2010.tif")
Rice = raster("G:/data/suitability/rice2010.tif")
Rubber= raster("G:/data/suitability/Rubber2010.tif")
Sugarcane= raster("G:/data/suitability/Sugarcane2010.tif")

Coconut04 = extract(Coconut, samp04)
Coconut07 = extract(Coconut, samp07)
Coconut10 = extract(Coconut, samp10)
Coconut13 = extract(Coconut, samp13)
Coconut16 = extract(Coconut, samp16)
Coconut19 = extract(Coconut, samp19)
Data$coconut = c(Coconut04,Coconut07,Coconut10,Coconut13,Coconut16,Coconut19)
Coffee04 = extract(Coffee, samp04)
Coffee07 = extract(Coffee, samp07)
Coffee10 = extract(Coffee, samp10)
Coffee13 = extract(Coffee, samp13)
Coffee16 = extract(Coffee, samp16)
Coffee19 = extract(Coffee, samp19)
Data$Coffee = c(Coffee04,Coffee07,Coffee10,Coffee13,Coffee16,Coffee19)
Olp04 = extract(Olp, samp04)
Olp07 = extract(Olp, samp07)
Olp10 = extract(Olp, samp10)
Olp13 = extract(Olp, samp13)
Olp16 = extract(Olp, samp16)
Olp19 = extract(Olp, samp19)
Data$Olp = c(Olp04,Olp07,Olp10,Olp13,Olp16,Olp19)
Rice04 = extract(Rice, samp04)
Rice07 = extract(Rice, samp07)
Rice10 = extract(Rice, samp10)
Rice13 = extract(Rice, samp13)
Rice16 = extract(Rice, samp16)
Rice19 = extract(Rice, samp19)
Data$Rice = c(Rice04,Rice07,Rice10,Rice13,Rice16,Rice19)
Rubber04 = extract(Rubber, samp04)
Rubber07 = extract(Rubber, samp07)
Rubber10 = extract(Rubber, samp10)
Rubber13 = extract(Rubber, samp13)
Rubber16 = extract(Rubber, samp16)
Rubber19 = extract(Rubber, samp19)
Data$Rubber = c(Rubber04,Rubber07,Rubber10,Rubber13,Rubber16,Rubber19)
Sugarcane04 = extract(Sugarcane, samp04)
Sugarcane07 = extract(Sugarcane, samp07)
Sugarcane10 = extract(Sugarcane, samp10)
Sugarcane13 = extract(Sugarcane, samp13)
Sugarcane16 = extract(Sugarcane, samp16)
Sugarcane19 = extract(Sugarcane, samp19)
Data$Sugarcane = c(Sugarcane04,Sugarcane07,Sugarcane10,Sugarcane13,Sugarcane16,Sugarcane19)

write.csv(Data, "./data2/data_raw.csv", row.names = F)


#数据异常值处理
Data = read.csv("./data2/data_raw.csv")
Data = na.omit(Data)
IQR_gdp = IQR(Data$gdp, na.rm = T)
IQR_hf = IQR(Data$hf, na.rm = T)
IQR_pop = IQR(Data$pop, na.rm = T)
IQR_city = IQR(Data$city, na.rm = T)
IQR_road = IQR(Data$road, na.rm = T)

Q_gdp = quantile(Data$gdp, probs = c(0.25, 0.75), na.rm = T)
Q_hf = quantile(Data$hf, probs = c(0.25, 0.75), na.rm = T)
Q_pop = quantile(Data$pop, probs = c(0.25, 0.75), na.rm = T)
Q_city = quantile(Data$city, probs = c(0.25, 0.75), na.rm = T)
Q_road = quantile(Data$road, probs = c(0.25, 0.75), na.rm = T)

Data = Data[ - which((Data$gdp > Q_gdp[2] + 5 * IQR_gdp) | (Data$gdp < Q_gdp[1] - 5 * IQR_gdp)),]
#Data = Data[ - which((Data$gdp_var > Q_gdp_var[2] + 9 * IQR_gdp_var) | (Data$gdp_var < Q_gdp_var[1] - 9 * IQR_gdp_var)),]
#Data = Data[ - which((Data$hf > Q_hf[2] + 5 * IQR_hf) | (Data$hf < Q_hf[1] - 5 * IQR_hf)),]
#Data = Data[ - which((Data$hf_var > Q_hf_var[2] + 5 * IQR_hf_var) | (Data$hf_var < Q_hf_var[1] - 5 * IQR_hf_var)),]
Data = Data[ - which((Data$pop > Q_pop[2] + 5 * IQR_pop) | (Data$pop < Q_pop[1] - 5 * IQR_pop)),]
#Data = Data[ - which((Data$pop_var > Q_pop_var[2] + 5 * IQR_pop_var) | (Data$pop_var < Q_pop_var[1] - 5 * IQR_pop_var)),]
#Data = Data[ - which((Data$city > Q_city[2] + 1.5 * IQR_city) | (Data$city < Q_city[1] - 1.5 * IQR_city)),]
#Data = Data[ - which((Data$road > Q_road[2] + 1.5 * IQR_road) | (Data$road < Q_road[1] - 1.5 * IQR_road)),]
write.csv(Data,"./data2/Data1.csv" )

#标准化
Data1 = read.csv("./data2/Data1.csv")
Data1$gdp = scale(Data1$gdp)[,1]
Data1$hf = scale(Data1$hf)[,1]
Data1$pop = scale(Data1$pop)[,1]
Data1$city = scale(Data1$city)[,1]
Data1$road = scale(Data1$road)[,1]
Data1$coconut = scale(Data1$coconut)[,1]
Data1$Coffee = scale(Data1$Coffee)[,1]
Data1$Olp = scale(Data1$Olp)[,1]
Data1$Rice = scale(Data1$Rice)[,1]
Data1$Rubber = scale(Data1$Rubber)[,1]
Data1$Sugarcane = scale(Data1$Sugarcane)[,1]
write.csv(Data1,"./data2/Data2.csv")


#模型构建
library(tidyverse)
library(terra)
library(raster)
library(randomForest)


Mangrove = read.csv("E:/Mangrove/data3/p01/Mangrove.csv")
Mangrove$loss = as.factor(Mangrove$loss)
Mangrove$PA = as.factor(Mangrove$PA)
Mangrove = Mangrove[c("loss","gdp", "hf", "pop", "city", "road","PA", "coconut", "Coffee", "Olp", "Rice", "Rubber", "Sugarcane")]
colnames(Mangrove) = c("loss","gdp", "hf", "pop", "city", "road","PA", "coconut", "coffee", "olp", "rice", "rubber", "sugarcane")

#01 data split
head(Mangrove)
set.seed(10)
randidx = sample(nrow(Mangrove), 0.8*nrow(Mangrove))
train_Mangrove= Mangrove[randidx,]
test_Mangrove = Mangrove[-randidx,]
#02采样
datarecipe <- recipe(loss ~ ., train_Mangrove) %>% 
  step_downsample(loss, under_ratio = 1.5) %>%
  #step_rose(loss, over_ratio = 1) %>%
  prep()
datarecipe

train_Mangrove <- bake(datarecipe, new_data = NULL) %>%
  select(loss, everything())
test_Mangrove <- bake(datarecipe, new_data = test_Mangrove) %>%
  select(loss, everything())

table(train_Mangrove$loss)
table(test_Mangrove$loss)


mod_f = as.formula(loss ~ gdp  + hf  + pop  + city + road + PA + coconut + coffee + olp + rice + rubber + sugarcane+river+coastline)
RFmodel = randomForest(mod_f,data=train_Mangrove, ntrees = 300, mtry = 15/3, importance = T, type = "classification")
saveRDS(RFmodel, "E:/Mangrove/data3/p01/RF1.rds")
predRF = predict(RFmodel, newdata = test_Mangrove)
mean(predRF == test_Mangrove$loss)

#模型预测
RFmodel = readRDS("E:/Mangrove/data2/RF1.rds")
pred_list = rast("E:/Mangrove/data2/SSP5pred_list.tif")
pred_list$city = lapp(pred_list$city, function(x){ifelse(x < 0, 0, x)})
layerIdx = c(1,2,3,4,5,7,8,10,11,12)
for (layer in layerIdx)
{
  pred_list[[layer]] = scale(pred_list[[layer]])
  print(layer)
}

pred_set = list()
for (layer in 1:12)
{
  pred_set[[layer]] = raster(pred_list[[layer]])
  print(layer)
}
pred_set = stack(pred_set)
writeRaster(pred_set, "E:/Mangrove/data2/SSP5pred_set.tif")
summary(pred_set)

pred_set = rast("E:/Mangrove/data3/pred_list.tif")
pred_set = stack(pred_set)
names(pred_set) = c("gdp", "hf", "pop", "city", "road", "PA", "coconut", "coffee", "olp", "rice", "rubber", "sugarcane")
PredRas = predict(pred_set,RFmodel,na.rm=T)
plot(PredRas)
writeRaster(PredRas, "E:/Mangrove/data3/mangrove16_19.tif")





Loss04 = vect("E:/NUS/mian/2004mian3.shp")
Loss07 = vect("E:/NUS/mian/2007mian3.shp")
Loss10 = vect("E:/NUS/mian/2010mian3.shp")
Loss13 = vect("E:/NUS/mian/2013mian3.shp")
Loss16 = vect("E:/NUS/mian/2016mian3.shp")
Loss19 = vect("E:/NUS/mian/2019mian3.shp")
PredRas = rast("E:/Mangrove/data2/mangrove16_19_pred.tif")

Loss_total = rbind(Loss04, Loss07, Loss10, Loss13, Loss16,Loss19)
R = rast(ext(PredRas1), ncol = ncol(PredRas1), nrow = nrow(PredRas1))
R = project(R, crs(PredRas1))
Loss19$val1 = 1
LossRas19 = rasterize(Loss19, R, field = "val1", background = NA, touches = T)
writeRaster(LossRas19 , "E:/Mangrove/data2/LossRas.tif", overwrite = T)
LossRas19= rast("E:/Mangrove/data2/LossRas.tif")
global(LossRas19, fun = "sum", na.rm = T)
global(PredRas, fun = "sum", na.rm = T)

test = function(x, y)
{
  ifelse(x == y, 1, 0)
}
test_list = lapp(sds(list(PredRas, LossRas19)), fun = test)
global(test_list, fun = "sum", na.rm = T)


plot(LossRas19)

plot(PredRas)

Mangrove16_19 = mask(PredRas,Loss_total, inverse = T)
plot(Mangrove16_19)
writeRaster(Mangrove16_19 , "E:/Mangrove/data2/Mangrove16_19 .tif", overwrite = T)



