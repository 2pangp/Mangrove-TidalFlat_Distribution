setwd("E:/NUS/")

#library
library(ranger)
library(terra)
library(raster)

#sample points
samp07 = vect("E:/NUS/sample/2007sample.shp")
samp10 = vect("E:/NUS/sample/2010sample.shp")
samp13 = vect("E:/NUS/sample/2013sample.shp")
samp16 = vect("E:/NUS/sample/2016sample.shp")
samp19 = vect("E:/NUS/sample/2019sample.shp") 

#data extraction
CityAccess = raster("G:/data/Accessibility to cities/accessibility_to_cities_2015_v1.0/accessibility_to_cities_2015_v1.0.tif")
RoadAccess = raster("G:/data/road/road22.tif")
RoadAccess = projectRaster(RoadAccess, crs = crs(CityAccess))

CityAccess = calc(CityAccess, fun = function(x){ifelse((x < 0) & (x > -9999), 0, x)})
CityAccess = calc(CityAccess, fun = function(x){ifelse((x == -9999), NA, x)})
Access = data.frame(matrix(ncol = 5, nrow = 200000))
Access[,1] = raster::extract(CityAccess, samp07, buffer=2000, fun=mean, na.rm = T)
Access[,2] = raster::extract(CityAccess, samp10, buffer=2000, fun=mean, na.rm = T)
Access[,3] = raster::extract(CityAccess, samp13, buffer=2000, fun=mean, na.rm = T)
Access[,4] = raster::extract(CityAccess, samp16, buffer=2000, fun=mean, na.rm = T)
Access[,5] = raster::extract(CityAccess, samp19, buffer=2000, fun=mean, na.rm = T)
colnames(Access) = c("07","10","13","16","19")
write.csv(Access, "./data/CityAccess.csv", row.names = F)

Access1 = data.frame(matrix(ncol = 5, nrow = 200000))
Access1[,1] = raster::extract(RoadAccess, samp07, buffer=2000, fun=mean, na.rm = T)
Access1[,2] = raster::extract(RoadAccess, samp10, buffer=2000, fun=mean, na.rm = T)
Access1[,3] = raster::extract(RoadAccess, samp13, buffer=2000, fun=mean, na.rm = T)
Access1[,4] = raster::extract(RoadAccess, samp16, buffer=2000, fun=mean, na.rm = T)
Access1[,5] = raster::extract(RoadAccess, samp19, buffer=2000, fun=mean, na.rm = T)
colnames(Access1) = c("07","10","13","16","19")
write.csv(Access1, "./data/RoadAccess.csv", row.names = F)

#HFP
human07 = raster("G:/data/human/human05_07.tif")
human10 = raster("G:/data/human/human08_10.tif")
human13 = raster("G:/data/human/human11_13.tif")
human16 = raster("G:/data/human/human14_16.tif")
human19 = raster("G:/data/human/human17_19.tif")
human071 = projectRaster(human07, crs = crs(samp07))
human101 = projectRaster(human10, crs = crs(samp07))
human131 = projectRaster(human13, crs = crs(samp07))
human161 = projectRaster(human16, crs = crs(samp07))
human191 = projectRaster(human19, crs = crs(samp07))
human = data.frame(matrix(ncol = 5, nrow = 200000))
human[,1] = raster::extract(human071, samp07, buffer=2000, fun=mean, na.rm = T)
human[,2] = raster::extract(human101, samp10, buffer=2000, fun=mean, na.rm = T)
human[,3] = raster::extract(human131, samp13, buffer=2000, fun=mean, na.rm = T)
human[,4] = raster::extract(human161, samp16, buffer=2000, fun=mean, na.rm = T)
human[,5] = raster::extract(human191, samp19, buffer=2000, fun=mean, na.rm = T)

colnames(human) = c("07","10","13","16","19")
write.csv(human, "./data/human.csv", row.names = F)

Pop07 = raster("G:/data/Pop/Pop05_07.tif")
Pop10 = raster("G:/data/Pop/Pop08_10.tif")
Pop13 = raster("G:/data/Pop/Pop11_13.tif")
Pop16 = raster("G:/data/Pop/Pop14_16.tif")
Pop19 = raster("G:/data/Pop/Pop17_19.tif")
Pop = data.frame(matrix(ncol = 5, nrow = 200000))
Pop[,1] = raster::extract(Pop07, samp07, buffer=2000, fun=mean, na.rm = T)
Pop[,2] = raster::extract(Pop10, samp10, buffer=2000, fun=mean, na.rm = T)
Pop[,3] = raster::extract(Pop13, samp13, buffer=2000, fun=mean, na.rm = T)
Pop[,4] = raster::extract(Pop16, samp16, buffer=2000, fun=mean, na.rm = T)
Pop[,5] = raster::extract(Pop19, samp19, buffer=2000, fun=mean, na.rm = T)

colnames(Pop) = c("07","10","13","16","19")
write.csv(Pop, "./data/Pop.csv", row.names = F)

GDP07 = raster("G:/data/GDP_real1/GDP05_07.tif")
GDP10 = raster("G:/data/GDP_real1/GDP08_10.tif")
GDP13 = raster("G:/data/GDP_real1/GDP11_13.tif")
GDP16 = raster("G:/data/GDP_real1/GDP14_16.tif")
GDP19 = raster("G:/data/GDP_real1/GDP17_19.tif")

GDP = data.frame(matrix(ncol = 5, nrow = 200000))
GDP[,1] = raster::extract(GDP07, samp07, buffer=2000, fun=mean, na.rm = T)
GDP[,2] = raster::extract(GDP10, samp10, buffer=2000, fun=mean, na.rm = T)
GDP[,3] = raster::extract(GDP13, samp13, buffer=2000, fun=mean, na.rm = T)
GDP[,4] = raster::extract(GDP16, samp16, buffer=2000, fun=mean, na.rm = T)
GDP[,5] = raster::extract(GDP19, samp19, buffer=2000, fun=mean, na.rm = T)

colnames(GDP) = c("07","10","13","16","19")
write.csv(GDP, "./data/GDP.csv", row.names = F)

# Combining data and create a data frame
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

#Get latitude and longitude data
samp07[,c('lon', 'lat')] = geom(samp07, df = T)[,c('x', 'y')]
samp10[,c('lon', 'lat')] = geom(samp10, df = T)[,c('x', 'y')]
samp13[,c('lon', 'lat')] = geom(samp13, df = T)[,c('x', 'y')]
samp16[,c('lon', 'lat')] = geom(samp16, df = T)[,c('x', 'y')]
samp19[,c('lon', 'lat')] = geom(samp19, df = T)[,c('x', 'y')]
samp = rbind(samp07, samp10, samp13, samp16, samp19)
Data = cbind(Data, data.frame(samp[,c('lon', 'lat')]))
write.csv(Data, "./data/Data.csv", row.names = F)

#Get Protected area data
PA = raster("./DataMap/PAs.tif")
PA07 = raster::extract(PA, samp07, ID = F)
PA10 = raster::extract(PA, samp10, ID = F)
PA13 = raster::extract(PA, samp13, ID = F)
PA16 = raster::extract(PA, samp16, ID = F)
PA19 = raster::extract(PA, samp19, ID = F)
Data$PA = c(PA07,PA10,PA13,PA16,PA19)

#Get crop suitability data
Coconut = raster("G:/data/suitability/Coconut2010.tif")
Coffee = raster("G:/data/suitability/Coffee2010.tif")
Olp = raster("G:/data/suitability/olp2010.tif")
Rice = raster("G:/data/suitability/rice2010.tif")
Rubber= raster("G:/data/suitability/Rubber2010.tif")
Sugarcane= raster("G:/data/suitability/Sugarcane2010.tif")

Coconut07 = raster::extract(Coconut, samp07)
Coconut10 = raster::extract(Coconut, samp10)
Coconut13 = raster::extract(Coconut, samp13)
Coconut16 = raster::extract(Coconut, samp16)
Coconut19 = raster::extract(Coconut, samp19)
Data$coconut = c(Coconut07,Coconut10,Coconut13,Coconut16,Coconut19)
Coffee07 = raster::extract(Coffee, samp07)
Coffee10 = raster::extract(Coffee, samp10)
Coffee13 = raster::extract(Coffee, samp13)
Coffee16 = raster::extract(Coffee, samp16)
Coffee19 = raster::extract(Coffee, samp19)
Data$coffee = c(Coffee07,Coffee10,Coffee13,Coffee16,Coffee19)
Olp07 = raster::extract(Olp, samp07)
Olp10 = raster::extract(Olp, samp10)
Olp13 = raster::extract(Olp, samp13)
Olp16 = raster::extract(Olp, samp16)
Olp19 = raster::extract(Olp, samp19)
Data$olp = c(Olp07,Olp10,Olp13,Olp16,Olp19)
Rice07 = raster::extract(Rice, samp07)
Rice10 = raster::extract(Rice, samp10)
Rice13 = raster::extract(Rice, samp13)
Rice16 = raster::extract(Rice, samp16)
Rice19 = raster::extract(Rice, samp19)
Data$rice = c(Rice07,Rice10,Rice13,Rice16,Rice19)
Rubber07 = raster::extract(Rubber, samp07)
Rubber10 = raster::extract(Rubber, samp10)
Rubber13 = raster::extract(Rubber, samp13)
Rubber16 = raster::extract(Rubber, samp16)
Rubber19 = raster::extract(Rubber, samp19)
Data$rubber = c(Rubber07,Rubber10,Rubber13,Rubber16,Rubber19)
Sugarcane07 = raster::extract(Sugarcane, samp07)
Sugarcane10 = raster::extract(Sugarcane, samp10)
Sugarcane13 = raster::extract(Sugarcane, samp13)
Sugarcane16 = raster::extract(Sugarcane, samp16)
Sugarcane19 = raster::extract(Sugarcane, samp19)
Data$sugarcane = c(Sugarcane07,Sugarcane10,Sugarcane13,Sugarcane16,Sugarcane19)

#Get loss previous data
Loss04 = raster("E:/Mangrove/crop/loss_prev_04.tif")
Loss07 = raster("E:/Mangrove/crop/loss_prev_07.tif")
Loss10 = raster("E:/Mangrove/crop/loss_prev_10.tif")
Loss13 = raster("E:/Mangrove/crop/loss_prev_13.tif")
Loss16 = raster("E:/Mangrove/crop/loss_prev_16.tif")
Data$loss_prev = c(raster::extract(Loss04, samp07),
                   raster::extract(Loss07, samp10),
                   raster::extract(Loss10, samp13),
                   raster::extract(Loss13, samp16),
                   raster::extract(Loss16, samp19))
#Get loss data
Loss07 = raster("E:/NUS/data2/Loss2007.tif")
Loss10 = raster("E:/NUS/data2/Loss2010.tif")
Loss13 = raster("E:/NUS/data2/Loss2013.tif")
Loss16 = raster("E:/NUS/data2/Loss2016.tif")
Loss19 = raster("E:/NUS/data2/Loss2019.tif")
Data$loss = c(raster::extract(Loss07, samp07),
              raster::extract(Loss10, samp10),
              raster::extract(Loss13, samp13),
              raster::extract(Loss16, samp16),
              raster::extract(Loss19, samp19))
#Get ISO3 data
Data$ISO3 = c(ISO307$ISO3,ISO310$ISO3,ISO313$ISO3,ISO316$ISO3,ISO319$ISO3)

#Get River and coastline data
river = raster("E:/NUS/new data/river.tif")
coastline = raster("E:/NUS/new data/coastline.tif")
Data$river =  c(raster::extract(river, samp07),
                raster::extract(river, samp10),
                raster::extract(river, samp13),
                raster::extract(river, samp16),
                raster::extract(river, samp19))
Data$coastline =  c(raster::extract(coastline, samp07),
                raster::extract(coastline, samp10),
                raster::extract(coastline, samp13),
                raster::extract(coastline, samp16),
                raster::extract(coastline, samp19))

write.csv(Data,"./data/Data_raw.csv", row.names = F)

#Data sorting and processing
Data = read.csv("./data/Data_raw.csv")
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
write.csv(Data,"./data/Data1.csv", row.names = F)

Data1 = read.csv("./data/Data1.csv")
colnames(Data1)
Data1$gdp = scale(Data1$gdp)[,1]
Data1$hf = scale(Data1$hf)[,1]
Data1$pop = scale(Data1$pop)[,1]
Data1$city = scale(Data1$city)[,1]
Data1$road = scale(Data1$road)[,1]
Data1$coconut = scale(Data1$coconut)[,1]
Data1$coffee = scale(Data1$coffee)[,1]
Data1$olp = scale(Data1$olp)[,1]
Data1$rice = scale(Data1$rice)[,1]
Data1$rubber = scale(Data1$rubber)[,1]
Data1$sugarcane = scale(Data1$sugarcane)[,1]
Data1$coconut_prod = scale(Data1$coconut_prod)[,1]
Data1$coffee_prod = scale(Data1$coffee_prod)[,1]
Data1$olp_prod = scale(Data1$olp_prod)[,1]
Data1$rice_prod = scale(Data1$rice_prod)[,1]
Data1$rubber_prod = scale(Data1$rubber_prod)[,1]
Data1$sugarcane_prod = scale(Data1$sugarcane_prod)[,1]
Data1$river = scale(Data1$river)[,1]
Data1$coastline = scale(Data1$coastline)[,1]

write.csv(Data1,"./data/Mangrove.csv", row.names = F)


# Modelling
library(tidyverse)
library(terra)
library(raster)
library(randomForest)
library(ranger)
library(gbm)
library(caret)
Mangrove = read.csv("./data/Mangrove.csv")
Mangrove = Mangrove[c("loss","gdp", "hf", "pop", "city", "road","PA", "coconut", "coffee", "olp", "rice", "rubber", "sugarcane","loss_prev")]
#colnames(Mangrove) = c("loss","gdp", "hf", "pop", "city", "road","PA", "coconut", "coffee", "olp", "rice", "rubber", "sugarcane","loss_prev")
#01 data split
head(Mangrove)
set.seed(10)
randidx = sample(nrow(Mangrove), 0.8*nrow(Mangrove))
train_Mangrove= Mangrove[randidx,]
test_Mangrove = Mangrove[-randidx,]
str(Mangrove)
mod_f = as.formula(loss ~ gdp  + hf  + pop  + city + road + PA + coconut + coffee + olp + 
                     rice + rubber + sugarcane+loss_prev+river+coastline)
RFmodel = ranger(mod_f,data= train_Mangrove, num.trees = 300, mtry = 15/3,importance = 'impurity', verbose = T)
saveRDS(RFmodel, "./data/RF1.rds")

# Variable importance
importance(RFmodel)
pdp = partial(RFmodel,
              train = train_Mangrove,
              pred.var = 'loss_prev')
line(pdp$loss_prev, pdp$yhat, type = 'l')


# Testing
RFmodel = readRDS("./data/RF1.rds")
RFPred = predict(RFmodel, data = test_Mangrove)
test_Mangrove$pred <- RFPred$predictions
RMSE(test_Mangrove$loss, RFPred$predictions) / (max(test_Mangrove$loss) - min(test_Mangrove$loss))

test_mod <- lm(test_Mangrove, formula = pred ~ loss - 1)
summary(test_mod)

cor(train_Mangrove$loss, predict(RFmodel, train_Mangrove, n.trees = 300)$predictions) #RSQ_Trn
cor(test_Mangrove$loss, predict(RFmodel, test_Mangrove, n.trees = 300)$predictions) #RSQ_Tst


# Prediction
RFmodel = readRDS("./data/RF1.rds")
pred_list = rast("./predict/2030ssp1pred_list.tif")
pred_list$city = lapp(pred_list$city, function(x){ifelse(x < 0, 0, x)})
layerIdx = c(1,2,3,4,5,7,8,9,10,11,12)
for (layer in layerIdx)
{
  pred_list[[layer]] = scale(pred_list[[layer]])
  print(layer)
}

pred_set = list()
for (layer in 1:13)
{
  pred_set[[layer]] = raster(pred_list[[layer]])
  print(layer)
}
pred_set = stack(pred_set)
writeRaster(pred_set, "./predict/2030ssp1pred_set.tif")
plot(pred_set)

pred_set = rast("./predict/2030ssp1pred_set.tif")
names(pred_set) = c("gdp", "hf", "pop", "city", "road", "PA", "coconut", "coffee", "olp", "rice", 
                    "rubber", "sugarcane","loss_prev","coconut_prod","coffee_prod","olp_prod", "rice_prod",
                    "rubber_prod", "sugarcane_prod","river","coastline")
PredRas = predict(pred_set,RFmodel,na.rm=T)
writeRaster(PredRas, "./predict/2030ssp1_pred.tif")

plot(PredRas)
summary(RFmodel)


#BRT modelling
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
saveRDS(brtMan, "./data/brtMan.rds")

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


