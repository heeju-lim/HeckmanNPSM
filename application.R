library(AER)
data("PSID1976")

cc         <- ifelse(PSID1976$participation=="yes", 1, 0) # cc == 0  missing data
wage      <- PSID1976$wage
hwage      <- PSID1976$hwage
youngkids  <- PSID1976$youngkids
tax        <- PSID1976$tax
feducation <- PSID1976$feducation
education  <- PSID1976$education
city       <- PSID1976$city

# Transformado as variáveis "factors" em quantitativas
city     <- ifelse(city=="yes", 1, 0)     # yes == 1 and no == 0

x1= hwage
x2= youngkids
x3= tax
x4= feducation
x5= education
x6= city

#library(HeckmanEM)

x <- cbind(1,x5,x6)
w <- cbind(1,x1,x2,x3,x4,x5,x6)
y <- ifelse(wage>0, log(wage), 0)


n<-length(cc)
npsm<- HeckmanNPSM(y, x, w, cc)

