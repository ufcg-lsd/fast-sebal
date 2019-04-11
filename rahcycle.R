library(R.utils)
library(raster)
library(rgdal)
library(maptools)
library(ncdf4)
library(sp)
library(snow)
library(here)

Rn<-raster("Testes/Teste04/Rn.tif")
TS<-raster("Testes/Teste04/TS.tif")
NDVI<-raster("Testes/Teste04/NDVI.tif")
G<-raster("Testes/Teste04/G.tif")
alb<-raster("Testes/Teste04/alb.tif")
fic.sw <- "input/station.csv"
table.sw <- (read.csv(fic.sw, sep=";", header=FALSE, stringsAsFactors=FALSE))

k <- 0.41		# Von Karman
g <- 9.81		# Gravity
rho <- 1.15		# Air density
cp <- 1004		# Specific heat of air
Gsc <- 0.082		# Solar constant (0.0820 MJ m-2 min-1)

x<-3                                    # Wind speed sensor Height (meters)
hc<-0.2                                 # Vegetation height (meters)
Lat<-table.sw$V4[1]    # Station Latitude
Long<-table.sw$V5[1]   # Station Longitude

    # Surface roughness parameters in station
azom <- -3              #Parameter for the Zom image
bzom <- 6.47    #Parameter for the Zom image
F_int <- 0.16

zom.est <- hc*0.12

	# Friction velocity at the station (ustar.est)
ustar.est <- k*table.sw$V6[2]/log((x)/zom.est)
	
	# Velocity 200 meters
u200 <- ustar.est/k*log(200/zom.est)
	
	# Zom for all pixels
zom <- NDVI
zom[] <- exp(azom+bzom*NDVI[])
print("Passou do zom")
	
	# Initial values
ustar<-NDVI
ustar[]<-k*u200/(log(200/zom[]))		# Friction velocity for all pixels #RASTER - VETOR 
	
ustar_before <- NDVI
ustar_before[] <- k*u200/(log(200/zom[])) #FIXME:
print("Passou ustar")
	
rah<-NDVI
rah[]<-(log(2/0.1))/(ustar[]*k) 		# Aerodynamic resistance for all pixels #RASTER - VETOR
	
rah_before <- NDVI
rah_before[] <-(log(2/0.1))/(ustar[]*k)  #FIXME:
print("Passou rah")

hot.pixel.rn <- Rn[1408, 1440]
hot.pixel.g <- G[1408, 1440]
hot.pixel.ts <- TS[1408, 1440]
cold.pixel.ts <- TS[4289, 5765]
print(c(hot.pixel.rn, hot.pixel.g, hot.pixel.ts, cold.pixel.ts))

H.hot<-hot.pixel.rn-hot.pixel.g 
value.pixel.rah<-rah[1408, 1440]

i<-1
Erro<-TRUE
print("Before rah correct")
	#print(proc.time())
	# Beginning of the cycle stability
while(Erro){

    print(paste("Loop", i))

    # print("Ustar before loop")
    # print(ustar[])

    # print("Rah before loop")
    # print(rah[])

    rah.hot.0<-value.pixel.rah[i] # Value
    
    # Hot and cold pixels      
    dt.hot<-H.hot*rah.hot.0/(rho*cp) # Value
    b<-dt.hot/(hot.pixel.ts-cold.pixel.ts) # Value
    a<- -b*(cold.pixel.ts-273.15) # Value
    
	# All pixels
    H<-rho*cp*(a+b*(TS[]-273.15))/rah[] # Changed from Raster to Vector
    L<- -1*((rho*cp*ustar[]^3*TS[])/(k*g*H)) # Changed from Raster to Vector
    y_0.1<-(1-16*0.1/L)^0.25 # Changed from Raster to Vector
    y_2<-(1-16*2/L)^0.25 # Changed from Raster to Vector
    x200<-(1-16*200/L)^0.25 # Changed from Raster to Vector
    psi_0.1<-2*log((1+y_0.1^2)/2) # Changed from Raster to Vector
    psi_0.1[L>0 &!is.na(L)]<--5*(0.1/L[L>0 &!is.na(L)]) # Changed from Raster to Vector
    psi_2<-2*log((1+y_2^2)/2)  # Changed from Raster to Vector
    psi_2[L>0 &!is.na(L) ]<--5*(2/L[L>0 &!is.na(L)]) # Changed from Raster to Vector
    psi_200<-2*log((1+x200)/2)+log((1+x200^2)/2)-2*atan(x200)+0.5*pi # Changed from Raster to Vector
    psi_200[L>0 &!is.na(L) ]<--5*(2/L[(L>0 &!is.na(L))]) # Changed from Raster to Vector
    ustar<-k*u200/(log(200/zom[])-psi_200) # Changed from Raster to Vector # Friction velocity for all pixels
    rah<-NDVI
    rah[]<-(log(2/0.1)-psi_2+psi_0.1)/(ustar*k) # Changed from Raster to Vector # Aerodynamic resistency for all pixels
    rah.hot <- rah[1408, 1440]
    value.pixel.rah<-c(value.pixel.rah,rah.hot) # Value
    
    # print("L")
    # print(L[])
    # print("")

    # print("y_0.1")
    # print(y_0.1[])
    # print("")

    # print("y_2")
    # print(y_2[])
    # print("")

    # print("x_200")
    # print(x200[])
    # print("")

    # print("psi_0.1")
    # print(psi_0.1[])
    # print("")

    # print("psi_2")
    # print(psi_2[])
    # print("")

    # print("psi_200")
    # print(psi_200[])
    # print("")

    print("Rahs do hot pixel")
    print(value.pixel.rah)

    print("Erro")
    print(abs(1-rah.hot.0/rah.hot))
    Erro<-(abs(1-rah.hot.0/rah.hot)>=0.05)
    i<-i+1
}
	
	#print(proc.time())
	# End sensible heat flux (H)
	
dt.hot<-H.hot*rah.hot.0/(rho*cp) # Value
b<-dt.hot/(hot.pixel.ts-cold.pixel.ts) # Value
a<- -b*(cold.pixel.ts-273.15) # Value                      

print("H HOT")
print(H.hot)
print("Rah hot")
print(rah.hot.0)
print("DT HOT")
print(dt.hot)
print("B")
print(b)
print("A")
print(a)

	#print(proc.time())
	
	# All pixels
H <-rho*cp*(a+b*(TS[]-273.15))/rah[] # Vector 
H[(H>(Rn[]-G[]) &!is.na(H))]<-(Rn[]-G[])[(H>(Rn[]-G[]) &!is.na(H))] # Vector
	
ustar_final <- NDVI
ustar_final[] <- ustar


H_final <- NDVI
H_final[] <- H


	# Instant latent heat flux (LE)
LE<-Rn[]-G[]-H
LE_final <- NDVI
LE_final[] <- Rn[]-G[]-H

Dia.juliano <- 3

d_sun_earth <- 0.98330
	# Upscalling temporal
dr<-(1/d_sun_earth)^2 		# Inverse square of the distance on Earth-SOL
sigma<-0.409*sin(((2*pi/365)*Dia.juliano)-1.39) # Declination Solar (rad)
phi<-(pi/180)*Lat 								# Solar latitude in degrees
omegas<-acos(-tan(phi)*tan(sigma)) 				# Angle Time for sunsets (rad)
Ra24h<-(((24*60/pi)*Gsc*dr)*(omegas*sin(phi)*
	        sin(sigma)+cos(phi)*cos(sigma)*sin(omegas)))*(1000000/86400)
	
	#print(proc.time())
	
	# Short wave radiation incident in 24 hours (Rs24h)
Rs24h<-F_int*sqrt(max(table.sw$V7[])-min(table.sw$V7[]))*Ra24h
	
FL<-110                                
Rn24h_dB<-(1-alb[])*Rs24h-FL*Rs24h/Ra24h		# Method of Bruin #VETOR
Rn24h_dB_final <- NDVI
Rn24h_dB_final[] <-(1-alb[])*Rs24h-FL*Rs24h/Ra24h
	# Evapotranspiration fraction Bastiaanssen
EF<-NDVI
EF[]<-LE/(Rn[]-G[])
	
	# Sensible heat flux 24 hours (H24h), nao foi usado
	# H24h_dB<-(1-EF[])*Rn24h_dB
	
	# Latent Heat Flux 24 hours (LE24h)
LE24h_dB<-EF[]*Rn24h_dB
LE24h_dB_final <- NDVI
LE24h_dB_final[] <- EF[]*Rn24h_dB
	
	# Evapotranspiration 24 hours (ET24h)
ET24h_dB<-NDVI
ET24h_dB[]<-LE24h_dB*86400/((2.501-0.00236* (max(table.sw$V7[])+min(table.sw$V7[]))/2)*10^6)
	#print(proc.time())
	
# print("ZOM")
# print(zom[])
# print("USTAR")
# print(ustar_final[])
# print("RAH")
# print(rah[])
# print("H")
# print(H_final[])
# print("LE")
# print(LE_final[])
# print("RN24H")
# print(Rn24h_dB_final[])
# print("LE24H")
# print(LE24h_dB_final[])
# print("EF")
# print(EF[])
# print("ET24H")
# print(ET24h_dB[])

output.evapo<-stack(zom, ustar_final, rah, H_final, LE_final, Rn24h_dB_final, LE24h_dB_final, EF, ET24h_dB)

output.names <- c('zom', 'ustar_after', 'Rah_after', 'H', 'LatentHF', 'Rn24h', 'LatentHF24h', 'EF', 'ET24h')

output.path <- "Testes/Teste11/R_OUTPUT"

names(output.evapo) <- output.names
writeRaster(output.evapo, output.path, overwrite=TRUE, format="GTiff", varname="R_OUTPUT", varunit="daily", longname="R_OUTPUT", xname="lon", yname="lat", bylayer=TRUE, suffix="names")
print("Depois do writeRaster")
	#print(proc.time())