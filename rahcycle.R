library(R.utils)
library(raster)
library(rgdal)
library(maptools)
library(ncdf4)
library(sp)
library(snow)
library(here)

Rn<-raster("Testes/Teste06/meuRn.tif")
TS<-raster("Testes/Teste06/meuTS.tif")
NDVI<-raster("Testes/Teste06/meuNDVI.tif")
G<-raster("Testes/Teste06/meuG.tif")
alb<-raster("Testes/Teste06/meuAlbedo.tif")
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

hot.pixel.rn <- Rn[5, 7]
hot.pixel.g <- G[5, 7]
hot.pixel.ts <- TS[5, 7]
cold.pixel.ts <- TS[3, 9]
print(c(hot.pixel.rn, hot.pixel.g, hot.pixel.ts, cold.pixel.ts))

H.hot<-hot.pixel.rn-hot.pixel.g 
value.pixel.rah<-rah[5, 7]

i<-1
Erro<-TRUE
print("Before rah correct")
	#print(proc.time())
	# Beginning of the cycle stability
while(Erro){
  rah.hot.0<-value.pixel.rah[i] # Value
  print("rah.hot.0")
  # Hot and cold pixels      
  dt.hot<-H.hot*rah.hot.0/(rho*cp) # Value
  b<-dt.hot/(hot.pixel.ts-cold.pixel.ts) # Value
  a<- -b*(cold.pixel.ts-273.15) # Value
  print("before H")
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
  rah.hot <- rah[5, 7]
  value.pixel.rah<-c(value.pixel.rah,rah.hot) # Value
  print("after rah")
  print(i)
  print(value.pixel.rah)
  Erro<-(abs(1-rah.hot.0/rah.hot)>=0.05)
  print(Erro)
  i<-i+1
}
	
	#print(proc.time())
print("After rah correct")
	# End sensible heat flux (H)
	
dt.hot<-H.hot*rah.hot.0/(rho*cp) # Value
b<-dt.hot/(hot.pixel.ts-cold.pixel.ts) # Value
a<- -b*(cold.pixel.ts-273.15) # Value                      
	
	#print(proc.time())
	
	# All pixels
H <-rho*cp*(a+b*(TS[]-273.15))/rah[] # Vector 
H[(H>(Rn[]-G[]) &!is.na(H))]<-(Rn[]-G[])[(H>(Rn[]-G[]) &!is.na(H))] # Vector
	
ustar_final <- NDVI
ustar_final[] <- ustar
print("Passou ustar final")

	#rah_final <- NDVI
	#rah_final[] <- rah
	#print("Passou rah final")

H_final <- NDVI
	#H_final[] <-rho*cp*(a+b*(TS[]-273.15))/rah[] # Vector 
	#H_final[(H>(Rn[]-G[]) &!is.na(H))]<-(Rn[]-G[])[(H>(Rn[]-G[]) &!is.na(H))] # Vector
H_final[] <- H
print("Passou H final")

	#print(proc.time())

	# Instant latent heat flux (LE)
LE<-Rn[]-G[]-H
LE_final <- NDVI
LE_final[] <- Rn[]-G[]-H
print("Passou LE")

MTL <- read.table("input/MTL.txt", skip=0, nrows=140, sep="=", quote = "''", as.is=TRUE) # Reading MTL File
fic <- substr(MTL$V2[MTL$V1 == grep(pattern="LANDSAT_SCENE_ID", MTL$V1, value=T)], 3, 23)
Dia.juliano <- as.numeric(substr(fic, 14, 16))	#Julian Day

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
print("Passou Rn24h")
	# Evapotranspiration fraction Bastiaanssen
EF<-NDVI
EF[]<-LE/(Rn[]-G[])
	
	# Sensible heat flux 24 hours (H24h), nao foi usado
	# H24h_dB<-(1-EF[])*Rn24h_dB
	
	# Latent Heat Flux 24 hours (LE24h)
LE24h_dB<-EF[]*Rn24h_dB
LE24h_dB_final <- NDVI
LE24h_dB_final[] <- EF[]*Rn24h_dB
print("Passou LE24h")
	
	# Evapotranspiration 24 hours (ET24h)
ET24h_dB<-NDVI
ET24h_dB[]<-LE24h_dB*86400/((2.501-0.00236* (max(table.sw$V7[])+min(table.sw$V7[]))/2)*10^6)
print("Antes de salvar")
	#print(proc.time())
	
print("ZOM")
print(zom[])
print("USTAR")
print(ustar_final[])
print("RAH")
print(rah[])
print("H")
print(H_final[])
print("LE")
print(LE_final[])
print("RN24H")
print(Rn24h_dB_final[])
print("LE24H")
print(LE24h_dB_final[])
print("EF")
print(EF[])
print("ET24H")
print(ET24h_dB[])

output.evapo<-stack(zom, ustar_final, rah, H_final, LE_final, Rn24h_dB_final, LE24h_dB_final, EF, ET24h_dB)

output.names <- c('zom', 'ustar_after', 'Rah_after', 'H', 'LatentHF', 'Rn24h', 'LatentHF24h', 'EF', 'ET24h')

output.path <- "Testes/Teste06/R_OUTPUT"

names(output.evapo) <- output.names
writeRaster(output.evapo, output.path, overwrite=TRUE, format="CDF", varname=fic, varunit="daily", longname=fic, xname="lon", yname="lat", bylayer=TRUE, suffix="names")
print("Depois do writeRaster")
	#print(proc.time())
fic <- "R_OUTPUT"	
	# Opening old EF NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_EF.nc",sep="")
nc<-nc_open(var_output, write=TRUE,readunlim=FALSE,verbose=TRUE,auto_GMT=FALSE,suppress_dimvals=FALSE)
	
	# New EF file name
file_output<-paste("Testes/Teste06","/",fic,"_EF.nc",sep="")
oldEFValues<-ncvar_get(nc,fic)
newEFValues<-ncvar_def("EF","daily",list(dimLonDef,dimLatDef,tdim),longname="EF",missval=NaN,prec="double")
nc_close(nc)
newEFNCDF4<-nc_create(file_output,newEFValues)
ncvar_put(newEFNCDF4,"EF",oldEFValues,start=c(1,1,1),count=c(raster.elevation@ncols,raster.elevation@nrows,1))
nc_close(newEFNCDF4)
	
	#print(proc.time())
	
	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_ET24h.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New ET24h file name
file_output<-paste("Testes/Teste06","/",fic,"_ET24h.nc",sep="")
oldET24hValues<-ncvar_get(nc,fic)
newET24hValues<-ncvar_def("ET24h","daily", list(dimLonDef, dimLatDef, tdim), longname="ET24h", missval=NaN, prec="double")
nc_close(nc)
newET24hNCDF4<-nc_create(file_output,newET24hValues)
ncvar_put(newET24hNCDF4, "ET24h", oldET24hValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newET24hNCDF4)
	
	#print(proc.time())
	#FIXME:
	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_zom.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New zom file name
file_output<-paste("Testes/Teste06","/",fic,"_zom.nc",sep="")
oldzomValues<-ncvar_get(nc,fic)
newzomValues<-ncvar_def("zom","daily", list(dimLonDef, dimLatDef, tdim), longname="zom", missval=NaN, prec="double")
nc_close(nc)
newzomNCDF4<-nc_create(file_output,newzomValues)
ncvar_put(newzomNCDF4, "zom", oldzomValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newzomNCDF4)


	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_ustar_after.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New ustar_after file name
file_output<-paste("Testes/Teste06","/",fic,"_ustar_after.nc",sep="")
oldustar_afterValues<-ncvar_get(nc,fic)
newustar_afterValues<-ncvar_def("ustar_after","daily", list(dimLonDef, dimLatDef, tdim), longname="ustar_after", missval=NaN, prec="double")
nc_close(nc)
newustar_afterNCDF4<-nc_create(file_output,newustar_afterValues)
ncvar_put(newustar_afterNCDF4, "ustar_after", oldustar_afterValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newustar_afterNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_Rah_after.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New Rah_after file name
file_output<-paste("Testes/Teste06","/",fic,"_Rah_after.nc",sep="")
oldRah_afterValues<-ncvar_get(nc,fic)
newRah_afterValues<-ncvar_def("Rah_after","daily", list(dimLonDef, dimLatDef, tdim), longname="Rah_after", missval=NaN, prec="double")
nc_close(nc)
newRah_afterNCDF4<-nc_create(file_output,newRah_afterValues)
ncvar_put(newRah_afterNCDF4, "Rah_after", oldRah_afterValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newRah_afterNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_H.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New H file name
file_output<-paste("Testes/Teste06","/",fic,"_H.nc",sep="")
oldHValues<-ncvar_get(nc,fic)
newHValues<-ncvar_def("H","daily", list(dimLonDef, dimLatDef, tdim), longname="H", missval=NaN, prec="double")
nc_close(nc)
newHNCDF4<-nc_create(file_output,newHValues)
ncvar_put(newHNCDF4, "H", oldHValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newHNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_LatentHF.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New LatentHF file name
file_output<-paste("Testes/Teste06","/",fic,"_LatentHF.nc",sep="")
oldLatentHFValues<-ncvar_get(nc,fic)
newLatentHFValues<-ncvar_def("LatentHF","daily", list(dimLonDef, dimLatDef, tdim), longname="LatentHF", missval=NaN, prec="double")
nc_close(nc)
newLatentHFNCDF4<-nc_create(file_output,newLatentHFValues)
ncvar_put(newLatentHFNCDF4, "LatentHF", oldLatentHFValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newLatentHFNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_Rn24h.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New Rn24h file name
file_output<-paste("Testes/Teste06","/",fic,"_Rn24h.nc",sep="")
oldRn24hValues<-ncvar_get(nc,fic)
newRn24hValues<-ncvar_def("Rn24h","daily", list(dimLonDef, dimLatDef, tdim), longname="Rn24h", missval=NaN, prec="double")
nc_close(nc)
newRn24hNCDF4<-nc_create(file_output,newRn24hValues)
ncvar_put(newRn24hNCDF4, "Rn24h", oldRn24hValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newRn24hNCDF4)
	
	#print(proc.time())

	# Opening old ET24h NetCDF
var_output<-paste("Testes/Teste06","/",fic,"_LatentHF24h.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)
	
	# New LatentHF24h file name
file_output<-paste("Testes/Teste06","/",fic,"_LatentHF24h.nc",sep="")
oldLatentHF24hValues<-ncvar_get(nc,fic)
newLatentHF24hValues<-ncvar_def("LatentHF24h","daily", list(dimLonDef, dimLatDef, tdim), longname="LatentHF24h", missval=NaN, prec="double")
nc_close(nc)
newLatentHF24hNCDF4<-nc_create(file_output,newLatentHF24hValues)
ncvar_put(newLatentHF24hNCDF4, "LatentHF24h", oldLatentHF24hValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newLatentHF24hNCDF4)