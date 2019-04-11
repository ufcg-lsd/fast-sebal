library(R.utils)
library(raster)
library(rgdal)
library(maptools)
library(ncdf4)
library(sp)
library(snow)
library(here)

hotPixelSelection <- function(Rn, G, TS, NDVI){

	HO<-Rn[]-G[] # Read as a Vector
	x<-TS[][(NDVI[]>0.15 &!is.na(NDVI[]))  & (NDVI[]<0.20 &!is.na(NDVI[])) ] # Returns a vector
	x<-x[x>273.16]
	TS.c.hot<-sort(x)[round(0.95*length(x))] # Returns one value
	HO.c.hot<-HO[][(NDVI[]>0.15 &!is.na(NDVI[])) & (NDVI[]<0.20 &!is.na(NDVI[])) & TS[]==TS.c.hot] # Returns one value
	print("NDVI ente 0.15 e 0.20")
	print(length(x))
	print("Temperatura escolhida")
	print(TS.c.hot)
	print("HO dos que tem essa temp")
	print(HO.c.hot)
	if (length(HO.c.hot)==1){
		ll.hot<-which(TS[]==TS.c.hot & HO[]==HO.c.hot)
		xy.hot <- xyFromCell(TS, ll.hot)
		ll.hot.f<-cbind(as.vector(xy.hot[1,1]), as.vector(xy.hot[1,2]))
	  }else{
		HO.c.hot.min<-sort(HO.c.hot)[ceiling(0.25*length(HO.c.hot))]
		HO.c.hot.max<-sort(HO.c.hot)[ceiling(0.75*length(HO.c.hot))]
		ll.hot<-which(TS[]==TS.c.hot & HO[]>HO.c.hot.min & HO[]<HO.c.hot.max)
                xy.hot <- xyFromCell(TS, ll.hot)
		print("Localizacao dos pixels")
		print(xFromCol(TS, xy.hot[,1]))
		print(yFromRow(TS, xy.hot[,2]))
		print(xy.hot)
		NDVI.hot<-extract(NDVI,xy.hot, buffer=105)
		NDVI.hot.2<-NDVI.hot[!sapply(NDVI.hot, is.null)]
		print("Extract result")
		print(NDVI.hot.2)
		NDVI.hot.cv <- sapply(NDVI.hot.2,sd, na.rm=TRUE)/sapply(NDVI.hot.2, mean, na.rm=TRUE)
                print("CVs")
		print(NDVI.hot.cv)
		i.NDVI.hot.cv<-which.min(NDVI.hot.cv)
		ll.hot.f<-cbind(as.vector(xy.hot[i.NDVI.hot.cv,1]), as.vector(xy.hot[i.NDVI.hot.cv,2]))
	}

	return(ll.hot.f)
}

coldPixelSelection <- function(Rn, G, TS, NDVI){
		
	HO<-Rn[]-G[] # Read as a Vector
	
	x<-TS[][(NDVI[]<0 &!is.na(NDVI[])) & !is.na(HO)]
	x<-x[x>273.16]
	
	TS.c.cold<-sort(x)[round(0.5*length(x))]
	
	HO.c.cold<-HO[(NDVI[]<0 & !is.na(NDVI[])) & TS[]==TS.c.cold & !is.na(HO)]
	print("NDVI menor que 0")
	print(length(x))
	print("Temperatura escolhida")
	print(TS.c.cold)
	print("HOs dos que tem essa temp")
	print(HO.c.cold)

	if (length(HO.c.cold)==1){
		ll.cold<-which(TS[]==TS.c.cold & HO==HO.c.cold)
		xy.cold <- xyFromCell(TS, ll.cold)
		ll.cold.f<-cbind(as.vector(xy.cold[1,1]), as.vector(xy.cold[1,2]))
	}else{
		HO.c.cold.min<-sort(HO.c.cold)[ceiling(0.25*length(HO.c.cold))]
		HO.c.cold.max<-sort(HO.c.cold)[ceiling(0.75*length(HO.c.cold))]
		
		ll.cold<-which(TS[]==TS.c.cold & (HO>HO.c.cold.min &!is.na(HO)) & (HO<HO.c.cold.max & !is.na(HO)))
		print("Localizacao")
		xy.cold <- xyFromCell(TS, ll.cold)
		print(xFromCol(TS, xy.cold[,1]))
		print(yFromRow(TS, xy.cold[,2]))
		print(xy.cold)
		NDVI.cold<-extract(NDVI,xy.cold, buffer=105)
		print("Extract result")
		print(NDVI.cold)
		NDVI.cold.2<-NDVI.cold[!sapply(NDVI.cold, is.null)]
		
		# Maximum number of neighboring pixels with $NVDI < 0$
		t<-function(x){ sum(x<0,na.rm = TRUE)}
		n.neg.NDVI<-sapply(NDVI.cold.2,t)
		print(n.neg.NDVI)
		i.NDVI.cold<-which.max(n.neg.NDVI)
		
		ll.cold.f<-cbind(as.vector(xy.cold[i.NDVI.cold,1]), as.vector(xy.cold[i.NDVI.cold,2]))
	}

	return(ll.cold.f)
}

windVelocity200 <- function(){
	#VER A DECLARAÇÃO DAS CONSTANTES HC, K, X
	zom.est <- hc*0.12

	# Friction velocity at the station (ustar.est)
	ustar.est <- k*table.sw$V6[2]/log((x)/zom.est)
	
	# Velocity 200 meters
	u200 <- ustar.est/k*log(200/zom.est)

	return(u200)
}

####################### Selection of reference pixels ###################################

# Getting the rasters output
Rn<-raster("Testes/Teste04/Rn.tif")
TS<-raster("Testes/Teste04/TS.tif")
NDVI<-raster("Testes/Teste04/NDVI.tif")
G<-raster("Testes/Teste04/G.tif")
alb<-raster("Testes/Teste04/alb.tif")
fic.sw <- "input/station.csv"
table.sw <- (read.csv(fic.sw, sep=";", header=FALSE, stringsAsFactors=FALSE))

####################################CONSTANTS################################################

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

ll.hot.f <- hotPixelSelection(Rn, G, TS, NDVI)
ll.cold.f <- coldPixelSelection(Rn, G, TS, NDVI)

# Location of reference pixels (hot and cold)
ll_ref<-rbind(ll.hot.f[1,],ll.cold.f[1,])
colnames(ll_ref)<-c("long", "lat")
rownames(ll_ref)<-c("hot","cold")
#print(proc.time())
print(ll_ref)
####################################################################################

# Velocity 200 meters
u200 <- windVelocity200()

# Zom for all pixels
#zom<-exp(azom+bzom*NDVI[]) # Changed from Raster to Vector
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

base_ref<-stack(NDVI,TS,Rn,G,ustar,rah) # Raster
nbase<-c("NDVI","TS","Rn","G")
names(base_ref)<-c(nbase,"ustar","rah")

value.pixels.ref<-extract(base_ref,ll_ref)
rownames(value.pixels.ref)<-c("hot","cold")
H.hot<-value.pixels.ref["hot","Rn"]-value.pixels.ref["hot","G"]  
value.pixel.rah<-value.pixels.ref["hot","rah"]
print(value.pixels.ref)

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
    b<-dt.hot/(value.pixels.ref["hot","TS"]-value.pixels.ref["cold","TS"]) # Value
    a<- -b*(value.pixels.ref["cold","TS"]-273.15) # Value
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
    rah.hot<-extract(rah,matrix(ll_ref["hot",],1,2)) # Value
    value.pixel.rah<-c(value.pixel.rah,rah.hot) # Value
    print("after rah")
    print(i)
    print(value.pixel.rah)
    Erro<-(abs(1-rah.hot.0/rah.hot)>=0.05)
    i<-i+1
}

L_final <- NDVI
L_final[] <- L

y_0.1_final <- NDVI
y_0.1_final[] <- y_0.1

y_2_final <- NDVI
y_2_final[] <- y_2

x200_final <- NDVI
x200_final[] <- x200

psi_0.1_final <- NDVI
psi_0.1_final[] <- psi_0.1

psi_2_final <- NDVI
psi_2_final[] <- psi_2

psi_200_final <- NDVI
psi_200_final[] <- psi_200


#print(proc.time())
print("After rah correct")
# End sensible heat flux (H)

# Hot and cold pixels
dt.hot<-H.hot*rah.hot/(rho*cp)                  
b<-dt.hot/(value.pixels.ref["hot","TS"]-value.pixels.ref["cold","TS"]) 
a<- -b*(value.pixels.ref["cold","TS"]-273.15)                          

#print(proc.time())

# All pixels
H <-rho*cp*(a+b*(TS[]-273.15))/rah[] # Vector 
H[(H>(Rn[]-G[]) &!is.na(H))]<-(Rn[]-G[])[(H>(Rn[]-G[]) &!is.na(H))] # Vector

ustar_final <- NDVI
ustar_final[] <- ustar
print("Passou ustar final")

H_final <- NDVI
H_final[] <- H
print("Passou H final")

#print(proc.time())

# Instant latent heat flux (LE)
LE<-Rn[]-G[]-H
LE_final <- NDVI
LE_final[] <- Rn[]-G[]-H
print("Passou LE")
# Upscalling temporal
dr<-(1/d_sun_earth$dist[Dia.juliano])^2 		# Inverse square of the distance on Earth-SOL
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

output.evapo<-stack(zom, ustar_final, rah, H_final, LE_final, Rn24h_dB_final, LE24h_dB_final, EF, ET24h_dB)

output.names <- c('zom', 'ustar_after', 'Rah_after', 'H', 'LatentHF', 'Rn24h', 'LatentHF24h', 'EF', 'ET24h')

output.path <- "Testes/Teste11/R_OUTPUT"

names(output.evapo) <- output.names
writeRaster(output.evapo, output.path, overwrite=TRUE, format="GTiff", varname="R_OUTPUT", varunit="daily", longname="R_OUTPUT", xname="lon", yname="lat", bylayer=TRUE, suffix="names")
print("Depois do writeRaster")