# Constants
mass 	<- 5.125	# mass of ball in oz
circ 	<- 9.125	# circumference of ball in inches
tau 	<- 10000
beta 	<- 0.0001217	# constant for computing corrected barometric pressure
Cd 	<- 0.33
lambda 	<- 1.512
grav	<- 32.174
dt.c	<- 0.001
t.max	<- 1.1		# Don't simulate beyond this many seconds of elapsed time
y.min	<- -3		# Don't simulate beyond this many feet behind home plate
y.flag	<- 17/12	# front of home plate
x.l	<- (-8.5-1.45)/12
x.r	<- (8.5+1.45)/12

# Air density in lb/ft^3 at temperature temp (in degrees C), elevation elev (in m), 
# pressure p (in mm Hg), relative humidity rh (%)
rho.f <- function(temp,elev,p,rh,svp) {
	1.2929*(273/(temp+273)*(p*exp(-1*beta*elev)-0.3783*rh*svp/100)/760)*0.06261
}

# Constant for ... something. In practice a constant multiple of rho, which is
# calculated as in the above function.
c0.f <- function(rho) {
	# Original calculation: 0.07182*rho*(5.125/mass)*(circ/9.125)^2
	# Since 5.125/mass = 1 and circ/9.125=1, simplifies to:
	0.07182*rho
}

# Saturation vapor pressure at temperature temp (in degrees C)
svp.f <- function(temp) {
	4.5841*exp((18.687-temp/234.5)*temp/(257.14+temp))
}

# Reynolds number for v=100mph at temperature temp (in degrees C), rho (in lb/ft^3)
re100.f <- function(temp,rho) {
	r 	<- rho/0.06261
	re100 	<- r*44.7*(circ/(pi*39.37))*(temp+273.16+120)/(0.000001512*(temp+273.16)^1.5)
	return(re100)
}

#
# pitch.trajectory:
#
# Function for describing the full trajectory of a pitch based on initial conditions.
#
# Takes the "9 parameters" from the "9 parameter fit model" (x0,y0,z0,v0x,v0y,v0z,wx,wy,wz)
# plus a constant c0 that reflects game conditions.
#
# Returns a data frame describing the pitch conditions at each time increment (dt.c).
#

pitch.trajectory <- function(x0,y0,z0,v0x,v0y,v0z,wx,wy,wz,c0) {
	
	# initialize some local variables
	x	<- x0
	y	<- y0
	z	<- z0
	vx	<- v0x
	vy	<- v0y
	vz	<- v0z
	v	<- sqrt(vx^2+vy^2+vz^2)
    omega   <- sqrt(wx^2+wy^2+wz^2)
    romega  <- (circ/2/pi)*omega/12
	S 	<- (romega/v0)
	Cl 	<- 1/(2.32+0.4/S)	# Cross model
#	Cl	<- 0.09+0.6*S		# Hubbard model
	fact.m	<- c0*v*Cl/omega
	aMagx	<- fact.m*(wy*v0z-wz*v0y)
	aMagy	<- fact.m*(wz*v0x-wx*v0z)
	aMagz	<- fact.m*(wx*v0y-wy*v0x)
	fact.d 	<- -1*c0*Cd
	adragx	<- fact.d*v*v0x
	adragy	<- fact.d*v*v0y
	adragz	<- fact.d*v*v0z
	w	<- omega*30/pi
	ax	<- adragx+aMagx
	ay	<- adragy+aMagy
	az	<- adragz+aMagz-grav
	xx	<- x0+vx*dt.c+ax*dt.c*dt.c/2
	yy	<- y0+vy*dt.c+ay*dt.c*dt.c/2
	zz	<- z0+vz*dt.c+az*dt.c*dt.c/2
	t	<- 0

	# initialize data frame of results with initial conditions
	traj.df <- data.frame(t=t,x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,v=v,S=S,Cl=Cl,adragx=adragx,adragy=adragy,
			adragz=adragz,w=w,aMagx=aMagx,aMagy=aMagy,aMagz=aMagz,ax=ax,ay=ay,az=az,
			xx=xx,yy=yy,zz=zz)
	
	# Loop to calculate state at each t+dt
	while((t<=t.max) & (y>=y.min)) {
		t 	<- t+dt.c
		x	<- xx
		y	<- yy
		z	<- zz
		vx	<- vx+ax*dt.c
		vy	<- vy+ay*dt.c
		vz	<- vz+az*dt.c
		v	<- sqrt(vx^2+vy^2+vz^2)
		S	<- (romega/v)*exp(-1*t/tau)
		Cl	<- 1/(2.32+0.4/S)	# Cross model
#		Cl	<- 0.09+0.6*S		# Hubbard model
		adragx	<- fact.d*v*vx
		adragy	<- fact.d*v*vy
		adragz	<- fact.d*v*vz
		w	<- omega*exp(-1*t/tau)*30/pi
		fact.m	<- c0*v*Cl/omega
		aMagx	<- fact.m*(wy*vz-wz*vy)
		aMagy	<- fact.m*(wz*vx-wx*vz)
		aMagz	<- fact.m*(wx*vy-wy*vx)
		ax	<- adragx+aMagx
		ay	<- adragy+aMagy
		az	<- adragz+aMagz-grav
		xx	<- xx+vx*dt.c+ax*dt.c*dt.c/2
		yy	<- yy+vy*dt.c+ay*dt.c*dt.c/2
		zz	<- zz+vz*dt.c+az*dt.c*dt.c/2
		
		traj.n	<- c(t,x,y,z,vx,vy,vz,v,S,Cl,adragx,adragy,adragz,w,aMagx,aMagy,aMagz,ax,ay,az,xx,yy,zz)
		traj.df <- rbind(traj.df,traj.n)
	}

	return(traj.df)
}

# given two sets of (x,y,z) coordinates, interpolate between them based on a target y value
pos.interp.y <- function(v1,v2,y.targ) {
	ratio	<- (v2[2]-y.targ)/(v2[2]-v1[2])
	v	<- v2-ratio*(v2-v1)
	return(v)
}

# given two sets of (x,y,z) coordinates, interpolate between them based on a target z value
pos.interp.z <- function(v1,v2,z.targ) {
	ratio	<- (v2[3]-z.targ)/(v2[3]-v1[3])
	v	<- v2-ratio*(v2-v1)
	return(v)
}

# given two sets of (x,y,z) coordinates, interpolate between them based on a target x value
pos.interp.x <- function(v1,v2,x.targ) {
	ratio	<- (v2[1]-x.targ)/(v2[1]-v1[1])
	v	<- v2-ratio*(v2-v1)
	return(v)
}

pitch.coords.at.y <- function(df,y.targ) {
	row 	<- df[df$y>=y.targ & df$yy<y.targ,]
	v1	<- c(row$x,row$y,row$z)
	v2	<- c(row$xx,row$yy,row$zz)
	v	<- pos.interp.y(v1,v2,y.targ)
	return(v)
}

# Function to determine if a given pair of x,y coordinates lie within the 2D boundaries of home plate.
# Returns true/false.
within.plate.boundaries <- function(x, y) {
	if(y<=(-1.45)) { return(FALSE) }			# coordinates are past lowest point of home plate
	if(y<=8.5/12) { 
		return(ifelse(abs(x)<=y+2.04,TRUE,FALSE))}	# coordinates are within/without the lower triangle; adjustment of 2.04 is for radius of ball: sqrt(1.45^2+1.45^2) gives the max x shift that still leaves a point of the ball touching the plate area.
	if(y<=17/12) { 
		return(ifelse(abs(x)<=x.r,TRUE,FALSE))}	# coordinates are within/without the upper rectangle
	return(FALSE)					# coordinates are above highest point of home plate
}

# Function to determine if pitch trajectory passes through 3D strike zone. Takes two sets of coordinates
# (x1,z1 and x2,z2) as arguments, and the z coordinates of the top and bottom of the strike zone 
# (zt and zb). x1,z1 and x2,z2 are presumed to be at y=17/12 and y=0 respectively.
in.zone.3d <- function(x1,z1,x2,z2,zt,zb) {
	# Note that if the pitch never goes above zb or below zt then it is out of the zone, period
	if((z1<=zb & z2<=zb) | (z1>zt & z2>zt)) { return(0) }
	
	# Also, if the x dimension is never within -8.5/12 and 8.5/12 then it is out of the zone, period
	if((x1<x.l & x2<x.l) | (x1>x.r & x2>x.r)) { return(0) }

	# If it went through front 2D zone plane, it's in the zone
	if(x1>=x.l & x1<=x.r & z1>=zb & z1<=zt) { return(1) }

	# If it went through middle 2D zone plane, it's in the zone
	vm <- pos.interp.y(c(x1,17/12,z1),c(x2,0,z2),8.5/12)
	if(vm[1]>=x.l & vm[1]<=x.r & vm[3]>=zb & vm[3]<=zt) { return(1) }

	# If it did neither but went through projection of home plate at bottom of zone, it's in the zone
	if(z1>=zb & z2<=zb) {
		vb <- pos.interp.z(c(x1,17/12,z1),c(x2,0,z2),zb)
		if(within.plate.boundaries(vb[1],vb[2])) { return(1) }
	}
	
	# If it did none of these but went through projection of home plate at top of zone, it's in the zone
	if(z1>=zt & z2<=zt) {
		vt <- pos.interp.z(c(x1,17/12,z1),c(x2,0,z2),zt)
		if(within.plate.boundaries(vt[1],vt[2])) { return(1) }
	}

	# Failed all tests; it's not in the zone
	return(0)
}

# Function to determine if pitch trajectory passes through traditional 2D strike zone. v is a vector
# of x,y,z coordinates, and y is presumed to be 17/12.
in.zone.2d <- function(x,z,zt,zb) {
	# simple check of x and z dimensions
	if(x>=x.l & x<=x.r & z>=zb & z<=zt) { return(1) }
	return(0)
}

# Return the (x,y,z) coords of the ball when it crosses front and back of plate, as a list of two vectors)
# This will greatly simplify the info that needs to be stored for analysis - instead of a whole data frame
# of coordinates for each individual pitch (themselves stored in a huge data frame), just tack on six columns 
# to the data frame to define these two points, and assume linear travel between them (good enough to fractions
# of an inch).
# NOTE: to unpack a list of vectors, use double-bracket ("[[]]" notation. E.g. if the list is called L,
# L[[1]] returns v1 as a vector and L[[2]] returns v2 as a vector; L[[2]][3] would give the z component of v2.
pitch.plate.crossings <- function(x0,y0,z0,v0x,v0y,v0z,wx,wy,wz,c0) {
	traj.df	<- pitch.trajectory(x0,y0,z0,v0x,v0y,v0z,wx,wy,wz,c0)
	v1	<- pitch.coords.at.y(traj.df,y.flag)
	v2	<- pitch.coords.at.y(traj.df,0)
#	return(list(v1,v2))
	return(c(v1,v2))	# Note: it is easier to work with results after the fact when in a single array rather than a list of arrays
}

# Variation: calculate the change in coordinates as the ball traverses the plate. So we can
# then fit that change to the ACTUAL px/pz values recorded by PITCHf/x (i.e. rather than calculate
# my own starting and ending coordinates, calculate the change and use THEIR starting coordinates)
pitch.plate.movement <- function(x0,y0,z0,v0x,v0y,v0z,v0,wx,wy,wz,c0) {
	traj.df	<- pitch.trajectory(x0,y0,z0,v0x,v0y,v0z,wx,wy,wz,omega,c0)
	v1	<- pitch.coords.at.y(traj.df,y.flag)
	v2	<- pitch.coords.at.y(traj.df,0)
	return(c(v2[1]-v1[1],v2[2]-v1[2],v2[3]-v1[3]))
}
