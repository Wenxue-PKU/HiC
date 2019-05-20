#########################################################################
#  Funtion "Mat": main algorithm using spline					#
#########################################################################
#	For given 4 points(P0,P1,P2,P3), 						#
#	calculate the interpolation curve: 						#
#	the interpolated curve lies only between P1 and P2;			#
## Input: 											#
# 	P0, P1, P2, P3 = 4 given points						#
#	T = tension, which controls the smoothness level (0~1)		#
#	    (0: linear interpolation, 1 = smoothest curve)			#
#	N = number of intervals (spline is evaluated at N+1 points)		#
## Output:											#
#	Mat = matrix of (N+1) coordinates 						#
#		which lies between P1 and P2, on the interpolated curve	#
#		where Mat[1,] = P1 and Mat[N+1,] = P2				#
#########################################################################

Mat = function(P0, P1, P2, P3, T, N){
	du = 1/N;
	temp = matrix(0, length(P0), N+1)
	for (k in 1:(N+1)){
		u = (k-1)*du;
		temp[,k] = evalcrdnd(P0, P1, P2, P3, T, u)
	}
	return(t(temp))
}
##########################################################################
#	Sub-funtion to implement "Mat"						 #
##########################################################################
evalcrdnd = function(P0, P1, P2, P3, T, u){
	s = (1-T)/2;
	temp1 = c(-s, 2-s, s-2, s);
	temp2 = c(2*s, s-3, 3-2*s, -s);
	temp3 = c(-s, 0, s, 0);
	temp4 = c(0,1,0,0);
	MC = rbind(temp1, temp2, temp3, temp4);
	
	G = matrix(0, 4, length(P0))
	for (i in  1:length(P0)){
		G[,i] = c(P0[i], P1[i], P2[i], P3[i]);
	}
	
	U = cbind(u^3, u^2, u, 1);
	PU = U%*%MC%*%G;
	return(PU);
}

##########################################################################
##   Final functions									 #
##########################################################################
# Input: Px, Py (2 coornidation vectors (x and y) with the same length	 #
#	   T = tension, the same as T in Mat					 #
#	   N = # of intervals								 #
#												 #
# Output: plot and coordinates of interpolated curve				 #
##########################################################################	

## 2D spline
Spl2.int = function(Px,Py, T, N){
	if (length(Px)!=length(Py)) {
    print("ERROR: Dimension doesn't match!")
	}else{ 
		n = length(Px)
		Px1 = c(Px[1], Px, Px[length(Px)]);
		Py1 = c(Py[1], Py, Py[length(Py)]);
		Out1 = NULL;
		for (k in 1:(n-1)){
			Out = Mat(c(Px1[k],Py1[k]),c(Px1[k+1],Py1[k+1]),c(Px1[k+2],Py1[k+2]),c(Px1[k+3],Py1[k+3]),T,N);
    	Out1 = rbind(Out1, Out[1:N,])
		}
		Out1 = rbind(Out1, c(Px[n], Py[n]))
	}
	#plot(Px, Py, xlab = "x", ylab= "y")
	#lines(Out1[,1], Out1[,2])
	return(Out1)
}

# Px = c(35, 16, 15, 25, 40, 65, 50, 60, 80);
# Py = c(47, 40, 15, 36, 15, 25, 40, 42, 27);
# temp = Spl2.int(Px, Py, 0, 10)


# ## 3D spline
# library(scatterplot3d)
# Pz = c(47, 20, 55, 46, 35, 25, 40, 32, 27);
# 
# Spl3.int = function(Px,Py,Pz,T, N){
# 	n = length(Px)
# 	Px1 = c(Px[1], Px, Px[length(Px)]);
# 	Py1 = c(Py[1], Py, Py[length(Py)]);
# 	Pz1 = c(Pz[1], Pz, Pz[length(Pz)]);
# 	#scatterplot3d(Px, Py, Pz)
# 		Out1 = NULL;
# 		for (k in 1:(n-1)){
# 			Out = Mat(c(Px1[k],Py1[k], Pz1[k]),c(Px1[k+1],Py1[k+1], Pz1[k+1]),
# 					c(Px1[k+2],Py1[k+2], Pz1[k+2]), c(Px1[k+3],Py1[k+3], Pz1[k+3]),T,N);
# 			Out1 = rbind(Out1, Out[1:N,])
# 		}
# 	Out1 = rbind(Out1, c(Px[n], Py[n], Pz[n]))
# 	temp = scatterplot3d(Out1[,1], Out1[,2], Out1[,3], type= "l", 
# 		xlab = "x", ylab = "y", zlab = "z", highlight.3d = T)
# 	temp$points3d(Px, Py, Pz)
# 	return(Out1)
# }
# temp = Spl3.int(Px, Py, Pz, 0, 20)
