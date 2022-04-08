#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rinternals.h>
#include<R_ext/Rdynload.h>

extern void SmoothEWS( 
	double *pData, 		//TxN matrix (vectorized) containing data to smooth
	int *pT, 		//Length of periodogram to smooth
	int *pN, 		//Number of periodograms to smooth
	int *pM, 		//Bandwidth parameter (# of weights either side of the centre)
	double *pWts, 		//Smoothing kernel weights
	int *pContribute,	//Vector length N: 1 -> calculate score, otherwise -> don't 
	double *pEps, 		//Tolerance to apply to smoothed periodogram
	double *pScore, 	//RETURN: GCV-score
	int *pErr 		//RETURN: error code
		);

void PsiJl(
	double *ppsi, 	//Wavelet function vector
	int *pL, 	//Length of ppsi
	int *pJ, 	//Number of levels
	double *pPsiJl	//RETURN: Wavelet autocor. fn between the fine level and all
		);

void Psijl(
	double *ppsi, 	//Wavelet function vector
	int *pL, 	//Length of ppsi
	int *pJ, 	//Number of levels
	int *pj, 	//Level index
	int *pl, 	//Level index
	double *pPsijl	//RETURN: Wavelet autocorrelation between levels j & l
		);

void PsiJL(
	double *pPsiJl,	//Wavelet autocorrelation function between the fine level and all
	int *pL, 	//length of wavelet function vector
	int *pJ, 	//Number of levels
	double *pPsiJL,	//RETURN:  Wavelet auto correlation funcition
	int *pErr	//RETURN: Error code
		);

extern  void AutoCorr(
	double *ppsi, 	//Wavelet function vector
	int *pL, 	//Length of ppsi
	int *pJ, 	//Number of levels
	double *pPsiJL, //RETURN: Wavelet AutoCor fn at all levels
	int *pErr	//RETURN: Error Code
		);

void A_lam_jlh(
	double *pPsi_jl,	//Wavelet autocorrelation between levels j & l
	int *pminjl,		//Location of lowest non-0 in pPsi_jl
	int *pmaxjl,		//Location of highest non-0 in pPsi_jl
	double *pPsi_h,		//Wavelet autocorrelation for level h & itself
	int *pminh, 		//Location of lowest non-0 in pPsi_h
	int *pmaxh,		//Location of highest non-0 in pPsi_h
	int *pL,		//2L+1 is the length of the vector pPsi_jl
	double *pA		//RETURN: Wavelet autocorr inner product, A^{\lambda}_{j,l,h}
		);
		
extern  void WaveCorrInnerProd(
	double *pPsi,		//Wavelet autocorrelation values (vectorized (2L+1)xJxJ array)
	int *pmin,		//Length J index vector place of lowest non-0 in pPsi per level pair
	int *pmax,		//Length J index vector place of higher non-0 in pPsi per level pair
	int *pL,		//2L+1 is the length of the vector pPsi_jl
	int *pJ,		//Number of levels
	double *pWACIP,		//RETURN: Wave. AutoCor. Inner Prod. (vectorized (2L+1)xJxJxJ array)
	int *pErr		//RETURN: error code
		);

extern  void SmoothCovEst(
	double *pSpq,		//Wavelet cross spectrum between series p&q (vectorized TxJ matix)
	double *pSpp,		//Wavelet spectrum of series p (vectorized TxJ matix)
	double *pSqq,		//Wavelet spectrum of series q (vectorized TxJ matix)
	double *pWACIP,		//Wavelet AutoCorr. Inner Product (vectorizd (2T+1)xJxJxJ array)
	int *pJ,		//Number of levels
	int *pT,		//Time series length
	int *pM,		//Smooth kernel bandwidth parameter
	double *pWts,		//Smoothing kernel weights
	double *pSmoothCov,	//RETURN: Smoothed covariance matrices (vectorized TxJxJ array)
	int *pErr		//RETURN: Error code
		);
	
void E_dpjk_dqlm(
	double *pSpec,	//Wavelet (Cross-) Spectrum
	double *pWACIP,	//Wavelet Autocorrelation Inner Product
	int *pj, 	//Level index
	int *pk, 	//Location index
	int *pl, 	//Level index
	int *pm, 	//Location index
	int *pJ, 	//Number of levels
	int *pT, 	//Time series length
	double *pSumAS	//RETURN: \sum_h A^{k-m}_{jlh} S^{(p,q)}_h ( k+m / 2T )
		);	

void CovRaw_Ijk_Ilm(
	double *pSpq,	//Wavelet Cross-Spectrum
	double *pSpp, 	//Wavelet Spectrum
	double *pSqq, 	//Wavelet Spectrum
	double *pWACIP,	//Wavelet Autocorrelation Inner Product
	int *pj, 	//Level index
	int *pk, 	//Location index
	int *pl, 	//Level index
	int *pm, 	//Location index
	int *pT, 	//Time series length
	int *pJ, 	//Number of levels
	double *pCov	//RETURN: cov( I^{(p,q)}_{j,k} , I^{(p,q)}_{l,m} )
		);

void VMV(
	double **ppMAT, 	//Square DxD matrix
	double *pVec,		//Vector of length D 
	int *pD, 		//Dimnesion of the matrix
	int *pTLcorner, 	//Offset indicating the first element of the matrix
	double *pVMV		//RETRUN: vector-matrix-vector product
		);

void SmoothCov_l1_l2(
	double *pWts,		        //Smoothing kernel weights
	int *pD,			//Width of kernel (length of pWts)
	int *pM,			//Kernel width parameter (width = 2M+1)
	int *pl1,			//Level index
	int *pl2,			//Level index
	int *pT,			//Time series length
	int *pJ,			//Number of levels
	double *pSpq,			//Wavelet cross-sectrum
	double *pSpp,			//Wavelet spectrum of p^th series
	double *pSqq,			//Wavelet spectrum of q^th series
	double *pWACIP,			//Wavelet Autocorrelation Inner Product
	double *pSmoothCov_l1_l2	//RETURN: Smoothed Covariance at levels l1 & l2
		);

void FillCovMatrix(
	double **ppCovMat,      //UPDATE: DxD square matrix
	int *pD,		//Dimension of ppCovMat, width of smoothing kernel
	int *pM,		//Smoothing kernel span, D=2M+1
	int *pTLcorner,		//ppCovMat offset to top left corner of matrix
	int *pFillAll,		//Indicator: 1 -> fill all of the matrix, otherwise -> fill edges
	int *pk,		//Location index
	int *pl1,		//Level index
	int *pl2,		//Level index
	int *pT,		//Time series length
	int *pJ,		//Number of levels
	double *pSpq, 		//Wavelet cross-spectrum
	double *pSpp, 		//Wavelet spectrum of p^th series
	double *pSqq, 		//Wavelet spectrum of q^th series
	double *pWACIP		//Wavelet Autocorrelation Inner Product
		);

//#############################################################################

//Smooth the Evolutionary wavelet spectrum and the Generalized Cross-Validation 
//  gamma deviance criterion(CALLED FROM R)
extern  void SmoothEWS(
	double *pData,		//TxN matrix (vectorized) containing data to smooth
	int *pT,		//Length of periodogram to smooth
	int *pN,		//Number of periodograms to smooth
	int *pM,		//Bandwidth parameter (number of weights either side of the centre)
	double *pWts,		//Smoothing kernel weights
	int *pContribute,	//Vector length N: 1 -> calculate score, otherwise -> don't 
	double *pEps,		//Tolerance to apply to smoothed periodogram
	double *pScore,		//RETURN: GCV-score
	int *pErr		//RETURN: error code
		){

	//ERROR CODES:
	//	0	:	No error
	//	101	:	Smoothing kernel too wind for the time series
	//	102	:	Smoothing kernel weights do not sum to 1 (within 1e-6 tolerance)
	//	201	:	There exists a 0 or -ve data value when calculating the GCV score
	//	202	:	There exists a 0 or -ve smoothed value when calculating the GCV score
	//	203	:	Central smooth weight is 1, resulting in a NaN GCV score
	
	int n, i, k, K, start, width;
	double q, Ratio, WtTotal, GCV, divisor, Smooth_tmp, Data_tmp;
	double *pSmooth;

	width = 2 * *pM + 1;
	*pScore = 0.0;
	
	if(width>=*pT){
		//Kernel too long for time series
		*pErr = 101;
		return;
	}
	
	WtTotal = -1.0;
	for(i=0;i<width;i++) WtTotal += pWts[i];
	if(WtTotal<0) WtTotal = -WtTotal;
	if(WtTotal > 0.000001){
		//weights do not sum to 1
		*pErr = 102;
		return;
	}
	
	if((1-pWts[*pM]) <= 0){
		*pErr = 203;
		return;
	}

	divisor = *pT * (1 - pWts[*pM]) * (1 - pWts[*pM]);
	
	for(n=0;n<*pN;n++){
		start = n * *pT;
		
		pSmooth = (double*)calloc(*pT,sizeof(double));
		GCV = 0.0;
		
		for(k=0;k<*pT;k++){
			pSmooth[k] = 0.0;
			for(i=0;i<width;i++){
				K = k+i-*pM;
				//Reflect boundaries
				if(K<0) K = -K;
				if(K>=*pT) K = 2*(*pT-1)-K;
				pSmooth[k] += pData[start+K] * pWts[i];
			}
			
			if(pContribute[n] == 1){
				//Calculate the GCV score
				Smooth_tmp = pSmooth[k];
				Data_tmp = pData[start+k];
				if(Smooth_tmp < *pEps) Smooth_tmp = *pEps;
				if(Data_tmp < *pEps) Data_tmp = *pEps;
				q = 1.0;
				if(k==0 || k==(*pT-1)) q = 0.5;
				
				if(Data_tmp <= 0){
					//Else GCV score would otherwise return a NaN
					*pErr = 201;
					free(pSmooth);
					return;
				}
				
				if(Smooth_tmp<=0){
					//GCV score would otherwise return a NaN
					*pErr = 202;
					free(pSmooth);
					return;
				}
				Ratio = Data_tmp/Smooth_tmp;
				GCV += q*(Ratio - log(Ratio)-1);
			}
		}
	
		if(pContribute[n] == 1)	*pScore += GCV/divisor;
		
		for(k=0;k<*pT;k++) pData[start+k] = pSmooth[k];
		free(pSmooth);
	}
	*pErr = 0;
	return;
}

//#############################################################################

// Wavelet Auto-correlation Function Between Level J & all cases for l
void PsiJl(
	double *ppsi, 	//Wavelet function vector
	int *pL, 	//Length of ppsi
	int *pJ, 	//Number of levels
	double *pPsiJl	//RETURN: Wavelet autocor. fn between the fine level and all
		){

	int l, j, place;
	j = *pJ - 1;
	for(l=0; l<*pJ; l++){
		place = (*pJ - l - 1) * (2 * *pL + 1);
		Psijl(ppsi, pL, pJ, &j, &l, &(pPsiJl[place])); // j>=l
	}
	return;
}

// Wavelet Auto-correlation Function Between Levels j & l
void Psijl(
	double *ppsi, 	//Wavelet function vector
	int *pL, 	//Length of ppsi
	int *pJ, 	//Number of levels
	int *pj, 	//Level index
	int *pl, 	//Level index
	double *pPsijl	//RETURN: Wavelet autocorrelation between levels j & l
		){

	double sq2;
	double val, scale;
	int i, unitl, unitj, C, Cl, Cj, Cjl, tau, k, kt;
  
	scale = 1;
	sq2 = sqrt(2);
	for(i=0;i<=*pj;i++) scale *= sq2;
	for(i=0;i<=*pl;i++) scale *= sq2;
	unitl = 1;
	for(i=0;i<(*pJ-*pl-1);i++) unitl *= 2;
	unitj = 1;
	for(i=0;i<(*pJ-*pj-1);i++) unitj *= 2;
	C = *pL/2;
	Cl = C/unitl;
	Cj = C/unitj;
	Cjl = Cl+Cj;
	for(tau = -Cjl; tau<=Cjl; tau++){
		val = 0;
		for(k = -Cj; k < Cj; k++){
			kt = k+tau;
			if(kt >= -Cl && kt<Cl){
				val += ppsi[(k+Cj)*unitj]*ppsi[(kt+Cl)*unitl];
			}
		}
		pPsijl[tau+*pL] = val/scale;
	}
	return;
}

// Wavelet Auto-correlation Function Between All Level Pairs 
//  - evaluation based on thining the coarse level wavelet
void PsiJL(
	double *pPsiJl,	//Wavelet autocorrelation function between the fine level and all
	int *pL, 	//length of wavelet function vector
	int *pJ, 	//Number of levels
	double *pPsiJL,	//RETURN: Wavelet auto correlation funcition
	int *pErr	//RETURN: Error code
		){
	//ERROR CODES
	//	0	:	No Error
	//	666	:	Attempting to extract a value that is outsize vector valid range

	int i, j, l, d, start_get, Cj, unitj, start_put_jl, start_put_lj, tau;
	double val;
	
	for(j=0;j<*pJ;j++){
		for(l=0;l<=j;l++){
			d = j-l;
			unitj = 1;
			for(i=0;i<(*pJ-j-1);i++) unitj *= 2;
			Cj = (int)(*pL / unitj);
			start_get = d * (2 * *pL + 1);
			start_put_jl = (l * *pJ + j) * (2 * *pL + 1); // j>=l
			start_put_lj = (j * *pJ + l) * (2 * *pL + 1); //j<=l
			for(tau=-Cj; tau<=Cj; tau++){
				if(tau*unitj+*pL < 0 || tau*unitj+*pL >= 2 * *pL + 1){
					*pErr = 666;
					return;
				}
				val = pPsiJl[start_get + tau*unitj + *pL];
				pPsiJL[start_put_jl + tau + *pL] = val;
				if(j!=l){
					pPsiJL[start_put_lj - tau + *pL] = val;
				}
			}
		}
	}
    *pErr = 0;
	return;
}

// Wavelet Auto-correlation Function Between All Level Pairs 
//  - envelope command (CALLED FROM R)
extern  void AutoCorr(
	double *ppsi, 	//Wavelet function vector
	int *pL, 	//Length of ppsi
	int *pJ, 	//Number of levels
	double *pPsiJL, //RETURN: Wavelet Autocor fn for all level pairs
	int *pErr	//RETURN: Error Code
		){
	//ERROR CODES
	//	0	:	No Error
	//	666	:	Attempting to extract a value that is outsize vector valid range
	
	double *pPsiJl;
	pPsiJl = (double*)calloc(*pJ * (2 * *pL + 1),sizeof(double));
	PsiJl(ppsi, pL, pJ, pPsiJl);
	
	PsiJL(pPsiJl, pL, pJ, pPsiJL, pErr);
	free(pPsiJl);
	
	if(*pErr!=0){
		return;
	}
	
	*pErr = 0;
	return;
}

// Wavelet Auto-correlation Inner Product - at levels j, l & h for all lambda cases
void A_lam_jlh(
	double *pPsi_jl,	//Wavelet autocorrelation between levels j & l
	int *pminjl,		//Location of lowest non-0 in pPsi_jl
	int *pmaxjl,		//Location of highest non-0 in pPsi_jl
	double *pPsi_h,		//Wavelet autocorrelation for level h & itself
	int *pminh, 		//Location of lowest non-0 in pPsi_h
	int *pmaxh,		//Location of highest non-0 in pPsi_h
	int *pL,		//2L+1 is the length of the vector pPsi_jl
	double *pA		//RETURN: Wavelet autocorr inner product, A^{\lambda}_{j,l,h}
		){
  
	int lambda, tau, lambda_tau, tau_min, tau_max, cont;
	
	for(lambda=-*pL; lambda<=*pL; lambda++){
		pA[lambda + *pL] = 0.0;
		
		//Constraint on tau to ensure -L <= lambda-tau <= L
        tau_min = lambda-*pL;
        if(tau_min < -*pL) tau_min = -*pL;
        tau_max = lambda + *pL;
        if(tau_max) tau_max = *pL;
	
		//Additional constrains, work with range where values are non-0.
        if(*pminh-*pL > tau_min) tau_min = *pminh - *pL;
        if(tau_min < lambda - (*pmaxjl - *pL)) tau_min = lambda - (*pmaxjl - *pL);

		cont = 1;
		for(tau=tau_min; cont==1; tau++){
			lambda_tau = lambda-tau;
			if(tau <= *pmaxh-*pL && lambda_tau >= *pminjl-*pL && tau<=tau_max){
				pA[lambda + *pL] += pPsi_jl[lambda_tau + *pL]*pPsi_h[tau + *pL];
			}else{
				cont=0;
			}
		}
	}
  
	return;
}

// Wavelet Auto-correlation Inner Product - All cases (CALLED FROM R)
extern  void WaveCorrInnerProd(
	double *pPsi,	//Wavelet autocorrelation values (vectorized (2L+1)xJxJ array)
	int *pmin,	//Length J index vector place of lowest non-0 in pPsi per level pair
	int *pmax,	//Length J index vector place of higher non-0 in pPsi per level pair
	int *pL,	//2L+1 is the length of the vector pPsi_jl
	int *pJ,	//Number of levels
	double *pWACIP,	//RETURN: Wave. AutoCor. Inner Prod. (vectorized (2L+1)xJxJxJ array)
	int *pErr	//RETURN: error code
		){

	//ERROR CODES:
	//	0	:	No error
	//	301	:	There exists a pmin value that is greater than the pmax pair

	int minh, maxh, minjl, maxjl, h,l,j, start;
	double *pA, *pPsi_h, *pPsi_jl;

	for(j=0;j<(*pJ * *pJ);j++){
		if(pmin[j]>pmax[j]){
			*pErr = 301;
			return;
		}
	}
  
	for(h=0; h<*pJ; h++){
		minh = pmin[h * *pJ + h];
		maxh = pmax[h * *pJ + h];
		pPsi_h = &(pPsi[(h * *pJ + h) * (2 * *pL +1)]);
		for(l=0; l<*pJ; l++){
			for(j=0; j<*pJ; j++){
				minjl = pmin[l * *pJ + j];
				maxjl = pmax[l * *pJ + j];
				pPsi_jl = &(pPsi[(l * *pJ + j) * (2 * *pL +1)]);
				start = ((h * *pJ + l) * *pJ + j) * (2 * *pL +1);
				pA = &(pWACIP[start]);
				A_lam_jlh(pPsi_jl,&minjl,&maxjl,pPsi_h,&minh,&maxh,pL,pA);
			}
		}
	}
	
	*pErr = 0;
	return;
}

//###########################################################################

//Sequence of smoothed periodogram covariance matries at all locations (CALLED FROM R)
extern void SmoothCovEst(
	double *pSpq,		//Wavelet cross spectrum between series p&q (vectorized TxJ matix)
	double *pSpp,		//Wavelet spectrum of series p (vectorized TxJ matix)
	double *pSqq,		//Wavelet spectrum of series q (vectorized TxJ matix)
	double *pWACIP,		//Wavelet AutoCorr. Inner Product (vectorizd (2T+1)xJxJxJ array)
	int *pJ,		//Number of levels
	int *pT,		//Time series length
	int *pM,		//Smooth kernel bandwidth parameter
	double *pWts,		//Smoothing kernel weights
	double *pSmoothCov,	//RETURN: Smoothed covariance matrices (vectorized TxJxJ array)
	int *pErr		//RETURN: Error code
		){

	//ERROR CODES:
	//	0	:	No error
	//	11	:	Incompatible dimensions, T!=2^J
	//	101	:	Smoothing kernel too wide for the time series
	//	102 	:	Smoothing kernel weights do not sum to 1 (within 1e-6 tolerance)
	
	int width, l1, l2, k, start12, start21;
	double total;
	
	l1=1;
	for(k=0; k<*pJ; k++) l1 *= 2;
	if(l1!=*pT){
	  *pErr = 11;
	  return;
	}
	
	width = 2 * *pM + 1;
	if(width>=*pT){
		*pErr = 101;
		return;
	}
	
	total = -1.0;
	for(k=0; k<width; k++) total += pWts[k];
	if(total > 0.000001 || total < -0.000001){
		*pErr = 102;
		return;
	}

	for(l2=0; l2<*pJ; l2++){
		for(l1=0; l1<=l2; l1++){
			start21 = (l2 * *pJ + l1) * *pT;
			SmoothCov_l1_l2(pWts, &width, pM, &l1, &l2, pT, pJ, pSpq, pSpp, pSqq, 
				pWACIP, &pSmoothCov[start21]);
			if(l1!=l2){
				start12 = (l1 * *pJ + l2) * *pT;
				for(k=0;k<*pT;k++){
					pSmoothCov[start12+k] = pSmoothCov[start21+k];
				}
			}
		}
	}
	
	*pErr = 0;
	return;
}

//Expectation of the product between wavelet coefficiets: 
//  first - channel p, level j & location k, second - channel q, level l & location m
void E_dpjk_dqlm(
	double *pSpec,	//Wavelet (Cross-) Spectrum between series p & q
	double *pWACIP,	//Wavelet Autocorrelation Inner Product
	int *pj, 	//Level index
	int *pk, 	//Location index
	int *pl, 	//Level index
	int *pm, 	//Location index
	int *pJ, 	//Number of levels
	int *pT, 	//Time series length
	double *pSumAS	//RETURN: \sum_h A^{k-m}_{jlh} S^{(p,q)}_h ( k+m / 2T )
		){

	int h, Av, loc, WACIP_place, Spq_place, lam;
	double Spq;
	
	loc = *pk + *pm;
	lam = *pm - *pk;
	if(loc % 2 == 0){
		Av=0; 
		Spq_place = 0.5*loc * *pJ;
	}else{
		Av=1;
		Spq_place = 0.5*(loc-1) * *pJ;
	}
	*pSumAS = 0.0;
	
	for(h=0; h<*pJ; h++){
		WACIP_place = ((h * *pJ + *pl) * *pJ + *pj) * (2 * *pT + 1) + lam + *pT;
		if(Av==1){
			Spq = 0.5*(pSpec[Spq_place+h]+pSpec[Spq_place+h+*pJ]);
		}else{
			Spq = pSpec[Spq_place + h];
		}
		*pSumAS += pWACIP[WACIP_place]*Spq;
	}
	return;
}

//Covariance between two raw periodograms:
//  first at level j & location k, second at level l and location m
void CovRaw_Ijk_Ilm(
	double *pSpq,	//Wavelet Cross-Spectrum
	double *pSpp, 	//Wavelet Spectrum
	double *pSqq, 	//Wavelet Spectrum
	double *pWACIP,	//Wavelet Autocorrelation Inner Product
	int *pj, 	//Level index
	int *pk, 	//Location index
	int *pl, 	//Level index
	int *pm, 	//Location index
	int *pT, 	//Time series length
	int *pJ, 	//Number of levels
	double *pCov	//RETURN: cov( I^{(p,q)}_{j,k} , I^{(p,q)}_{l,m} )
		){
	
	double ASpq, ASpp, ASqq;
	E_dpjk_dqlm(pSpq, pWACIP, pj, pk, pl, pm, pJ, pT, &ASpq);
	E_dpjk_dqlm(pSpp, pWACIP, pj, pk, pl, pm, pJ, pT, &ASpp);
	E_dpjk_dqlm(pSqq, pWACIP, pj, pk, pl, pm, pJ, pT, &ASqq);
	*pCov = ASpq*ASpq + ASpp*ASqq;
	return;
}

//Vector-Matrix-Vector multiplication
void VMV(
	double **ppMAT, 	//Square DxD matrix
	double *pVec,		//Vector of length D 
	int *pD, 		//Dimnesion of the matrix
	int *pTLcorner, 	//Offset indicating the first element of the matrix
	double *pVMV		//RETRUN: vector-matrix-vector product
		){
  
	double tmp;
	double *pMcol;
	int d1,d2,l1,l2;
	
	*pVMV = 0.0;
	for(d1=0; d1<*pD; d1++){
		tmp = 0.0;
		l1 = (d1+*pTLcorner) % *pD;
		pMcol = ppMAT[l1];
		for(d2=0; d2<*pD; d2++){
			l2 = (d2+*pTLcorner) % *pD;
			tmp += pMcol[l2] * pVec[d2];
		}
		*pVMV += tmp * pVec[d1];
	}
	
	return;
}

//Covariance between the smoothed periodogram at levels l1 & l2, both at location k
void SmoothCov_l1_l2(
    double *pWts,			//Smoothing kernel weights
	int *pD,			//Width of kernel (length of pWts)
	int *pM,			//Kernel width parameter (width = 2M+1)
	int *pl1,			//Level index
	int *pl2,			//Level index
	int *pT,			//Time series length
	int *pJ,			//Number of levels
	double *pSpq,			//Wavelet cross-sectrum
	double *pSpp,			//Wavelet spectrum of p^th series
	double *pSqq,			//Wavelet spectrum of q^th series
	double *pWACIP,			//Wavelet Autocorrelation Inner Product
	double *pSmoothCov_l1_l2	//RETURN: Smoothed Covariance at levels l1 & l2
		){
	
	int TopLeft, FillAll, k;
	
	//work space matrix ( 2M+1 x 2M+1 )
	double **ppCOV;
	ppCOV = (double**)calloc(*pD,sizeof(double*));
	for(k=0;k<*pD;k++) ppCOV[k] = (double*)calloc(*pD,sizeof(double));
	
	TopLeft = 0; 
	FillAll = 1;
	k=0;
	FillCovMatrix(ppCOV, pD, pM, &TopLeft, &FillAll, &k, pl1, pl2, 
		pT, pJ, pSpq, pSpp, pSqq, pWACIP);
	VMV(ppCOV, pWts, pD, &TopLeft, &(pSmoothCov_l1_l2[k]));
	FillAll = 0;
	
	for(k=1;k<*pT;k++){
		FillCovMatrix(ppCOV, pD, pM, &TopLeft, &FillAll, &k, pl1, pl2, 
			pT, pJ, pSpq, pSpp, pSqq, pWACIP);
		TopLeft++;
		if(TopLeft==*pD) TopLeft = 0;
		VMV(ppCOV, pWts, pD, &TopLeft, &(pSmoothCov_l1_l2[k]));
	}
	
	//free ppCOV
	for(k=0;k<*pD;k++) free(ppCOV[k]);
	free(ppCOV);
    return;
}

//Fills in the raw periodogram matrix ready for smoothing
void FillCovMatrix(
	double **ppCovMat,	//UPDATE: Square DxD matrix
	int *pD,		//Dimension of ppCovMat, width of smoothing kernel
	int *pM,		//Smoothing kernel span, D=2M+1
	int *pTLcorner,		//ppCovMat offset to top left corner of matrix
	int *pFillAll,		//Indicator: 1 -> fill all of the matrix, otherwise -> fill edges
	int *pk,		//Location index
	int *pl1,		//Level index
	int *pl2,		//Level index
	int *pT,		//Time series length
	int *pJ,		//Number of levels
	double *pSpq, 		//Wavelet cross-spectrum
	double *pSpp, 		//Wavelet spectrum of p^th series
	double *pSqq, 		//Wavelet spectrum of q^th series
	double *pWACIP		//Wavelet Autocorrelation Inner Product
		){
	
	int m1, m2, K1, K2, place1, place2;
	double *pCovVec;
	
	if(*pFillAll==1){
		for(m2 = *pk-*pM; m2 <= *pk+*pM; m2++){
		    place2 = *pTLcorner + m2 + *pM - *pk;
			if(place2 >= *pD) place2 -= *pD;
			pCovVec = ppCovMat[place2];
			K2 = m2;
			if(K2<0) K2 = -K2;
			if(K2>=*pT) K2 = 2*(*pT-1)-K2;
			for(m1 = *pk-*pM; m1<= *pk+*pM; m1++){
				K1 = m1;
				if(K1<0) K1 = -K1;
				if(K1>=*pT) K1 = 2*(*pT-1)-K1;
				place1 = *pTLcorner + m1 + *pM - *pk;
				if(place1 >= *pD) place1 -= *pD;
				CovRaw_Ijk_Ilm(pSpq, pSpp, pSqq, pWACIP, pl1, &K1, pl2, &K2, 
					pT, pJ, &(pCovVec[place1]));
			}
		}
	}else{
		pCovVec = ppCovMat[*pTLcorner];
		K2 = *pM+*pk;
		if(K2>=*pT) K2 = 2*(*pT-1)-K2;
		for(m1 = *pk-*pM; m1 <= *pk+*pM; m1++){
			place1 = (*pTLcorner + m1 + *pM - *pk +1) % *pD;
			K1 = m1;
			if(K1<0) K1 = -K1;
			if(K1>=*pT) K1 = 2*(*pT-1)-K1;
			CovRaw_Ijk_Ilm(pSpq, pSpp, pSqq, pWACIP, pl1, &K1, pl2, &K2, 
				pT, pJ, &(pCovVec[place1]));
		}

		K1 = *pM + *pk;
		if(K1>=*pT) K1 = 2 * (*pT-1) - K1;
		for(m2=*pk-*pM; m2<*pk+*pM; m2++){
			place2 = (*pTLcorner + m2 + *pM - *pk + 1) % *pD;
			pCovVec = ppCovMat[place2];
			K2 = m2;
			if(K2<0) K2 = -K2;
			if(K2>=*pT) K2 = 2*(*pT-1)-K2;
			CovRaw_Ijk_Ilm(pSpq, pSpp, pSqq, pWACIP, pl1, &K1, pl2, &K2, 
				pT, pJ, &(pCovVec[*pTLcorner]));
		}
	}
	return;
}

/*static R_NativePrimitiveArgType SmoothEWS_t[] = {
  REALSXP,INTSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP
};

static R_NativePrimitiveArgType AutoCorr_t[] = {
  REALSXP,INTSXP,INTSXP,REALSXP,INTSXP
};
static R_NativePrimitiveArgType WaveCorrInnerProd_t[] = {
  REALSXP,INTSXP,INTSXP,INTSXP,INTSXP,REALSXP,INTSXP
};
static R_NativePrimitiveArgType SmoothCovEst_t[] = {
  REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,INTSXP
};

static R_CMethodDef cMethods[] = {
  {"SmoothEWS", (DL_FUNC) &SmoothEWS, 9, SmoothEWS_t},
  {"AutoCorr", (DL_FUNC) &AutoCorr, 5, AutoCorr_t},
  {"WaveCorrInnerProd", (DL_FUNC) &WaveCorrInnerProd, 7, WaveCorrInnerProd_t},
  {"SmoothCovEst", (DL_FUNC) &SmoothCovEst, 10, SmoothCovEst_t},
  {NULL,NULL,0}
};*/

void R_init_mvLSW_Cfunctions(DllInfo *info){
     
    R_CMethodDef cMethods[] = {
        {"SmoothEWS", (DL_FUNC) &SmoothEWS, 9},
        {"AutoCorr", (DL_FUNC) &AutoCorr, 5},
        {"WaveCorrInnerProd", (DL_FUNC) &WaveCorrInnerProd, 7},
        {"SmoothCovEst", (DL_FUNC) &SmoothCovEst, 10},
        {NULL,NULL,0}
    };
    
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
}
