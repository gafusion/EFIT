/*
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "norm.h"
#include "macro.h"

extern __constant__ float cu_radix[mesh];
extern __constant__ int cu_num_block[mesh];
extern __constant__ int cu_blo_accu[mesh];
extern __constant__ int cu_num_se[1949];
extern __constant__ float cu_currentInver[mesh]; 
*/
__global__ void norm1(float* leas_flux)
{
	float c, d;
	if(blockIdx.x == 0)
	{
		__shared__ float flux[122];
		__shared__ float fflux[2][122];
		flux[threadIdx.x] = leas_flux[threadIdx.x+3000];        // si
		fflux[0][threadIdx.x] = leas_flux[threadIdx.x+3000+122];    //  z
		fflux[1][threadIdx.x] = float(threadIdx.x);  		// r   
//		printf("%d %d\n", threadIdx.x, bpp[0][threadIdx.x]);
		__syncthreads();

		for(unsigned int s=64;s>0;s>>=1)
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<122) 
			{
				c=flux[threadIdx.x];
				d=flux[k];
				if(c<d)
				{
					flux[threadIdx.x]=d;
					fflux[0][threadIdx.x] = fflux[0][k];
					fflux[1][threadIdx.x] = fflux[1][k];
				}
				__syncthreads();
			}
		}
		leas_flux[4] = fflux[0][0];       // z
		leas_flux[5] = fflux[1][0];       // r
	}
	else if(blockIdx.x == 1)
	{
		__shared__ float bp[122];
		__shared__ float bpp[2][122];
		bp[threadIdx.x] = leas_flux[threadIdx.x];        // br^2+bz^2
		bpp[0][threadIdx.x] = leas_flux[threadIdx.x+122];    //  z
		bpp[1][threadIdx.x] = threadIdx.x;  		// r   
//		printf("%d %d\n", threadIdx.x, bpp[0][threadIdx.x]);	      	
		__syncthreads();
		for(unsigned int s=64;s>0;s>>=1)
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<122) 
			{
				c=bp[threadIdx.x];
				d=bp[k];
				if(c>d)
				{
					bp[threadIdx.x]=d;
					bpp[0][threadIdx.x] = bpp[0][k];
					bpp[1][threadIdx.x] = bpp[1][k];
				}
				__syncthreads();
			}	
		}
		leas_flux[0] = bpp[0][0];       // z
		leas_flux[1] = bpp[1][0];       // r
//		printf("%f %f\n", leas_flux[0],leas_flux[1]);
	}
	else
	{
		__shared__ float bp[122];
		__shared__ float bpp[2][122];
		bp[threadIdx.x] = leas_flux[2000+threadIdx.x];
		bpp[0][threadIdx.x] = leas_flux[2000+threadIdx.x+122];
		bpp[1][threadIdx.x] = threadIdx.x;
		__syncthreads();
		for(unsigned int s=64;s>0;s>>=1)
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<122) 
			{
				c=bp[threadIdx.x];
				d=bp[k];
				if(c>d)
				{
					bp[threadIdx.x]=d;
					bpp[0][threadIdx.x] = bpp[0][k];
					bpp[1][threadIdx.x] = bpp[1][k];
				}
				__syncthreads();
			}
		}
		leas_flux[2] = bpp[0][0];
		leas_flux[3] = bpp[1][0];
//		printf("%f %f\n", leas_flux[2],leas_flux[3]);
	}	
}


__global__ void norm3(float* xpoint,float*cuLimFlux)
{
	__shared__ float slimFlux[LIMITER_NUM];
	float c,d;
	slimFlux[threadIdx.x]=cuLimFlux[threadIdx.x];
	__syncthreads();
	for(unsigned int s=64;s>0;s>>=1) 
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<LIMITER_NUM) 
		{
			c=slimFlux[threadIdx.x];
			d=slimFlux[k];
			if(c<d)
				slimFlux[threadIdx.x]=d;			
		}
		__syncthreads();
	}
	if(threadIdx.x ==0)
	{
		cuLimFlux[0] = slimFlux[0];
//		printf("\n\n%f\n\n",xpoint[5]);
		if(xpoint[2]<xpoint[0])
		{
			xpoint[1] = xpoint[0];
		}
		else
		{
			xpoint[1] = xpoint[2];
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
__global__ void limiterFlux(float* limiter, float* cu_result, float* xpoint,float* cuLimFlux)
{
	int rID,zID;
	float flux1,flux2,flux3,flux4;
	float m1,m2,m3,m4;
	float r,z,zx1,zx2,r0,z0;	
	
	r = limiter[threadIdx.x];
	z = limiter[threadIdx.x + LIMITER_NUM];

	__syncthreads();
	
	zx1 = xpoint[8];    // up xpoint
	zx2 = xpoint[6];		// bottom xpoint

	if (z<zx1 && z>zx2)     // the limiter points between the two x-points.
	{
		rID = int((r-igridr-delR)/delR);
		zID = int((z-igridz-delZ)/delZ);

		r0 = (rID+1)*delR+igridr;
		z0 = (zID+1)*delZ+igridz;

		flux1 = cu_result[rID*mesh1+zID];
		flux2 = cu_result[(rID+1)*mesh1+zID];
		flux3 = cu_result[rID*mesh1+zID+1];
		flux4 = cu_result[(rID+1)*mesh1+zID+1];

		m1 = (r0+delR-r)*(z0+delZ-z);
		m2 = (r-r0)*(z0+delZ-z);
		m3 = (r0+delR-r)*(z-z0);
		m4 = (r-r0)*(z-z0);

		cuLimFlux[threadIdx.x] = (m1*flux1+m2*flux2+m3*flux3+m4*flux4)/darea;
	}
	else                   
	{
		cuLimFlux[threadIdx.x] = xpoint[0];
	}

//	printf("%d %f \n", threadIdx.x, cuLimFlux[threadIdx.x]);
}

__global__ void c_error_comput(float* d_cInters, int* d_cPos, float* cu_result,float* cu_bound_flux,float* c_error)
{
	float bound = cu_bound_flux[0];
	int rID,zID;
	float flux1,flux2,flux3,flux4;
	float m1,m2,m3,m4;
	rID=d_cPos[blockIdx.x*2];
	zID=d_cPos[blockIdx.x*2+1];
	if(rID > 0)
	{
		flux1=cu_result[rID*64+zID];
		flux2=cu_result[(rID+1)*64+zID];
		flux3=cu_result[rID*64+zID+1];
		flux4=cu_result[(rID+1)*64+zID+1];
		m1=d_cInters[blockIdx.x*4];
		m2=d_cInters[blockIdx.x*4+1];
		m3=d_cInters[blockIdx.x*4+2];
		m4=d_cInters[blockIdx.x*4+3];
//if(blockIdx.x == 0)
//printf("rID:%d, zID:%d\n", rID, zID);
		c_error[blockIdx.x]=m1*flux1+m2*flux2+m3*flux3+m4*flux4-bound;
	}
	else
		c_error[blockIdx.x] = 0;
}



//////////////////////////////////////////////////////////////////////////////
/*
Funtion data_li:
It is used to reduce the number of eliments. And it also develops the current vector based 
on the flux vector;
fluxbound is the flux at the reference point of the plasma bound, and redu[0] is the flux 
at the magent axis; 
buffer[*][0] stores the original flux vector;
*/
__global__ void data_eli(float* cu_result, float* effect_flux,float* xpoint,float* cuLimFlux,float* right_side,float* cu_bound_flux,int* cu_num_se,float* current_representation)
{
	int kppcur,kffcur,pbdry,ffbdry,num,offset,a,b,c,d,i,j,e,f;
	float fluxbound,normfactor;
	float x[20],k,kk;

	kppcur = (int)current_representation[1];
	kffcur = (int)current_representation[2];
	pbdry = (int)current_representation[3]; 
	ffbdry = (int)current_representation[4];
	num = (int)current_representation[25];

//	float dz;
//	if(threadIdx.x==0&&blockIdx.x == 0)
//		printf("lim%f  bdry%f  ax%f  norm-factor%f\n",cuLimFlux[0],xpoint[2],xpoint[4],xpoint[4]-xpoint[2]);
	if(cuLimFlux[0]>xpoint[1])
		fluxbound = cuLimFlux[0];
	else
		fluxbound = xpoint[1];
	if(threadIdx.x==0&&blockIdx.x == 0)
	{
//		fluxbound = 0.3375518918;
		cu_bound_flux[0] = fluxbound;
		xpoint[11] = xpoint[4]-fluxbound;
	}

//	if(threadIdx.x==0&&blockIdx.x == 0)              
//		printf("lim%20.18lf  bdry%20.18lf  ax%20.18lf  diff%20.18lf\n",cuLimFlux[0],cu_bound_flux[0],xpoint[4],xpoint[11]);

	normfactor=1/(fluxbound-xpoint[4]);
	__shared__ float buffer[mesh][2];
	__shared__ float normedflux[mesh];
//	__shared__ float dzflux[mesh][2];
    	unsigned int index=blockIdx.x*(mesh+1)+threadIdx.x;
	buffer[threadIdx.x][1]=cu_result[index];
	__syncthreads();

//	if(threadIdx.x>0&&threadIdx.x<(mesh-1))
//		dzflux[threadIdx.x][0] = (buffer[threadIdx.x+1][4] - buffer[threadIdx.x-1][4])*INVZ;
//	else 
//		dzflux[threadIdx.x][0] = 0;

	c = (xpoint[6]-igridz-delZ)/delZ;
	d = (xpoint[8]-igridz-delZ)/delZ;
	e = (1.016 - igridr - delR)/delR;
	f = (2.3540 - igridr - delR)/delR;

	if(threadIdx.x<cu_num_block[blockIdx.x])
	{
		offset=cu_blo_accu[blockIdx.x];

		//get the positions of effective flux data
		a=cu_num_se[offset+threadIdx.x];

		//get the realative positon of effective flux data		
		b=threadIdx.x;

		//gather the values of effective fluxs into a new vector buffer[*][2] 
		buffer[b][0]=buffer[a][1];
//		dzflux[b][1] = dzflux[a][0];
		normedflux[b]=(buffer[b][0]-xpoint[4])*normfactor;
		if(normedflux[b]>=0&&normedflux[b]<=1&&a>(c-1)&&a<(d+1)&&blockIdx.x>e&&blockIdx.x<(f+1))
		{
			k = 1;
			kk = 1;
			for (j=0;j<kppcur;j++)
			{
				kk = kk*normedflux[b];	
			}
			for (i=0;i<kppcur;i++)
			{
				x[i] = cu_radix[blockIdx.x]*(k-pbdry*kk);
				k = k*normedflux[b];
//				buffer[b][0]=cu_radix[blockIdx.x];		
////				buffer[b][1]=cu_radix[blockIdx.x]*normedflux[b];
			}

			k = 1;
			kk = 1;
			for (j=0;j<kffcur;j++)
			{
				kk = kk*normedflux[b];	
			}			
			for (i=kppcur;i<kppcur+kffcur;i++)
			{
				x[i] = 1/cu_radix[blockIdx.x]*(k-ffbdry*kk);
				k = k*normedflux[b];		
			}
/*			
//			dz = dzflux[b][1]*normfactor; 
			//			printf("%f  ",dz);

			k=normedflux[b]*normedflux[b];
//			kk=normedflux[b]*normedflux[b]*normedflux[b];
			buffer[b][0]=cu_radix[blockIdx.x];
			buffer[b][1]=cu_radix[blockIdx.x]*normedflux[b];
			buffer[b][2]=1/cu_radix[blockIdx.x];
			buffer[b][3]=normedflux[b]/cu_radix[blockIdx.x];
			buffer[b][4]=k/cu_radix[blockIdx.x];
//			buffer[b][5] = cu_radix[blockIdx.x]*dz;
*/
		}
		else
		{
			for (i=0;i<num;i++)
			{
				x[i] = 0;
			}
		}
		//begin to store the transposed data into global memory
		for (i=0;i<num;i++)
		{
//			effect_flux[offset+b]=buffer[b][0];
//			effect_flux[ELI_NUM+offset+b]=buffer[b][1];
//			effect_flux[2*ELI_NUM+offset+b]=buffer[b][2];
//			effect_flux[3*ELI_NUM+offset+b]=buffer[b][3];
//			effect_flux[4*ELI_NUM+offset+b]=buffer[b][4];
//			effect_flux[5*ELI_NUM+offset+b]=buffer[b][5];
			effect_flux[i*ELI_NUM+offset+b] = x[i];
		}			
	}
}

__global__ void data_eli_init(float* cu_result, float* effect_flux, int* cu_num_se, float* current_representation)
{	
	int kppcur,kffcur,pbdry,ffbdry,num,offset,a,b,i,j;
	__shared__ float buffer[mesh][2];
	__shared__ float normedflux[mesh];
	float x[20],k,kk;

	kppcur = (int)current_representation[1];
	kffcur = (int)current_representation[2];
	pbdry = (int)current_representation[3]; 
	ffbdry = (int)current_representation[4];
	num = (int)current_representation[25];

    	unsigned int index=blockIdx.x*mesh+threadIdx.x;
	buffer[threadIdx.x][1]=cu_result[index];
	__syncthreads();

//	if(threadIdx.x==0&&blockIdx.x==0)
//	{
//		printf("ok!");
//	}

	if(threadIdx.x<cu_num_block[blockIdx.x])
	{
		offset=cu_blo_accu[blockIdx.x];
		//get the positions of effective flux data
		a=cu_num_se[offset+threadIdx.x];
		b=threadIdx.x;
		//gather the values of effective fluxs into a new vector buffer[*][4] 
		buffer[b][0]=buffer[a][1];
		normedflux[b]=buffer[b][0];
		if(normedflux[b]>0&&normedflux[b]<1)
		{
			k = 1;
			kk = 1;
			for (j=0;j<kppcur;j++)
			{
				kk = kk*normedflux[b];	
			}
			for (i=0;i<kppcur;i++)
			{
				x[i]=cu_radix[blockIdx.x]*(k-pbdry*kk);
				k = k*normedflux[b];
//				buffer[b][0]=cu_radix[blockIdx.x];		
////				buffer[b][1]=cu_radix[blockIdx.x]*normedflux[b];
			}
//			kk = pow(normedflux[b],kffcur);
			k = 1;
			kk = 1;
			for (j=0;j<kffcur;j++)
			{
				kk = kk*normedflux[b];	
			}			
			for (i=kppcur;i<kppcur+kffcur;i++)
			{
				x[i]=1/cu_radix[blockIdx.x]*(k-ffbdry*kk);
				k = k*normedflux[b];		
			}
//			for (i=kppcur+kffcur;i<num;i++)
//			{
//				x[i] = 0;
//			}

//			k=normedflux[b]*normedflux[b];
//			kk=normedflux[b]*normedflux[b]*normedflux[b];
//			buffer[b][0]=cu_radix[blockIdx.x];
//			buffer[b][0]=cu_radix[blockIdx.x];
//			buffer[b][1]=cu_radix[blockIdx.x]*normedflux[b];
//			buffer[b][1]=cu_radix[blockIdx.x]*normedflux[b];
//			buffer[b][2]=1/cu_radix[blockIdx.x];
//			buffer[b][2]=1/cu_radix[blockIdx.x];
//			buffer[b][3]=normedflux[b]/cu_radix[blockIdx.x];
//			buffer[b][4]=k/cu_radix[blockIdx.x];
		}
		else
		{
			for (i=0;i<num;i++)
			{
				x[i] = 0;
			}
//			buffer[b][0]=0;		
//			buffer[b][1]=0;
//			buffer[b][2]=0;
//			buffer[b][3]=0;		
//			buffer[b][4]=0;
		}
		//begin to store the transposed data into global memory
		for (i=0;i<num;i++)
		{	
//			effect_flux[offset+b]=buffer[b][0];
//			effect_flux[ELI_NUM+offset+b]=buffer[b][1];
//			effect_flux[2*ELI_NUM+offset+b]=buffer[b][2];
//			effect_flux[3*ELI_NUM+offset+b]=buffer[b][3];
			effect_flux[i*ELI_NUM+offset+b] = x[i];
//			printf("effect_flux %d = %f \n",i*ELI_NUM+offset+b,x[i]);
		}
	}
}

__global__ void compute_converge_init(float* cu_result,float* cu_result_compare,int* converge)
{
	unsigned int index1 = blockIdx.x * mesh + threadIdx.x;
	unsigned int index2 = blockIdx.x * (mesh+1) + threadIdx.x;
	cu_result_compare[index2] = cu_result[index1];

	if(threadIdx.x == 0 && blockIdx.x == 0)
	{
		converge[0] = 0;	
	}
}

__global__ void compute_converge(float* cu_result,float* cu_result_compare)
{
	__shared__ float buffer[mesh];
	float c,d;
	unsigned int index = blockIdx.x*(mesh+1) + threadIdx.x;

	buffer[threadIdx.x] = cu_result[index] - cu_result_compare[index];
	__syncthreads();

	for(unsigned int s=128;s>0;s>>=1) 
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s && k<255) 
		{
			c=buffer[threadIdx.x];
			d=buffer[k];
			if(c<d)
				buffer[threadIdx.x]=d;			
		}
		__syncthreads();
	}
	
	cu_result_compare[index] = cu_result[index];
	
	cu_result_compare[70000 + blockIdx.x] = buffer[0];
}

__global__ void compute_converge_post(float* cu_result_compare,float* xpoint,float* cu_bound_flux,int* converge,float* leas_flux)
{
	__shared__ float buffer[mesh];
	float c,d;
	float error;

	buffer[threadIdx.x] = cu_result_compare[70000 + threadIdx.x];
	if (threadIdx.x < 127)
	{
		buffer[threadIdx.x + 128] = cu_result_compare[70000 + threadIdx.x +128];
	}
	__syncthreads();

	for(unsigned int s=128;s>0;s>>=1) 
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s && k<255) 
		{
			c=buffer[threadIdx.x];
			d=buffer[k];
			if(c<d)
				buffer[threadIdx.x]=d;			
		}
		__syncthreads();
	}
	
	error = 0.001*(xpoint[4] - cu_bound_flux[0]);

	if(threadIdx.x == 0)
	{
		if(buffer[0] < error)
		{
			converge[0] = 1;
			leas_flux[0] = buffer[0];
//			printf("error %f ",buffer[0]);
//			printf("%f %f\n",error,leas_flux[907]);
		}
		else
		{
			converge[0] = 0;

			printf("error %f \n",buffer[0]/error);
//			printf("%d \n",converge[0]);
//			printf("%f %f\n",error,buffer[0]);
		}
	}
}

__global__ void result_regroup(float* cu_result,float* result,float* cu_bound_resut,float* cu_icflux,float* cu_ecflux)
{
	int index1 = (blockIdx.x-1)*256+threadIdx.x-1;
	int index2 = threadIdx.x*257+blockIdx.x;
	if(blockIdx.x>0&&blockIdx.x<256&&threadIdx.x>0&&threadIdx.x<256)
	{
		result[index2] = cu_result[index1];
	}
	else if (blockIdx.x == 0 && threadIdx.x>0 && threadIdx.x<256)
	{
		result[index2] = cu_bound_resut[threadIdx.x-1]+cu_icflux[blockIdx.x*257+threadIdx.x]+cu_ecflux[blockIdx.x*257+threadIdx.x];
	}
	else if (blockIdx.x == 256 && threadIdx.x>0 && threadIdx.x<256)
	{
		result[index2] = cu_bound_resut[3*mesh+threadIdx.x-1]+cu_icflux[blockIdx.x*257+threadIdx.x]+cu_ecflux[blockIdx.x*257+threadIdx.x];
	}
	else if (threadIdx.x == 0 && blockIdx.x>0 && blockIdx.x<256)
	{
		result[index2] = cu_bound_resut[mesh+(blockIdx.x-1)*2]+cu_icflux[blockIdx.x*257+threadIdx.x]+cu_ecflux[blockIdx.x*257+threadIdx.x];
	}
	else if (threadIdx.x == 256 && blockIdx.x>0 && blockIdx.x<256)
	{
		result[index2] = cu_bound_resut[mesh+(blockIdx.x-1)*2+1]+cu_icflux[blockIdx.x*257+threadIdx.x]+cu_ecflux[blockIdx.x*257+threadIdx.x];
	}
	else
	{
		result[index2] = cu_result[blockIdx.x*254+threadIdx.x*254/256];
	}
}

__global__ void jtrepresnt(float* cu_dgfc, float* cu_pro_diag, float* current_representation)
{
	int kppcur,kffcur,pbdry,ffbdry,n;
//	float pa,fa;

	kppcur = (int)current_representation[1];
	kffcur = (int)current_representation[2];
	pbdry = 1 - (int)current_representation[3]; 
	ffbdry = 1 - (int)current_representation[4];

	n = (int)current_representation[26];

//	pa = (float)current_representation[5];
//	fa = (float)current_representation[6];	

	if (blockIdx.x < kppcur*pbdry)
	{
		cu_dgfc[(blockIdx.x+EF_NFCOIL)*(EF_DIAG+n) + EF_DIAG] = current_representation[5+blockIdx.x];
	}
	else if (blockIdx.x < kppcur*pbdry + kffcur*ffbdry)
	{
		cu_dgfc[(blockIdx.x+EF_NFCOIL)*(EF_DIAG+n) + EF_DIAG + pbdry] = current_representation[5+10+blockIdx.x-kppcur*pbdry];
	}
	else
	{
		cu_pro_diag[blockIdx.x + EF_DIAG - (kppcur*pbdry + kffcur*ffbdry)] = 0;		
	}

//	printf("jt constraint added!! \n");
//	if (blockIdx.x == 0)
//	{
//	for (int i = 0;i<EF_DIAG+n;i++)
//	{
//		printf("%f\n",cu_pro_diag[i]);
//	}
//	}
}

__global__ void data_eli_init_spline(double* cu_result, float* effect_flux, int* cu_num_se, float* current_representation,double* cu_psiaxis_gafile,double *cu_psibdry_gafile)
{
	int num,offset,a,b,i,j,num_knot,fluxaxis,fluxbound,normfactor;
	__shared__ float buffer[mesh][2];
	__shared__ float normedflux[mesh];
	float x[10],y[10];

	num = 8;
	num_knot = current_representation[1];

	for (j=0;j<num_knot;j++)
	{
		y[j] = current_representation[2+j];
	}

 //   	unsigned int index=blockIdx.x*mesh+threadIdx.x;
//	buffer[threadIdx.x][1] = cu_result[index];
	unsigned int index = (blockIdx.x+1)*mesh2+threadIdx.x+1;
	buffer[threadIdx.x][1] = cu_result[index];

	fluxaxis = cu_psiaxis_gafile[0];
	fluxbound = cu_psibdry_gafile[0];
	normfactor=1/(fluxbound-fluxaxis);
	
	__syncthreads();

	if(threadIdx.x<cu_num_block[blockIdx.x])
	{
		offset=cu_blo_accu[blockIdx.x];
		//get the positions of effective flux data
		a=cu_num_se[offset+threadIdx.x];
		b=threadIdx.x;
		//gather the values of effective fluxs into a new vector buffer[*][4] 
		buffer[b][0]=buffer[a][1];
//		normedflux[b]=buffer[b][0];
		normedflux[b]=(buffer[b][0]-fluxaxis)*normfactor;
		if(normedflux[b]>=0 && normedflux[b]<=1)
		{
			x[0] = cu_radix[blockIdx.x];
			x[1] = cu_radix[blockIdx.x]*normedflux[b];
			x[2] = cu_radix[blockIdx.x]*sinh(tension*normedflux[b]);
			x[3] = cu_radix[blockIdx.x]*cosh(tension*normedflux[b]);
			x[4] = 1/cu_radix[blockIdx.x];
			x[5] = 1/cu_radix[blockIdx.x]*normedflux[b];
			x[6] = 1/cu_radix[blockIdx.x]*sinh(tension*normedflux[b]);
			x[7] = 1/cu_radix[blockIdx.x]*cosh(tension*normedflux[b]);
		}
		else
		{
			for (i=0;i<num;i++)
			{
				x[i] = 0;
			}
//			buffer[b][0]=0;		
//			buffer[b][1]=0;
//			buffer[b][2]=0;
//			buffer[b][3]=0;		
//			buffer[b][4]=0;
		}
		//begin to store the transposed data into global memory
		for (i=0;i<num;i++)
		{	
			for (j=0;j<num_knot-1;j++)
			{
				if (normedflux[b] >= y[j] && normedflux[b] < y[j+1])
				{
					effect_flux[j*num*ELI_NUM + i*ELI_NUM + offset+b] = x[i];
				}
			}
		}
	}
}

__global__ void spline_constrain(float* current_representation,float* x,float* y,float* z)
{
	float psi,factor;
	int index1,index0,num_knots;

	psi = x[blockIdx.x+1];
	num_knots = current_representation[1];
	index0 = EF_NFCOIL+(num_knots-1)*8;	
	index1 = blockIdx.x*index0*6 + index0*3*threadIdx.x;
	factor = powf(-1,threadIdx.y);

	y[index1 + index0*0 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 0] = factor*1;
	y[index1 + index0*0 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 1] = factor*psi;
	y[index1 + index0*0 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 2] = factor*sinh(tension*psi);
	y[index1 + index0*0 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 3] = factor*cosh(tension*psi);

	y[index1 + index0*1 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 0] = factor*0;
	y[index1 + index0*1 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 1] = factor*1;
	y[index1 + index0*1 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 2] = factor*tension*cosh(tension*psi);
	y[index1 + index0*1 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 3] = factor*tension*sinh(tension*psi);

	y[index1 + index0*2 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 0] = factor*0;
	y[index1 + index0*2 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 1] = factor*0;
	y[index1 + index0*2 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 2] = factor*tension*tension*sinh(tension*psi);
	y[index1 + index0*2 + EF_NFCOIL + blockIdx.x*8 + threadIdx.x*4 + threadIdx.y*8 + 3] = factor*tension*tension*cosh(tension*psi);
}

__global__ void pprime_flux_function_integral_constrain(float* current_representation,float* pres_psi,float* pres_left,float* y)
{
	float x;
	int i,j,num_knot,index;

	x = pres_psi[threadIdx.x];
	num_knot = current_representation[1];
	index = (EF_NFCOIL+(num_knot-1)*8)*threadIdx.x;

	for (j=0;j<num_knot-1;j++)  // integral value on each knot segment
	{
		if ((x >= y[j]) && (x < y[j+1]))
		{
			pres_left[index + EF_NFCOIL + j*8 + 0] = (y[j+1] - x)/darea;
			pres_left[index + EF_NFCOIL + j*8 + 1] = 0.5*(y[j+1]*y[j+1]-x*x)/darea;
			pres_left[index + EF_NFCOIL + j*8 + 2] = (cosh(tension*y[j+1]) - cosh(tension*x))/(tension*darea);
			pres_left[index + EF_NFCOIL + j*8 + 3] = (sinh(tension*y[j+1]) - sinh(tension*x))/(tension*darea);

			for (i=j+1;i<num_knot-1;i++)
			{
				pres_left[index + EF_NFCOIL + i*8 + 0] = (y[i+1] - y[i])/darea;
				pres_left[index + EF_NFCOIL + i*8 + 1] = 0.5*(y[i+1]*y[i+1]-y[i]*y[i])/darea;
				pres_left[index + EF_NFCOIL + i*8 + 2] = (cosh(tension*y[i+1]) - cosh(tension*y[i]))/(tension*darea);
				pres_left[index + EF_NFCOIL + i*8 + 3] = (sinh(tension*y[i+1]) + sinh(tension*y[i]))/(tension*darea);
			}
		}
	}
}

__global__ void jphi_constrain_1(float* cu_result,float* location,float* tani,float* coti,float* xpoint,int* num,float* nphi,float* xy_surface)
{
	// search flux surface for j_phi constrain
	int ra,rb,za,zb,r_axi,z_axi,i,c,d,e,f;
	float si1,si2,flux,dflux,aflux,factor;

	z_axi = int(location[4]);
	r_axi = int(location[5])+68;

	aflux = xpoint[4];
	dflux = xpoint[11];

	c = (xpoint[6]-igridz-delZ)/delZ;
	d = (xpoint[8]-igridz-delZ)/delZ;

	e = (1.016 - igridr - delR)/delR;
	f = (2.3540 - igridr - delR)/delR;
	
	__syncthreads();

	if (blockIdx.x<(1+numbdry/8))
	{
		if (threadIdx.x <= f+1-r_axi)
		{
			ra = r_axi+threadIdx.x;
			za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
			si1 = cu_result[ra*mesh1+za];

			rb = r_axi+threadIdx.x+1;
			zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();		
	}
	else if (blockIdx.x>=(1+numbdry/8) && blockIdx.x<(2*numbdry/8))
	{
		if (threadIdx.x <= d+1-z_axi)
		{
			za = z_axi+threadIdx.x;
			ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
			si1 = cu_result[ra*mesh1+za];

			zb = z_axi+threadIdx.x+1;
			rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(2*numbdry/8) && blockIdx.x<(3*numbdry/8))
        {
		if (threadIdx.x <= d+1-z_axi)
		{
                	za = z_axi+threadIdx.x;
       	      		ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
			si1 = cu_result[ra*mesh1+za];

           	        zb = z_axi+threadIdx.x+1;
        	        rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
               		si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(3*numbdry/8) && blockIdx.x<(4*numbdry/8))
	{
		if (threadIdx.x <= r_axi-e)
		{
			ra = r_axi-threadIdx.x;
			za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
			si1 = cu_result[ra*mesh1+za];

			rb = r_axi-threadIdx.x-1;
			zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
			si2 =cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(4*numbdry/8) && blockIdx.x<(5*numbdry/8))
        {
		if (threadIdx.x <= r_axi-e)
		{
       	        	ra = r_axi-threadIdx.x;
               		za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
			si1 = cu_result[ra*mesh1+za];

                	rb = r_axi-threadIdx.x-1;
                	zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
                	si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(5*numbdry/8) && blockIdx.x<(6*numbdry/8))
	{
		if (threadIdx.x <= z_axi-c)
		{
			za = z_axi-threadIdx.x;
			ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
			si1 = cu_result[ra*mesh1+za];
			
			zb = z_axi-threadIdx.x-1;
			rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(6*numbdry/8) && blockIdx.x<(7*numbdry/8))
        {
		if (threadIdx.x <= z_axi-c)
		{
			za = z_axi-threadIdx.x;
               		ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
			si1 = cu_result[ra*mesh1+za];

                	zb = z_axi-threadIdx.x-1;
                	rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
                	si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else
	{
		if (threadIdx.x <= f+1-r_axi)
		{
			ra = r_axi+threadIdx.x;
			za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
			si1 = cu_result[ra*mesh1+za];

			rb = r_axi+threadIdx.x+1;
			zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)<0 && (flux-si2)>0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ; 
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
}

__global__ void jphi_constrain_2(float*cu_result,float* nphi,float* xy_surface,float* cu_current_representation,float* y,float* jphi_left)
{
	int j,i,num_knots,index,ra,za;
	float r,z,dl,a,b,phi;
	float factor1,factor2,bp,br,bz;
	float flux1,flux2,flux3,flux4,flux_a,flux_b,flux_c,flux_d;

	__shared__ float left[numbdry][8];
	__shared__ float sum[numbdry];

	num_knots = cu_current_representation[1];
	phi = nphi[blockIdx.x];
	r = xy_surface[blockIdx.x*2*numbdry+threadIdx.x];
	z = xy_surface[(blockIdx.x*2+1)*numbdry+threadIdx.x];

	ra = int((r-igridr)/delR)-1;
	za = int((z-igridz)/delZ)-1;

	factor1 = (r-igridr)/delR - (ra+1);
	factor2 = (z-igridz)/delZ - (za+1);

	flux1 = cu_result[ra*mesh1+za];
	flux2 = cu_result[(ra+1)*mesh1+za];
	flux3 = cu_result[(ra+1)*mesh1+za+1];
	flux4 = cu_result[ra*mesh1+za+1];

	flux_a = factor1*flux2+(1-factor1)*flux1;	
	flux_b = factor2*flux3+(1-factor2)*flux2;
	flux_c = factor1*flux3+(1-factor1)*flux4; 
	flux_d = factor2*flux4+(1-factor2)*flux1;

	br = (flux_c-flux_a)/(-r*delZ);
	bz = (flux_b-flux_d)/(r*delR);
	bp = sqrt(br*br+bz*bz);

	if (threadIdx.x < (numbdry-1))
	{
		a = xy_surface[blockIdx.x*2*numbdry+threadIdx.x+1];
		b = xy_surface[(blockIdx.x*2+1)*numbdry+threadIdx.x+1];
		sum[threadIdx.x] = sqrt((r-a)*(r-a)+(z-b)*(z-b));
	}
	else if (threadIdx.x == (numbdry-1))
	{
		a = xy_surface[blockIdx.x*2*numbdry];
		b = xy_surface[(blockIdx.x*2+1)*numbdry];
		sum[threadIdx.x] = sqrt((r-a)*(r-a)+(z-b)*(z-b)); 
	}
	__syncthreads();

	if (threadIdx.x == 0)
	{
		dl = (sum[threadIdx.x] + sum[numbdry-1])/2;
	}
	else if (threadIdx.x > 0)
	{	
		dl = (sum[threadIdx.x] + sum[threadIdx.x-1])/2;
	}
	__syncthreads();

	for (j=0;j<num_knots-1;j++)
	{
		if ((phi >= y[j]) && (phi < y[j+1]))
		{
			left[threadIdx.x][0] = dl/bp;
			left[threadIdx.x][1] = phi * dl/bp;
			left[threadIdx.x][2] = sinh(tension*phi) * dl/bp;
			left[threadIdx.x][3] = cosh(tension*phi) * dl/bp;
			left[threadIdx.x][4] = dl/(r*r*bp);
			left[threadIdx.x][5] = phi * dl/(r*r*bp);
			left[threadIdx.x][6] = sinh(tension*phi) * dl/(r*r*bp);
			left[threadIdx.x][7] = cosh(tension*phi) * dl/(r*r*bp);
			sum[threadIdx.x] = dl/(r*bp);
			i = j;
		}
	}	
	__syncthreads();

	for(unsigned int s=64;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<numbdry)
                {
                        left[threadIdx.x][0] += left[k][0];
			left[threadIdx.x][1] += left[k][1];
			left[threadIdx.x][2] += left[k][2];
			left[threadIdx.x][3] += left[k][3];
			left[threadIdx.x][4] += left[k][4];
			left[threadIdx.x][5] += left[k][5];
			left[threadIdx.x][6] += left[k][6];
			left[threadIdx.x][7] += left[k][7];
			sum[threadIdx.x] += sum[k];
                }
		__syncthreads();
	}

	if (threadIdx.x < 8)
	{
		index = blockIdx.x * (EF_NFCOIL+(num_knots-1)*8);
//		x[blockIdx.x*(num_knot[0]-1)*8 + i*8 + threadIdx.x] = left[0][threadIdx.x]/sum[0];
		jphi_left[index + EF_NFCOIL + i*8 + threadIdx.x] = left[0][threadIdx.x]/(sum[0]*darea);
	}
}

__global__ void change_Pos_2(float* cu_dgfc, float* d_A, float* current_representation,int* diag_num_all,float* cu_pro_diag,float* right_side)
{
	int index,index2;
       	index = blockIdx.x * diag_num_all[0] + threadIdx.x;
	index2 = blockIdx.x * EF_DIAG + threadIdx.x;
	d_A[index] = cu_dgfc[index2];

	if (blockIdx.x == 0)
	{
		right_side[threadIdx.x] = cu_pro_diag[threadIdx.x];
	}
}

__global__ void change_Pos_mse(float* cu_dgfc, float* d_A, float* current_representation,int* diag_num_all)
{
	int index;
       	index = blockIdx.x * diag_num_all[0] + threadIdx.x;
	d_A[index] = cu_dgfc[index];
}

__global__ void change_Pos_jphi(float* jphi_left, float* d_A, float* j_phi, float* j_phi_uncer,float* diag_right,float* current_representation,int* diag_num_all,int* offset)
{
	int index1,index2;

	float m,w,k;

       	index1 = blockIdx.x * diag_num_all[0] + offset[0] + threadIdx.x;
	index2 = threadIdx.x * (EF_NFCOIL+(int)current_representation[25]) + blockIdx.x;

	m = j_phi[threadIdx.x]*1000;
	w = j_phi_uncer[threadIdx.x];
	k = 1/sqrt(m*w*m*w);

	d_A[index1] = jphi_left[index2]*k;

	if (blockIdx.x == 0)
	{
		diag_right[index1] = m*k;
//		printf("TH %d M %f W %f K %f \n",threadIdx.x+1,m,w,k);
	}
}

__global__ void change_Pos_pres_init(float* pres_left, float* d_A, float* pressure, float* pres_uncer, float* diag_right, float* current_representation,int* diag_num_all,int* offset)
{
	int index1,index2;
	
	float m,w,k;
	
       	index1 = blockIdx.x * diag_num_all[0] + offset[0] + threadIdx.x;
	index2 = threadIdx.x * (EF_NFCOIL+(int)current_representation[25]) + blockIdx.x;

	m = pressure[threadIdx.x];
	w = pres_uncer[threadIdx.x];
	k = 1/sqrt(m*w*m*w);

	d_A[index1] = pres_left[index2]*k;

	if (blockIdx.x == 0)
	{
		diag_right[index1] = m*k;
	}
}

__global__ void change_Pos_pres(float* pres_left, float* d_A, float* pressure, float* pres_uncer, float* diag_right, float* current_representation,int* diag_num_all,int* offset,float* xpoint)
{
	int index1,index2;
	
	float m,w,k;
	
       	index1 = blockIdx.x * diag_num_all[0] + offset[0] + threadIdx.x;
	index2 = threadIdx.x * (EF_NFCOIL+(int)current_representation[25]) + blockIdx.x;

	m = pressure[threadIdx.x];
	w = pres_uncer[threadIdx.x];
	k = 1/sqrt(m*w*m*w);

	d_A[index1] = pres_left[index2]*k*xpoint[11];

	if (blockIdx.x == 0)
	{
		diag_right[index1] = m*k;
	}
}

__global__ void change_Pos_spline(float* spline_left, float* d_A, float* current_representation,int* diag_num_all,int* offset,float* cu_pro_diag)
{
	int index1,index2;
       	index1 = blockIdx.x * diag_num_all[0] + offset[0] + threadIdx.x;
	index2 = threadIdx.x * (EF_NFCOIL+(int)current_representation[25]) + blockIdx.x;

	d_A[index1] = spline_left[index2];
}

__global__ void spline_svd_d_a_init(float* d_A,int* diag_num_all)
{
	int index1;
       	index1 = blockIdx.x * diag_num_all[0] + threadIdx.x;
	d_A[index1] = 0;
}

__global__ void data_eli_spline(float* cu_result, float* effect_flux,float* xpoint,float* cuLimFlux,float* right_side,float* cu_bound_flux,int* cu_num_se,float* current_representation)
{
	int offset,a,b,c,d,i,j,e,f,num_knot,num;
	float fluxbound,normfactor;
	float x[10],y[10];

	num_knot = current_representation[1];
	num = 8;

	for (j=0;j<num_knot;j++)
	{
		y[j] = current_representation[2+j];
	}

//	float dz;
	if(threadIdx.x==0&&blockIdx.x == 0)
		printf("lim%f  xpoint%f  ax%f  norm-factor%f\n",cuLimFlux[0],xpoint[2],xpoint[4],xpoint[4]-xpoint[2]);
//
	if(cuLimFlux[0]>xpoint[1])
		fluxbound = cuLimFlux[0];
	else
		fluxbound = xpoint[1];

	if(threadIdx.x==0&&blockIdx.x == 0)
	{
		cu_bound_flux[0] = fluxbound;
		xpoint[11] = xpoint[4]-fluxbound;
	}

//	if(threadIdx.x==0&&blockIdx.x == 0)              
//		printf("lim%20.18lf  bdry%20.18lf  ax%20.18lf  diff%20.18lf\n",cuLimFlux[0],cu_bound_flux[0],xpoint[4],xpoint[11]);

	normfactor=1/(fluxbound-xpoint[4]);
	__shared__ float buffer[mesh][2];
	__shared__ float normedflux[mesh];
    	unsigned int index = blockIdx.x*(mesh+1)+threadIdx.x;
	buffer[threadIdx.x][1] = cu_result[index];
	__syncthreads();

	c = (xpoint[6]-igridz-delZ)/delZ;
	d = (xpoint[8]-igridz-delZ)/delZ;
	e = (1.016 - igridr - delR)/delR;
	f = (2.3540 - igridr - delR)/delR;

	if(threadIdx.x<cu_num_block[blockIdx.x])
	{
		offset=cu_blo_accu[blockIdx.x];

		//get the positions of effective flux data
		a=cu_num_se[offset+threadIdx.x];

		//get the realative positon of effective flux data		
		b=threadIdx.x;

		//gather the values of effective fluxs into a new vector buffer[*][2] 
		buffer[b][0]=buffer[a][1];
//		dzflux[b][1] = dzflux[a][0];
		normedflux[b]=(buffer[b][0]-xpoint[4])*normfactor;
		if(normedflux[b]>=0&&normedflux[b]<=1&&a>(c-1)&&a<(d+1)&&blockIdx.x>e&&blockIdx.x<(f+1))
		{
			x[0] = cu_radix[blockIdx.x];
			x[1] = cu_radix[blockIdx.x]*normedflux[b];
			x[2] = cu_radix[blockIdx.x]*sinh(tension*normedflux[b]);
			x[3] = cu_radix[blockIdx.x]*cosh(tension*normedflux[b]);
			x[4] = 1/cu_radix[blockIdx.x];
			x[5] = 1/cu_radix[blockIdx.x]*normedflux[b];
			x[6] = 1/cu_radix[blockIdx.x]*sinh(tension*normedflux[b]);
			x[7] = 1/cu_radix[blockIdx.x]*cosh(tension*normedflux[b]);
		}
		else
		{
			for (i=0;i<num;i++)
			{
				x[i] = 0;
			}
		}
		//begin to store the transposed data into global memory
		
		for (i=0;i<num;i++)
		{
			for (j=0;j<num_knot-1;j++)
			{
				if (normedflux[b] >= y[j] && normedflux[b] < y[j+1])
				{
					effect_flux[j*num*ELI_NUM+i*ELI_NUM+offset+b] = x[i];
				}
			}
		}			
	}
}

__global__ void jphi_constrain_init_1(float* cu_result,float* location,float* tani,float* coti,float* xpoint,int* num,float* nphi,float* xy_surface)
{
	// search flux surface for j_phi constrain
	int ra,rb,za,zb,r_axi,z_axi,i;
	float si1,si2,flux,factor;

	z_axi = 127;
	r_axi = 121;
	
	__syncthreads();

	if (blockIdx.x<(1+numbdry/8))
	{
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
			ra = r_axi+threadIdx.x;
			za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
			si1 = cu_result[ra*mesh+za];

			rb = r_axi+threadIdx.x+1;
			zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
			si2 = cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();		
	}
	else if (blockIdx.x>=(1+numbdry/8) && blockIdx.x<(2*numbdry/8))
	{
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
			za = z_axi+threadIdx.x;
			ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
			si1 = cu_result[ra*mesh+za];

			zb = z_axi+threadIdx.x+1;
			rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
			si2 = cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(2*numbdry/8) && blockIdx.x<(3*numbdry/8))
        {
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
                	za = z_axi+threadIdx.x;
       	      		ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
			si1 = cu_result[ra*mesh+za];

           	        zb = z_axi+threadIdx.x+1;
        	        rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
               		si2 = cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(3*numbdry/8) && blockIdx.x<(4*numbdry/8))
	{
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
			ra = r_axi-threadIdx.x;
			za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
			si1 = cu_result[ra*mesh+za];

			rb = r_axi-threadIdx.x-1;
			zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
			si2 =cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(4*numbdry/8) && blockIdx.x<(5*numbdry/8))
        {
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
       	        	ra = r_axi-threadIdx.x;
               		za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
			si1 = cu_result[ra*mesh+za];

                	rb = r_axi-threadIdx.x-1;
                	zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
                	si2 = cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(5*numbdry/8) && blockIdx.x<(6*numbdry/8))
	{
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
			za = z_axi-threadIdx.x;
			ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
			si1 = cu_result[ra*mesh+za];
			
			zb = z_axi-threadIdx.x-1;
			rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
			si2 = cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(6*numbdry/8) && blockIdx.x<(7*numbdry/8))
        {
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
			za = z_axi-threadIdx.x;
               		ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
			si1 = cu_result[ra*mesh+za];

                	zb = z_axi-threadIdx.x-1;
                	rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
                	si2 = cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else
	{
		if (1)//(threadIdx.x <= num[blockIdx.x])
		{
			ra = r_axi+threadIdx.x;
			za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
			si1 = cu_result[ra*mesh+za];

			rb = r_axi+threadIdx.x+1;
			zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
			si2 = cu_result[rb*mesh+zb];

			for(i=0;i<num[0];i++)
			{
				flux = nphi[i];
				if ((flux-si1)>0 && (flux-si2)<0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ; 
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
}

/*


__global__ void flux_function_integral_constrain(float* cu_result,float* current_representation,float* r_press,float* z_press,float* cu_press,float* test_right)
{
	float x;
	int i,j;

	x = [threadIdx.x];

	for (j=0;j<num_knot-1;j++)  // integral value on each knot segment
	{
		if ((x >= y[j]) && (x <= y[j+1]))
		{
			integral[threadIdx.x*(num_knot-1)*4 + j*4] = (y[j+1] - x);
			integral[threadIdx.x*(num_knot-1)*4 + j*4 + 1] = 0.5*(y[j+1]*y[j+1]-x*x);
			integral[threadIdx.x*(num_knot-1)*4 + j*4 + 2] = cosh(tension*y[j+1]) - cosh(tension*x);
			integral[threadIdx.x*(num_knot-1)*4 + j*4 + 3] = sinh(tension*y[j+1]) - sinh(tension*x);

			for (i=j+1;i<num_knot-1;i++)
			{
				integral[threadIdx.x*(num_knot-1)*4 + i*4] = (y[i+1] - y[i]);
				integral[threadIdx.x*(num_knot-1)*4 + i*4 + 1] = 0.5*(y[i+1]*y[i+1]-y[i]*y[i]);
				integral[threadIdx.x*(num_knot-1)*4 + i*4 + 2] = cosh(tension*y[i+1]) - cosh(tension*y[i]);
				integral[threadIdx.x*(num_knot-1)*4 + i*4 + 3] = sinh(tension*y[i+1]) - sinh(tension*y[i]);
			}
		}
	}	
}

__global__ void jphi_constrain_1(float* cu_result,float* location,float* tani,float* coti,float* xpoint,int* num,float* nphi,float* xy_surface)
{
	// search flux surface for j_phi constrain
	int ra,rb,za,zb,r_axi,z_axi,i,N1,N2;
	float si1,si2,flux,dflux,aflux,factor;

	z_axi = int(location[4]);
	r_axi = int(location[5])+68;

	aflux = xpoint[4];
	dflux = xpoint[11];
	
	__syncthreads();

	if (blockIdx.x<(1+numbdry/8))
	{
		if(threadIdx.x <= num[blockIdx.x])
		{
			ra = r_axi+threadIdx.x;
			za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
			si1 = cu_result[ra*mesh1+za];

			rb = r_axi+threadIdx.x+1;
			zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();		
	}
	else if (blockIdx.x>=(1+numbdry/8) && blockIdx.x<(2*numbdry/8))
	{
		if(threadIdx.x <= num[blockIdx.x])
		{
			za = z_axi+threadIdx.x;
			ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
			si1 = cu_result[ra*mesh1+za];

			zb = z_axi+threadIdx.x+1;
			rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(2*numbdry/8) && blockIdx.x<(3*numbdry/8))
        {
		if(threadIdx.x <= num[blockIdx.x])
		{
                	za = z_axi+threadIdx.x;
       	      		ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
			si1 = cu_result[ra*mesh1+za];

           	        zb = z_axi+threadIdx.x+1;
        	        rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
               		si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(3*numbdry/8) && blockIdx.x<(4*numbdry/8))
	{
		if(threadIdx.x <= num[blockIdx.x])
		{
			ra = r_axi-threadIdx.x;
			za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
			si1 = cu_result[ra*mesh1+za];

			rb = r_axi-threadIdx.x-1;
			zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
			si2 =cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(4*numbdry/8) && blockIdx.x<(5*numbdry/8))
        {
		if(threadIdx.x <= num[blockIdx.x])
		{
       	        	ra = r_axi-threadIdx.x;
               		za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
			si1 = cu_result[ra*mesh1+za];

                	rb = r_axi-threadIdx.x-1;
                	zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
                	si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(5*numbdry/8) && blockIdx.x<(6*numbdry/8))
	{
		if(threadIdx.x <= num[blockIdx.x])
		{
			za = z_axi-threadIdx.x;
			ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
			si1 = cu_result[ra*mesh1+za];
			
			zb = z_axi-threadIdx.x-1;
			rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else if (blockIdx.x>=(6*numbdry/8) && blockIdx.x<(7*numbdry/8))
        {
		if(threadIdx.x <= num[blockIdx.x])
		{
			za = z_axi-threadIdx.x;
               		ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
			si1 = cu_result[ra*mesh1+za];

                	zb = z_axi-threadIdx.x-1;
                	rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
                	si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        		xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
	else
	{
		if(threadIdx.x <= num[blockIdx.x])
		{
			ra = r_axi+threadIdx.x;
			za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
			si1 = cu_result[ra*mesh1+za];

			rb = r_axi+threadIdx.x+1;
			zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
			si2 = cu_result[rb*mesh1+zb];

			for(i=0;i<num[0];i++)
			{
				flux = aflux - nphi[i]*dflux;
				if ((flux-si1)>=0 && (flux-si2)<=0)
				{
					factor = (flux-si1)/(si2-si1);
					xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
					xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ; 
				}
			}
			__syncthreads();
		}
		__syncthreads();
	}
}

__global__ void jphi_constrain_2(float* nphi,float* xy_surface,int* num_knot,float* y,float* x)
{
	int j,i;
	float r,z,dl,a,b,jphi,phi;
	__shared__ float left[numbdry][8];
	__shared__ float sum[numbdry];

	num_knot = num[0];
	phi = nphi[blockIdx.x];
	r = xy_surface[blockIdx.x*2*numbdry+threadIdx.x];
	z = xy_surface[(blockIdx.x*2+1)*numbdry+threadIdx.x];

	if (threadIdx.x < (numbdry-1))
	{
		a = xy_surface[blockIdx.x*2*numbdry+threadIdx.x+1];
		b = xy_surface[(blockIdx.x*2+1)*numbdry+threadIdx.x+1];
		sum[threadIdx.x] = sqrt((r-a)*(r-a)+(z-b)*(z-b));
	}
	else if (threadIdx.x == (numbdry-1))
	{
		a = xy_surface[blockIdx.x*2*numbdry];
		b = xy_surface[(blockIdx.x*2+1)*numbdry];
		sum[threadIdx.x] = sqrt((r-a)*(r-a)+(z-b)*(z-b)); 
	}
	__syncthreads();

	if (threadIdx.x == 0)
	{
		dl = (sum[threadIdx.x] + sum[numbdry-1])/2;
	}
	else if (threadIdx.x > 0)
	{	
		dl = (sum[threadIdx.x] + sum[threadIdx.x-1])/2;
	}

	for (j=0;j<num_knot[0]-1;j++)
	{
		if ((phi-y[j])*(phi-y[j+1])<0)
		{
			left[threadIdx.x][0] = r * dl;
			left[threadIdx.x][1] = r * phi * dl;
			left[threadIdx.x][2] = r * sinh(tension*phi) * dl;
			left[threadIdx.x][3] = r * cosh(tension*phi) * dl;
			left[threadIdx.x][4] = 1/r * dl;
			left[threadIdx.x][5] = 1/r * phi * dl;
			left[threadIdx.x][6] = 1/r * sinh(tension*phi) * dl;
			left[threadIdx.x][7] = 1/r * cosh(tension*phi) * dl;
			i = j;
		}
	}	
	__syncthreads();

	for(unsigned int s=64;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<numbdry)
                {
                        left[threadIdx.x][0] += left[k][0];
			left[threadIdx.x][1] += left[k][1];
			left[threadIdx.x][2] += left[k][2];
			left[threadIdx.x][3] += left[k][3];
			left[threadIdx.x][4] += left[k][4];
			left[threadIdx.x][5] += left[k][5];
			left[threadIdx.x][6] += left[k][6];
			left[threadIdx.x][7] += left[k][7];
			sum[threadIdx.x] += sum[k];
                }
		__syncthreads();
	}

	if (threadIdx.x < 8)
	{
		x[blockIdx.x*(num_knot[0]-1)*8 + i*8 + threadIdx.x] = left[0][threadIdx.x]/sum[0];	
	}
}


*/
