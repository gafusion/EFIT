/*
#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<cutil.h>
#include<cublas_v2.h>
#include <sys/time.h>
#include "xlocate.h"
#include "macro.h"

extern __constant__ float cu_radix[mesh];
extern __constant__ int cu_num_block[mesh];
extern __constant__ int cu_blo_accu[mesh];
extern __constant__ int cu_num_se[1949];
extern __constant__ float cu_currentInver[mesh]; 
*/
__global__ void expansion_finder(float *cu_result,float* leas_flux)
{
	float br,bz,c,d,invr;
	__shared__ float sDiff[4][63];
	__shared__ float position[63];
	position[threadIdx.x] = float(threadIdx.x);
	sDiff[3][threadIdx.x] = 9999;
	__syncthreads();
	if (blockIdx.x<122)
	{
		invr=1/cu_radix[blockIdx.x+23];
		sDiff[0][threadIdx.x] = cu_result[(blockIdx.x+22)*mesh1+threadIdx.x+7];
		sDiff[1][threadIdx.x] = cu_result[(blockIdx.x+23)*mesh1+threadIdx.x+7];
		sDiff[2][threadIdx.x] = cu_result[(blockIdx.x+24)*mesh1+threadIdx.x+7];
		bz = invr*(sDiff[2][threadIdx.x]-sDiff[0][threadIdx.x])*INVR;        // bz at blockIdx.x+23
		__syncthreads();
		if(threadIdx.x>0&&threadIdx.x<62)
		{
			br = -invr*(sDiff[1][threadIdx.x+1]-sDiff[1][threadIdx.x-1])*INVZ;  // br at threadIdx.x+7+1
			sDiff[3][threadIdx.x-1] = br*br + bz*bz;   // bp at r = blockIdx.x+23, z = br at threadIdx.x+7+1 
		}
		__syncthreads();
		for(unsigned int s=32;s>0;s>>=1)
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<61) 
			{
				c=sDiff[3][threadIdx.x];
				d=sDiff[3][k];
				if(c>d)
				{
					sDiff[3][threadIdx.x]=d;
					position[threadIdx.x] = position[k];
				}
			}
		}
		leas_flux[blockIdx.x] = sDiff[3][0];
		leas_flux[122+blockIdx.x] = position[0]+1+7;
//		printf("%f  %f %f\n",leas_flux[blockIdx.x],leas_flux[(blockIdx.x)+25],blockIdx.x);
//		printf("%f \n",leas_flux[(blockIdx.x)+25]);
	}
	else if (blockIdx.x<244)             
	{
		invr=1/cu_radix[blockIdx.x+23-122];
		sDiff[0][threadIdx.x] = cu_result[(blockIdx.x-122+22)*mesh1+threadIdx.x+187];
		sDiff[1][threadIdx.x] = cu_result[(blockIdx.x-122+23)*mesh1+threadIdx.x+187];
		sDiff[2][threadIdx.x] = cu_result[(blockIdx.x-122+24)*mesh1+threadIdx.x+187];
		bz = invr*(sDiff[2][threadIdx.x]-sDiff[0][threadIdx.x])*INVR;   // br at blockIdx.x+23
		__syncthreads();
		if(threadIdx.x>0&&threadIdx.x<62)
		{
			br = -invr*(sDiff[1][threadIdx.x+1]-sDiff[1][threadIdx.x-1])*INVZ; // br at threadIdx.x+187+1
			sDiff[3][threadIdx.x-1] = br*br + bz*bz; // bp at r = blockIdx.x+23, z = br at threadIdx.x+187+1 
		}
		__syncthreads();
		for(unsigned int s=32;s>0;s>>=1)
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<61) 
			{
				c=sDiff[3][threadIdx.x];
				d=sDiff[3][k];
				if(c>d)
				{
					sDiff[3][threadIdx.x]=d;
					position[threadIdx.x] = position[k];
				}
			}
		}
		leas_flux[2000+blockIdx.x-122] = sDiff[3][0];
		leas_flux[2000+blockIdx.x] = position[0]+1+187;
//		printf("%f %f %f\n",leas_flux[2000+(blockIdx.x-25)],leas_flux[2000+(blockIdx.x-25)+25],(blockIdx.x-25));
//		printf("%f \n",leas_flux[2000+(blockIdx.x-25)+25]);
	}
//	leas_flux[]
	else	
	{
		sDiff[0][threadIdx.x] = cu_result[(blockIdx.x-244+68)*mesh1+threadIdx.x+99];
		__syncthreads();
		for(unsigned int s=32;s>0;s>>=1) 
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<57) 
			{
				c=sDiff[0][threadIdx.x];
				d=sDiff[0][k];
				if(c<d)     // for plasma current >0
//				if(c>d)		// for plasma current < 0
				{
					sDiff[0][threadIdx.x] = d;
					position[threadIdx.x] = position[k];	
				}
			}
		}
		// write result for this block to global mem
		if (threadIdx.x == 0)
		{ 
			leas_flux[3000+blockIdx.x-244] = sDiff[0][0];
			leas_flux[3000+blockIdx.x-244+122] = position[0]+99;
		}	
	}	
}
/*

 |e  f| |dr| = |-d|       (Bz equation)
 |b  c| |dz|   |-a|       (Br equation)

and solved using determinants:

 dr = (af - cd)/(ec - bf)
 dz = (bd - ea)/(ec - bf)

where  a = Br_expansion, b = dBr/dr, c = dBr/dz
       d = Bz_expansion, e = dBz/dr, f = dBz/dz

*/
__global__ void xpoint_brbz(float* leas_flux, float* cu_result)
{
	int rID, zID;
	float br,bz,invr;
	__shared__ float tempbr[5];
	__shared__ float tempbz[5];
	if(blockIdx.z == 0)
	{
		rID = int(leas_flux[1])+23-1+blockIdx.x;
		zID = int(leas_flux[0])-1+blockIdx.y;
		invr = 1/cu_radix[rID];

		br = -invr*((cu_result[(rID+threadIdx.x)*mesh1+(zID+1+threadIdx.y)])-(cu_result[(rID+threadIdx.x)*mesh1+(zID-1+threadIdx.y)]))*INVZ;
		bz = invr*((cu_result[(rID+1+threadIdx.x)*mesh1+(zID+threadIdx.y)])-(cu_result[(rID-1+threadIdx.x)*mesh1+(zID+threadIdx.y)]))*INVR;

		tempbr[threadIdx.x*2+threadIdx.y] = br;
		tempbz[threadIdx.x*2+threadIdx.y] = bz;
		
		__syncthreads();

		tempbr[4] = (tempbr[0]+tempbr[1]+tempbr[2]+tempbr[3])/4;
		tempbz[4] = (tempbz[0]+tempbz[1]+tempbz[2]+tempbz[3])/4;

		leas_flux[4000+100*blockIdx.z+blockIdx.y*2+blockIdx.x] = tempbr[4];
		leas_flux[5000+100*blockIdx.z+blockIdx.y*2+blockIdx.x] = tempbz[4];
	}
	if(blockIdx.z == 1)
	{
		rID = int(leas_flux[3])+23-1+blockIdx.x;
		zID = int(leas_flux[2])-1+blockIdx.y;
		invr = 1/cu_radix[rID];

		br = -invr*((cu_result[(rID+threadIdx.x)*mesh1+(zID+1+threadIdx.y)])-(cu_result[(rID+threadIdx.x)*mesh1+(zID-1+threadIdx.y)]))*INVZ;
		bz = invr*((cu_result[(rID+1+threadIdx.x)*mesh1+(zID+threadIdx.y)])-(cu_result[(rID-1+threadIdx.x)*mesh1+(zID+threadIdx.y)]))*INVR;

		tempbr[threadIdx.x*2+threadIdx.y] = br;
		tempbz[threadIdx.x*2+threadIdx.y] = bz;
		
		__syncthreads();

		tempbr[4] = (tempbr[0]+tempbr[1]+tempbr[2]+tempbr[3])/4;
		tempbz[4] = (tempbz[0]+tempbz[1]+tempbz[2]+tempbz[3])/4;

		leas_flux[4000+100*blockIdx.z+blockIdx.y*2+blockIdx.x] = tempbr[4];
		leas_flux[5000+100*blockIdx.z+blockIdx.y*2+blockIdx.x] = tempbz[4];
	}
	if(blockIdx.z == 2)
	{
		rID = int(leas_flux[5])+68-1+blockIdx.x;
		zID = int(leas_flux[4])-1+blockIdx.y;
		invr = 1/cu_radix[rID];

		br = -invr*((cu_result[(rID+threadIdx.x)*mesh1+(zID+1+threadIdx.y)])-(cu_result[(rID+threadIdx.x)*mesh1+(zID-1+threadIdx.y)]))*INVZ;
		bz = invr*((cu_result[(rID+1+threadIdx.x)*mesh1+(zID+threadIdx.y)])-(cu_result[(rID-1+threadIdx.x)*mesh1+(zID+threadIdx.y)]))*INVR;

		tempbr[threadIdx.x*2+threadIdx.y] = br;
		tempbz[threadIdx.x*2+threadIdx.y] = bz;
		
		__syncthreads();

		tempbr[4] = (tempbr[0]+tempbr[1]+tempbr[2]+tempbr[3])/4;
		tempbz[4] = (tempbz[0]+tempbz[1]+tempbz[2]+tempbz[3])/4;

		leas_flux[4000+100*blockIdx.z+blockIdx.y*2+blockIdx.x] = tempbr[4];
		leas_flux[5000+100*blockIdx.z+blockIdx.y*2+blockIdx.x] = tempbz[4];
	}
}
/*
__global__ void xpoint_brbz1(float* leas_flux, float* cu_result)
{
	int rID, zID, dr, dz;
	float br,bz,invr;
	
	if(blockIdx.x == 0)
	{
		if(threadIdx.y == 0)
		{
			dr = threadIdx.x*2 - 1;
			dz = 0;
			rID = int(leas_flux[1])+23+dr;
			zID = int(leas_flux[0])+dz;
			invr = 1/float(cu_radix[rID]);

			br = -invr*(float(cu_result[rID*mesh1+(zID+1)])-float(cu_result[(rID)*mesh1+(zID-1)]))*INVZ;
			bz = invr*(float(cu_result[(rID+1)*mesh1+(zID)])-float(cu_result[(rID-1)*mesh1+(zID)]))*INVR;
			__syncthreads();

			leas_flux[4000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = br;
			leas_flux[5000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = bz;
		}
		if(threadIdx.y == 1)
		{
			dz = threadIdx.x*2 - 1;
			dr = 0;
			rID = int(leas_flux[1])+23+dr;
			zID = int(leas_flux[0])+dz;
			invr = 1/float(cu_radix[rID]);

			br = -invr*(float(cu_result[(rID)*mesh1+(zID+1)])-float(cu_result[(rID)*mesh1+(zID-1)]))*INVZ;
			bz = invr*(float(cu_result[(rID+1)*mesh1+(zID)])-float(cu_result[(rID-1)*mesh1+(zID)]))*INVR;
			__syncthreads();

			leas_flux[4000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = br;
			leas_flux[5000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = bz;
		}
	}
	if(blockIdx.x == 1)
	{
		if(threadIdx.y == 0)
		{
			dr = threadIdx.x*2 - 1;
			dz = 0;
			rID = int(leas_flux[2])+23+dr;
			zID = int(leas_flux[3])+dz;
			invr = 1/float(cu_radix[rID]);

			br = -invr*(float(cu_result[rID*mesh1+(zID+1)])-float(cu_result[(rID)*mesh1+(zID-1)]))*INVZ;
			bz = invr*(float(cu_result[(rID+1)*mesh1+(zID)])-float(cu_result[(rID-1)*mesh1+(zID)]))*INVR;
			__syncthreads();

			leas_flux[4000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = br;
			leas_flux[5000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = bz;
		}
		if(threadIdx.y == 1)
		{
			dz = threadIdx.x*2 - 1;
			dr = 0;
			rID = int(leas_flux[2])+23+dr;
			zID = int(leas_flux[3])+dz;
			invr = 1/float(cu_radix[rID]);

			br = -invr*(float(cu_result[(rID)*mesh1+(zID+1)])-float(cu_result[(rID)*mesh1+(zID-1)]))*INVZ;
			bz = invr*(float(cu_result[(rID+1)*mesh1+(zID)])-float(cu_result[(rID-1)*mesh1+(zID)]))*INVR;
			__syncthreads();

			leas_flux[4000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = br;
			leas_flux[5000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = bz;
		}
	}
	if(blockIdx.x == 2)
	{
		rID = int(leas_flux[4])+68-1+blockIdx.x;
		zID = int(leas_flux[5])-1+blockIdx.y;

		if(threadIdx.y == 0)
		{
			dr = threadIdx.x*2 - 1;
			dz = 0;
			rID = int(leas_flux[4])+23+dr;
			zID = int(leas_flux[5])+dz;
			invr = 1/float(cu_radix[rID]);

			br = -invr*(float(cu_result[rID*mesh1+(zID+1)])-float(cu_result[(rID)*mesh1+(zID-1)]))*INVZ;
			bz = invr*(float(cu_result[(rID+1)*mesh1+(zID)])-float(cu_result[(rID-1)*mesh1+(zID)]))*INVR;
			__syncthreads();

			leas_flux[4000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = br;
			leas_flux[5000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = bz;
		}
		if(threadIdx.y == 1)
		{
			dz = threadIdx.x*2 - 1;
			dr = 0;
			rID = int(leas_flux[4])+23+dr;
			zID = int(leas_flux[5])+dz;
			invr = 1/float(cu_radix[rID]);

			br = -invr*(float(cu_result[(rID)*mesh1+(zID+1)])-float(cu_result[(rID)*mesh1+(zID-1)]))*INVZ;
			bz = invr*(float(cu_result[(rID+1)*mesh1+(zID)])-float(cu_result[(rID-1)*mesh1+(zID)]))*INVR;
			__syncthreads();

			leas_flux[4000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = br;
			leas_flux[5000+100*blockIdx.x+threadIdx.y*2+threadIdx.x] = bz;
		}
	}
}
*/
__global__ void post_xpoint_brbz(float* leas_flux)
{
	__shared__ float br[4];
	__shared__ float bz[4];

	br[threadIdx.x] = leas_flux[4000+100*blockIdx.x+threadIdx.x];
	bz[threadIdx.x] = leas_flux[5000+100*blockIdx.x+threadIdx.x];

	float temp;

	if (threadIdx.x ==0) //br
	{
		temp = (br[0]+br[1]+br[2]+br[3])/4;	
	}
	if (threadIdx.x ==1) //bz
	{
		temp = (bz[0]+bz[1]+bz[2]+bz[3])/4;	
	}
	if (threadIdx.x ==2) //brr
	{
		temp = ((br[3]+br[1])-(br[2]+br[0]))/(2*delR);
	}
	if (threadIdx.x ==3) //brz
	{
		temp = ((br[3]+br[2])-(br[0]+br[1]))/(2*delZ);
	}
	if (threadIdx.x ==4) //bzr
	{
		temp = ((bz[3]+bz[1])-(bz[2]+bz[0]))/(2*delR);
	}
	if (threadIdx.x ==5) //bzz
	{
		temp = ((bz[3]+bz[2])-(bz[0]+bz[1]))/(2*delZ);
	}

	leas_flux[1000+blockIdx.x*6+threadIdx.x] = temp;
}


__global__ void xpoint_pos(float* leas_flux, float* cu_result, float* xpoint, float* cu_xlocation)
{
	float ebr, ebz, ebrr, ebrz, ebzr, ebzz, dr, dz, er, ez, xr, xz, xflux;
	int rID, zID;
	if(blockIdx.x == 0)
	{
		zID = leas_flux[0];
		rID = leas_flux[1];
		ebr = leas_flux[1000+0];
		ebz = leas_flux[1000+1];
		ebrr = leas_flux[1000+2];
		ebrz = leas_flux[1000+3];
		ebzr = leas_flux[1000+4];
		ebzz = leas_flux[1000+5];
//		printf("\ndown br %f,bz %f,brr %f,brz %f,bzr %f,bzz %f\n", ebr, ebz, ebrr, ebrz, ebzr, ebzz);
		dr = (ebr*ebzz - ebrz*ebz)/(ebzr*ebrz - ebrr*ebzz);
	 	dz = (ebrr*ebz - ebzr*ebr)/(ebzr*ebrz - ebrr*ebzz);
		er = igridr+delR+delR*(rID+23);
		ez = igridz+delZ+delZ*zID;
		xr = er + dr;
		xz = ez + dz;
//		printf("\ndown r er %f,ez %f,dr %f,dz %f,xr %f,xz %f\n", er, ez, dr, dz, xr, xz);
//		printf("\ndown d  %f, %f \n", dr ,dz);
//		printf("\ndown x %f, %f \n", xr ,xz);
		xflux = cu_result[(rID+23)*mesh1+zID];
//		printf("\ndown flux old %f \n", xflux);
		xflux += cu_radix[rID+23]*(-dz*ebr+ebz*dr);
		printf("\ndown flux new %f \n", xflux);
		xpoint[0] = xflux;
		xpoint[5] = xr ;	
		xpoint[6] = xz ;
//		cu_xlocation[0] = xr;
//		cu_xlocation[1] = xz;
	}
	if(blockIdx.x == 1)
	{
		zID = leas_flux[2];
		rID = leas_flux[3];
		ebr = leas_flux[1000+6];
		ebz = leas_flux[1000+7];
		ebrr = leas_flux[1000+8];
		ebrz = leas_flux[1000+9];
		ebzr = leas_flux[1000+10];
		ebzz = leas_flux[1000+11];
//		printf("\nup br %f,bz %f,brr %f,brz %f,bzr %f,bzz %f\n", ebr, ebz, ebrr, ebrz, ebzr, ebzz);
		dr = (ebr*ebzz - ebrz*ebz)/(ebzr*ebrz - ebrr*ebzz);
	 	dz = (ebrr*ebz - ebzr*ebr)/(ebzr*ebrz - ebrr*ebzz);
		er = igridr+delR+delR*(rID+23);
		ez = igridz+delZ+delZ*zID;
		xr = er + dr;
		xz = ez + dz;
//		printf("\nup r er %f,ez %f,dr %f,dz %f,xr %f,xz %f\n", er, ez, dr, dz, xr, xz);
//		printf("\nup  %f, %f \n", er ,ez);
//		printf("\nup d %f, %f \n", dr ,dz);
//		printf("\nup x %f, %f \n", xr ,xz);
		xflux = cu_result[(rID+23)*mesh1+zID];
//		printf("\nup flux old %f \n", xflux);
		xflux += cu_radix[rID+23]*(-dz*ebr+ebz*dr);
		xpoint[2] = xflux;
		xpoint[7] = xr ;	
		xpoint[8] = xz ;	
		printf("\nup flux new %f \n", xflux);
//		cu_xlocation[2] = xr;
//		cu_xlocation[3] = xz;
	}
	if(blockIdx.x == 2)
	{
		zID = leas_flux[4];
		rID = leas_flux[5];
		ebr = leas_flux[1000+12];
		ebz = leas_flux[1000+13];
		ebrr = leas_flux[1000+14];
		ebrz = leas_flux[1000+15];
		ebzr = leas_flux[1000+16];
		ebzz = leas_flux[1000+17];
//		printf("\nmid br %f,bz %f,brr %f,brz %f,bzr %f,bzz %f\n", ebr, ebz, ebrr, ebrz, ebzr, ebzz);
		dr = (ebr*ebzz - ebrz*ebz)/(ebzr*ebrz - ebrr*ebzz);
	 	dz = (ebrr*ebz - ebzr*ebr)/(ebzr*ebrz - ebrr*ebzz);
		er = igridr+delR+delR*(rID+68);
		ez = igridz+delZ+delZ*zID;
		xr = er + dr;
		xz = ez + dz;
//		printf("\nmid r r er %f,ez %f,dr %f,dz %f,xr %f,xz %f\n", er, ez, dr, dz, xr, xz);
//		printf("\nmid  %f, %f \n", er ,ez);
//		printf("\nmid  %f, %f \n", dr ,dz);
//		printf("\nmid  %f, %f \n", xr ,xz);
		xflux = cu_result[(rID+68)*mesh1+zID];
//		printf("\nmid flux old %f \n", xflux);
		xflux += cu_radix[rID+68]*(-dz*ebr+ebz*dr);
		xpoint[4] = xflux;
		xpoint[9] = er ;	
		xpoint[10] = xz ;	
		printf("\nmid flux new %f \n", xflux);
		cu_xlocation[4] = float(zID);
		cu_xlocation[5] = float(rID);
	}
}


