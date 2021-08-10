__global__ void searchbdry(float* cu_result,float* leas_flux,float* xpoint,float* cu_bound_flux,float* tani,float* coti,float* xybdry,int* maxgridbdry)
{
	int ra,rb,za,zb,r_axi,z_axi;
	float si1,si2,dif,bdry,factor;
	__shared__ float xbdry[128];
	__shared__ float ybdry[128];
	__shared__ int nbdry[128];
	nbdry[threadIdx.x] = 128;
	z_axi = int(leas_flux[4]);
	r_axi = int(leas_flux[5])+68;
	bdry = cu_bound_flux[0];
//	printf("  %f  ",bdry);
//	printf("%d  %d  \n",r_axi,z_axi);
//	printf("%15.7e  %15.7e\n",tani[blockIdx.x],coti[blockIdx.x]);
	__syncthreads();
	if (blockIdx.x<(1+numbdry/8))
	{
		ra = r_axi+threadIdx.x;
		za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
		rb = r_axi+threadIdx.x+1;
		zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
		si1 = cu_result[ra*mesh1+za];
		si2 = cu_result[rb*mesh1+zb];
		dif = (si1-bdry)*(si2-bdry);
		if (dif<=0)
		{
			factor = (bdry-si1)/(si2-si1);
			xbdry[threadIdx.x] = igridr+delR+delR*ra+factor*delR;
			ybdry[threadIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
			nbdry[threadIdx.x] = threadIdx.x;		
		}
		__syncthreads();
		for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
				int c = nbdry[threadIdx.x];
				int d = nbdry[k];
				if (c>d)
				{
					nbdry[threadIdx.x] = d; 
				}
				__syncthreads();
			}	
		}
//		xybdry[blockIdx.x] = xbdry[nbdry[0]];
//		xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
	}		
	else if (blockIdx.x>=(1+numbdry/8) && blockIdx.x<(2*numbdry/8))
	{
		za = z_axi+threadIdx.x;
		ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
		zb = z_axi+threadIdx.x+1;
		rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
		si1 = cu_result[ra*mesh1+za];
		si2 = cu_result[rb*mesh1+zb];
		dif = (si1-bdry)*(si2-bdry);
		if (dif<=0)
		{
			factor = (bdry-si1)/(si2-si1);
			xbdry[threadIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
			ybdry[threadIdx.x] = igridz+delZ+delZ*za+factor*delZ;
		       	nbdry[threadIdx.x] = threadIdx.x;	
		}
		__syncthreads();
		for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
				int c = nbdry[threadIdx.x];
				int d = nbdry[k];
				if (c>d)
				{
					nbdry[threadIdx.x] = d; 
				}
			}
			__syncthreads();	
		}
//		xybdry[blockIdx.x] = xbdry[nbdry[0]];
//		xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
		__syncthreads();
		if (nbdry[0] > (235-z_axi))
		{
			nbdry[0] = int((xpoint[8]-igridz)/delZ) - 1 - z_axi;
			xbdry[nbdry[0]] = xpoint[7];
			ybdry[nbdry[0]] = xpoint[8];
		}
	}
	else if (blockIdx.x>=(2*numbdry/8) && blockIdx.x<(3*numbdry/8))
        {
                za = z_axi+threadIdx.x;
                ra = r_axi+int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
                zb = z_axi+threadIdx.x+1;
                rb = r_axi+int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
                si1 = cu_result[ra*mesh1+za];
                si2 = cu_result[rb*mesh1+zb];
                dif = (si1-bdry)*(si2-bdry);
                if (dif<=0)
                {
                        factor = (bdry-si1)/(si2-si1);
                        xbdry[threadIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        ybdry[threadIdx.x] = igridz+delZ+delZ*za+factor*delZ;
                        nbdry[threadIdx.x] = threadIdx.x;
                }
		__syncthreads();
                for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
                                int c = nbdry[threadIdx.x];
                                int d = nbdry[k];
                                if (c>d)
                                {
                                        nbdry[threadIdx.x] = d;
                                }
                        }
			__syncthreads();
                }
		__syncthreads();
//              xybdry[blockIdx.x] = xbdry[nbdry[0]];
//              xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
                if (nbdry[0] > (235-z_axi))
		{
			nbdry[0] = int((xpoint[8]-igridz)/delZ) - 1 - z_axi;
			xbdry[nbdry[0]] = xpoint[7];
			ybdry[nbdry[0]] = xpoint[8];
//			printf("add up x-point at search\n");
		}
        }
	else if (blockIdx.x>=(3*numbdry/8) && blockIdx.x<(4*numbdry/8))
	{
		ra = r_axi-threadIdx.x;
		za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ);
		rb = r_axi-threadIdx.x-1;
		zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ);
		si1 = cu_result[ra*mesh1+za];
		si2 = cu_result[rb*mesh1+zb];
		dif = (si1-bdry)*(si2-bdry);
		if (dif<=0)
		{
			factor = (bdry-si1)/(si2-si1);
			xbdry[threadIdx.x] = igridr+delR+delR*ra-factor*delR;
			ybdry[threadIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
			nbdry[threadIdx.x] = threadIdx.x;
		}
		__syncthreads();
		for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
				int c = nbdry[threadIdx.x];
				int d = nbdry[k];
				if (c>d)
				{
					nbdry[threadIdx.x] = d; 
				}
			}
			__syncthreads();
		}
//		xybdry[blockIdx.x] = xbdry[nbdry[0]];
//		xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
	}
	else if (blockIdx.x>=(4*numbdry/8) && blockIdx.x<(5*numbdry/8))
        {
                ra = r_axi-threadIdx.x;
                za = z_axi-int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
                rb = r_axi-threadIdx.x-1;
                zb = z_axi-int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
                si1 = cu_result[ra*mesh1+za];
                si2 = cu_result[rb*mesh1+zb];
                dif = (si1-bdry)*(si2-bdry);
                if (dif<=0)
                {
                        factor = (bdry-si1)/(si2-si1);
                        xbdry[threadIdx.x] = igridr+delR+delR*ra-factor*delR;
                        ybdry[threadIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
                        nbdry[threadIdx.x] = threadIdx.x;
                }
		__syncthreads();
                for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
                                int c = nbdry[threadIdx.x];
                                int d = nbdry[k];
                                if (c>d)
                                {
                                        nbdry[threadIdx.x] = d;
                                }
                        }
			__syncthreads();
                }
//              xybdry[blockIdx.x] = xbdry[nbdry[0]];
//              //              xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
        }
	else if (blockIdx.x>=(5*numbdry/8) && blockIdx.x<(6*numbdry/8))
	{
		za = z_axi-threadIdx.x;
		ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR)-1;
		zb = z_axi-threadIdx.x-1;
		rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR)-1;
		si1 = cu_result[ra*mesh1+za];
		si2 = cu_result[rb*mesh1+zb];
		dif = (si1-bdry)*(si2-bdry);
		if (dif<=0)
		{
			factor = (bdry-si1)/(si2-si1);
			xbdry[threadIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
			ybdry[threadIdx.x] = igridz+delZ+delZ*za-factor*delZ;
			nbdry[threadIdx.x] = threadIdx.x;
		}
		__syncthreads();
		for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
				int c = nbdry[threadIdx.x];
				int d = nbdry[k];
				if (c>d)
				{
					nbdry[threadIdx.x] = d; 
				}
			}
			__syncthreads();	
		}
//		xybdry[blockIdx.x] = xbdry[nbdry[0]];
//		xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
		__syncthreads();
		if (nbdry[0] > (z_axi-18))
		{
			nbdry[0] = z_axi - int((xpoint[6]-igridz)/delZ);
			xbdry[nbdry[0]] = xpoint[5];
			ybdry[nbdry[0]] = xpoint[6];
		}
	}
	else if (blockIdx.x>=(6*numbdry/8) && blockIdx.x<(7*numbdry/8))
        {
                za = z_axi-threadIdx.x;
                ra = r_axi-int(threadIdx.x*delZ*coti[blockIdx.x]/delR);
                zb = z_axi-threadIdx.x-1;
                rb = r_axi-int((threadIdx.x+1)*delZ*coti[blockIdx.x]/delR);
                si1 = cu_result[ra*mesh1+za];
                si2 = cu_result[rb*mesh1+zb];
                dif = (si1-bdry)*(si2-bdry);
                if (dif<=0)
                {
                        factor = (bdry-si1)/(si2-si1);
                        xbdry[threadIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        ybdry[threadIdx.x] = igridz+delZ+delZ*za-factor*delZ;
                        nbdry[threadIdx.x] = threadIdx.x;
                }
		__syncthreads();
                for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
                                int c = nbdry[threadIdx.x];
                                int d = nbdry[k];
                                if (c>d)
                                {
                                        nbdry[threadIdx.x] = d;
                                }
                        }
			__syncthreads();
                }
//              xybdry[blockIdx.x] = xbdry[nbdry[0]];
//              //              xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
		__syncthreads();
		if (nbdry[0] > (z_axi-18))
		{
			nbdry[0] = z_axi - int((xpoint[6]-igridz)/delZ);
			xbdry[nbdry[0]] = xpoint[5];
			ybdry[nbdry[0]] = xpoint[6];
		}
	}
	else
	{
		ra = r_axi+threadIdx.x;
		za = z_axi+int(threadIdx.x*delR*tani[blockIdx.x]/delZ)-1;
		rb = r_axi+threadIdx.x+1;
		zb = z_axi+int((threadIdx.x+1)*delR*tani[blockIdx.x]/delZ)-1;
		si1 = cu_result[ra*mesh1+za];
		si2 = cu_result[rb*mesh1+zb];
		dif = (si1-bdry)*(si2-bdry);
		if (dif<=0)
		{
			factor = (bdry-si1)/(si2-si1);
			xbdry[threadIdx.x] = igridr+delR+delR*ra+factor*delR;
			ybdry[threadIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ; 
			nbdry[threadIdx.x] = threadIdx.x;
		}
		__syncthreads();
		for(unsigned int s=64;s>0;s>>=1)
                {
                        int k=threadIdx.x+s;
                        if (threadIdx.x<s&&k<128)
                        {
				int c = nbdry[threadIdx.x];
				int d = nbdry[k];
				if (c>d)
				{
					nbdry[threadIdx.x] = d; 
				}
			}
			__syncthreads();	
		}
//		xybdry[blockIdx.x] = xbdry[nbdry[0]];
//		xybdry[blockIdx.x+72] = ybdry[nbdry[0]];
	}
	__syncthreads();
	if (threadIdx.x==0)
	{
 		xybdry[blockIdx.x] = xbdry[nbdry[0]];
        	xybdry[blockIdx.x+numbdry] = ybdry[nbdry[0]];
		maxgridbdry[blockIdx.x] = nbdry[0];
//		printf("block %d thread %d\n", blockIdx.x, maxgridbdry[blockIdx.x]);
//		printf("block %d x %f y %f thread %d\n", blockIdx.x, xbdry[nbdry[0]], ybdry[nbdry[0]], nbdry[0]);
	}
}

__global__ void postsearchbdry(float* xybdry,float* leas_flux,float* xpoint,float* cu_bound_flux)
{
	__shared__ float rzbdry[numbdry*2];
	int z_axi,r_axi,num_u,num_b;
	float ra,za,rxu,zxu,rxb,zxb,agu,agb,rtempu,ztempu,agtempu,rtempb,ztempb,agtempb;

	z_axi = int(leas_flux[4]);
	r_axi = int(leas_flux[5])+68;
	
	ra = igridr+delR+delR*r_axi;	
	za = igridz+delZ+delZ*z_axi;

	if (threadIdx.x < numbdry)
	{
		rzbdry[threadIdx.x] = xybdry[threadIdx.x];
		rzbdry[threadIdx.x+numbdry] = xybdry[threadIdx.x+numbdry];
	}
	__syncthreads();

	rxu = xpoint[7];
	zxu = xpoint[8];
	agu = atan2(zxu-za,rxu-ra);
	num_u = int(agu/(2*pi/numbdry));

	rtempu = rzbdry[num_u];
	ztempu = rzbdry[num_u+numbdry];
	agtempu = atan2(ztempu-za,rtempu-ra);

	if (agu < agtempu)
	{
		num_u = num_u -1;	
	}

	rxb = xpoint[5];
	zxb = xpoint[6];
	agb = atan2(zxb-za,rxb-ra);
	num_b = numbdry + int(agb/(2*pi/numbdry)) - 1;

	rtempb = rzbdry[num_b];
	ztempb = rzbdry[num_b+numbdry];
	agtempb = atan2(ztempb-za,rtempb-ra);
	
	if (agb < agtempb)
	{
		num_b = num_b +1;	
	}
	__syncthreads();

	if (threadIdx.x <= num_u)
	{
		xybdry[threadIdx.x] = rzbdry[threadIdx.x];
		xybdry[threadIdx.x+numbdry+3] = rzbdry[threadIdx.x+numbdry];
	}
	else if (threadIdx.x == num_u+1)
	{
		if (cu_bound_flux[0] == xpoint[2])
		{
			xybdry[threadIdx.x] = rxu;
			xybdry[threadIdx.x+numbdry+3] = zxu;
//			printf("add up x-point at post-search\n");
		}
		else
		{
			xybdry[threadIdx.x] = rzbdry[threadIdx.x-1];
			xybdry[threadIdx.x+numbdry+3] = rzbdry[threadIdx.x+numbdry-1];
		}
	}
	else if (threadIdx.x <= num_b+1)
	{
		xybdry[threadIdx.x] = rzbdry[threadIdx.x-1];
		xybdry[threadIdx.x+numbdry+3] = rzbdry[threadIdx.x+numbdry-1];
	}
	else if (threadIdx.x == num_b+2)
	{
		if (cu_bound_flux[0] == xpoint[0])
		{
			xybdry[threadIdx.x] = rxb;
			xybdry[threadIdx.x+numbdry+3] = zxb;
//			printf("add down x-point at post-search\n");
		}
		else
		{
			xybdry[threadIdx.x] = rzbdry[threadIdx.x-2];
			xybdry[threadIdx.x+numbdry+3] = rzbdry[threadIdx.x+numbdry-2];
		}
	}
	else if (threadIdx.x < numbdry+2)
	{	
		xybdry[threadIdx.x] = rzbdry[threadIdx.x-2];
		xybdry[threadIdx.x+numbdry+3] = rzbdry[threadIdx.x+numbdry-2];
	}
	else
	{
		xybdry[threadIdx.x] = rzbdry[0];
		xybdry[threadIdx.x+numbdry+3] = rzbdry[numbdry];
	}
	__syncthreads();
/*
	if (threadIdx.x==0)
	{
		printf("agu %f agb %f ",agu,agb);
		printf("num_u %d num_b %d ",num_u,num_b);
	}
*/
}

__global__ void forbetapli(float* cu_result,float* leas_flux,float* profile_result,float* cu_bound_flux,float* xpoint)
{
	float invr,fluxaxi,fluxbound,br,bz,normflux,factor,interfactor;
	int a,b,c,d,e;
	__shared__ float sDiff[7][mesh];
	__shared__ float profile[mesh2];

	invr = 1/cu_radix[blockIdx.x];
	sDiff[0][threadIdx.x] = 0;
        sDiff[1][threadIdx.x] = 0;
        sDiff[2][threadIdx.x] = 0;
	sDiff[3][threadIdx.x] = 0;
        sDiff[4][threadIdx.x] = 0;
        sDiff[5][threadIdx.x] = 0;
	sDiff[6][threadIdx.x] = 0;

	profile[threadIdx.x] = profile_result[mesh2 + threadIdx.x];
//	printf("thread %d p %f \n ",threadIdx.x,profile[threadIdx.x]);
	if(threadIdx.x<2)
	{
		profile[threadIdx.x+mesh] = profile_result[mesh2 + threadIdx.x + mesh];
	}
//	__syncthreads();

	if (blockIdx.x>0 && threadIdx.x>0 && blockIdx.x<254 && threadIdx.x<254)
	{	
		sDiff[0][threadIdx.x] = cu_result[(blockIdx.x-1)*mesh1+threadIdx.x];
		sDiff[1][threadIdx.x] = cu_result[blockIdx.x*mesh1+threadIdx.x];
		sDiff[2][threadIdx.x] = cu_result[(blockIdx.x+1)*mesh1+threadIdx.x];
	}
	__syncthreads();
	
	a = (1.016 - igridr - delR)/delR;
	b = (2.3540 - igridr - delR)/delR;
	c = (xpoint[6]-igridz-delZ)/delZ;
	d = (xpoint[8]-igridz-delZ)/delZ;	
	fluxaxi = xpoint[4];
	fluxbound = cu_bound_flux[0];
	
	if(sDiff[1][threadIdx.x]>=fluxbound && threadIdx.x>(c-1) && threadIdx.x<(d+1) && blockIdx.x>a && blockIdx.x<(b+1))
	{
		factor = 1;
//		if (sDiff[0][threadIdx.x] < fluxbound && sDiff[2][threadIdx.x] > fluxbound && sDiff[1][threadIdx.x]>fluxbound)
//		{
//			factor = 0.5 + (sDiff[1][threadIdx.x]-fluxbound)/(sDiff[1][threadIdx.x]-sDiff[0][threadIdx.x]);
//			printf("%f ",factor);
//		}
//		else if (sDiff[2][threadIdx.x] < fluxbound && sDiff[0][threadIdx.x] > fluxbound && sDiff[1][threadIdx.x]>fluxbound)
//		{
//			factor = 0.5 + (sDiff[1][threadIdx.x]-fluxbound)/(sDiff[1][threadIdx.x]-sDiff[2][threadIdx.x]);
//			printf("%f ",factor);
//		}
//		else if (sDiff[0][threadIdx.x] < fluxbound && sDiff[2][threadIdx.x] < fluxbound && sDiff[1][threadIdx.x]>fluxbound)
//		{
//			factor = (sDiff[1][threadIdx.x]-fluxbound)/(sDiff[1][threadIdx.x]-sDiff[0][threadIdx.x]) + (sDiff[1][threadIdx.x]-fluxbound)/(sDiff[1][threadIdx.x]-sDiff[2][threadIdx.x]);
//			printf("%f ",factor);
//		}
//		else if (sDiff[1][threadIdx.x] == fluxbound)
//		{
//			factor = 0.5;
//			printf("%f ",factor);
//		}
		bz = invr*(sDiff[2][threadIdx.x]-sDiff[0][threadIdx.x])*INVR;
		br = -invr*(sDiff[1][threadIdx.x+1]-sDiff[1][threadIdx.x-1])*INVZ;
		normflux = (sDiff[1][threadIdx.x]-fluxaxi)/(fluxbound-fluxaxi); 
		e = int(normflux*256);
		interfactor = float(normflux*256 - e); 

//		printf("e %d interfactor %f p %f \n ",e,interfactor,(interfactor*profile[e+1]+(1-interfactor)*profile[e]));

//		sDiff[3][threadIdx.x] = factor*2*pi*float(cu_radix[blockIdx.x])*(br*br + bz*bz);  

		sDiff[3][threadIdx.x] = factor*2*pi*cu_radix[blockIdx.x]*(br*br + bz*bz); 
		sDiff[4][threadIdx.x] = factor*2*pi*(cu_radix[blockIdx.x]) * (interfactor*profile[e+1]+(1-interfactor)*profile[e]);
		sDiff[5][threadIdx.x] = factor*2*pi*cu_radix[blockIdx.x];
		sDiff[6][threadIdx.x] = factor*1;
	}
	__syncthreads();

	for(unsigned int s=128;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<255) 
		{
			sDiff[3][threadIdx.x] += sDiff[3][k];      //bp2
			sDiff[4][threadIdx.x] += sDiff[4][k];	//p
			sDiff[5][threadIdx.x] += sDiff[5][k];	//volume
			sDiff[6][threadIdx.x] += sDiff[6][k];    //cross area
		}
		__syncthreads();
	}
	
	if(threadIdx.x==0)
	{
		leas_flux[500000+blockIdx.x]=sDiff[3][0];
		leas_flux[500000+blockIdx.x+mesh2]=sDiff[4][0];
		leas_flux[500000+blockIdx.x+mesh2*2]=sDiff[5][0];
		leas_flux[500000+blockIdx.x+mesh2*3]=sDiff[6][0];
	}	
}

__global__ void postbetapli(float* leas_flux,float* xybdry,float* cu_compute_diagnotics)
{
	__shared__ float rzbdry[(numbdry+3)*2];
	__shared__ float sumbdry[numbdry+2];
	__shared__ float sum[mesh][4];
	float a,b;
	int c;

	if (threadIdx.x < numbdry+3)
	{
		rzbdry[threadIdx.x] = xybdry[threadIdx.x];
		rzbdry[threadIdx.x+numbdry+3] = xybdry[threadIdx.x+numbdry+3];
		__syncthreads();
	
//		c = numbdry/2+numbdry*int(threadIdx.x/(1+numbdry/2))-threadIdx.x;
		c = threadIdx.x;
		xybdry[c*2] = rzbdry[threadIdx.x];
		xybdry[c*2+1] = rzbdry[threadIdx.x+numbdry+3];
	}
	__syncthreads();


	sum[threadIdx.x][0] = leas_flux[500000+threadIdx.x] ;     //bp2
	sum[threadIdx.x][1] = leas_flux[500000+threadIdx.x+mesh2];		//p
	sum[threadIdx.x][2] = leas_flux[500000+threadIdx.x+mesh2*2];	//volume
	sum[threadIdx.x][3] = leas_flux[500000+threadIdx.x+mesh2*3];	//cross area

	__syncthreads();	
	
	if (threadIdx.x < (numbdry+2))
	{
		a = rzbdry[threadIdx.x]-rzbdry[threadIdx.x+1];
		b = rzbdry[threadIdx.x+numbdry+3]-rzbdry[threadIdx.x+numbdry+3+1];
		sumbdry[threadIdx.x] = sqrt(a*a+b*b);
	}
/*
	else if (threadIdx.x == (numbdry-1))
	{
		a = rzbdry[threadIdx.x]-rzbdry[0];
		b = rzbdry[threadIdx.x+numbdry]-rzbdry[numbdry];
		sumbdry[threadIdx.x] = sqrt(a*a+b*b); 
	}
*/
	__syncthreads();

	for(unsigned int s=128;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<255) 
		{
			sum[threadIdx.x][0] += sum[k][0];   //bp2
			sum[threadIdx.x][1] += sum[k][1];	//p
			sum[threadIdx.x][2] += sum[k][2];	//volume
			sum[threadIdx.x][3] += sum[k][3];     //cross area
		}
		__syncthreads();
	}
	__syncthreads();

	for(unsigned int s=64;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<numbdry+2)
                {
                        sumbdry[threadIdx.x] += sumbdry[k];  //l
                }
		__syncthreads();
	}
//	sumbdry[0] = 5.364000;
	__syncthreads();
	
	if (threadIdx.x==0)
	{
		leas_flux[590000] = sum[0][2]*darea;   //volume
		leas_flux[590001] = tmu * sum[0][1] * sumbdry[0] * sumbdry[0] /(sum[0][2] * cu_compute_diagnotics[120]*cu_compute_diagnotics[120] * 2*pi);       //betap
		leas_flux[590002] = sum[0][3]*darea;   //cross area
		leas_flux[590003] = 10000*10000 * (sum[0][0]/sum[0][2]) * sumbdry[0]*sumbdry[0] / (16 * pi*pi * (cu_compute_diagnotics[120]/1000)*(cu_compute_diagnotics[120]/1000));    //li
		leas_flux[590004] = sumbdry[0];        //boundary
		leas_flux[590005] = sum[0][1]*1.5*darea;     //energy

		printf("volume %f betap %f cross_area %f li %f sumbdry %f energy %f ",leas_flux[590000],leas_flux[590001],leas_flux[590002],leas_flux[590003],sumbdry[0],leas_flux[590005]);

//		printf("%f %f %f ",leas_flux[590001],leas_flux[590003],leas_flux[590005]);
	}
/*
	if (threadIdx.x==0)
	{
		leas_flux[900+threadIdx.x] = tmu * sum[0][1] * sumbdry[0] * sumbdry[0] /(sum[0][2]*darea * ip*ip * 2*pi);	//betap
	}
	
	if (threadIdx.x==1)
	{
		leas_flux[900+threadIdx.x] = sum[0][2]*darea;	//cross area
	}

	if (threadIdx.x==2)
	{
		leas_flux[900+threadIdx.x] =  tmu * sum[0][0] * sumbdry[0]/(sum[0][2] * darea * ip);	//li
	}
*/
}

__global__ void prime_profile_compute(float* test_right, float* result, float* current_representation)
{
	float normflux,x[20],sum,y[21];
	int kppcur,kffcur,pbdry,ffbdry,num,i;

	kppcur = (int)current_representation[1];
	kffcur = (int)current_representation[2];
	pbdry = (int)current_representation[3]; 
	ffbdry = (int)current_representation[4];
	num = (int)current_representation[25];

	normflux = float(threadIdx.x)/mesh1;
//	normfactor = xpoint[4] - boundflux[0]; 
	sum = 0;
	y[0] = 1;

	for (i=0; i<num; i++)
	{
		x[i] = float(test_right[EF_NFCOIL + i]);
		y[i+1] = y[i]*normflux; 		
	}

	if (blockIdx.x == 0)   // ffprime
	{
		for (i=0; i<kffcur; i++)
		{
			sum = sum + x[kppcur+i]*(y[i]-(ffbdry)*y[kffcur]);
		}
//		result[blockIdx.x*mesh2 + threadIdx.x] = 4*pi*(ff1 - ff1*normflux)/(darea*10000000);
		result[(blockIdx.x+2)*mesh2 + threadIdx.x] = 4*pi*sum/(darea*10000000);

//		printf("FF' thread %d flux %f \n",threadIdx.x,result[(blockIdx.x+2)*mesh2 + threadIdx.x]);
	}
	else   // pprime
	{
		for (i=0; i<kppcur; i++)
		{
			sum = sum + x[i]*(y[i]-(pbdry)*y[kppcur]);
		}
//		result[blockIdx.x*mesh2 + threadIdx.x] = (p1 + p2*normflux - (p1+p2)*normflux*normflux)/darea;
		result[(blockIdx.x+2)*mesh2 + threadIdx.x] = sum/darea;	

//		printf("P' thread %d flux %f \n",threadIdx.x,result[(blockIdx.x+2)*mesh2 + threadIdx.x]);
	}
}

__global__ void prime_profile_compute_spline(float* test_right, float* result, float* current_representation, float* y)
{
	float normflux,sum;
	int num_knots,i,index;

	num_knots = (int)current_representation[1];

	normflux = float(threadIdx.x)/mesh1;

	if (blockIdx.x == 0)   // ffprime
	{
		for (i=0; i<num_knots-1; i++)
		{
			if (normflux >= y[i] && normflux <= y[i+1])
			{
				index = EF_NFCOIL + i*8 + 4;
				sum = test_right[index] + test_right[index+1]*normflux + test_right[index+2]*sinh(tension*normflux) + test_right[index+3]*cosh(tension*normflux);
			}
		}
//		result[blockIdx.x*mesh2 + threadIdx.x] = 4*pi*(ff1 - ff1*normflux)/(darea*10000000);
		result[(blockIdx.x+2)*mesh2 + threadIdx.x] = 4*pi*sum/(darea*10000000);

//		printf("FF' thread %d flux %f \n",threadIdx.x,result[(blockIdx.x+2)*mesh2 + threadIdx.x]);
	}
	else   // pprime
	{
		for (i=0; i<num_knots-1; i++)
		{
			if (normflux >= y[i] && normflux<=y[i+1])
			{
				index = EF_NFCOIL + i*8;
				sum = test_right[index] + test_right[index+1]*normflux + test_right[index+2]*sinh(tension*normflux) + test_right[index+3]*cosh(tension*normflux);
			}
		}
//		result[blockIdx.x*mesh2 + threadIdx.x] = (p1 + p2*normflux - (p1+p2)*normflux*normflux)/darea;
		result[(blockIdx.x+2)*mesh2 + threadIdx.x] = sum/darea;	

//		printf("P' thread %d flux %f \n",threadIdx.x,result[(blockIdx.x+2)*mesh2 + threadIdx.x]);
	}
}

__global__ void profile_compute(float* result, float* xpoint, float* boundflux, float* btmr, float* current_representation)
{
	float normfactor,btr,del;
	int i;
	__shared__ float prime[mesh2];
	__shared__ float profile[mesh2];

//	kppcur = (int)current_representation[1];
//	kffcur = (int)current_representation[2];
//	pbdry = (int)current_representation[3]; 
//	ffbdry = (int)current_representation[4];
//	num = (int)current_representation[25];

	btr = 1.69550002*btmr[0];
	del = float(1)/float(mesh1);
	normfactor = xpoint[4] - boundflux[0]; 

	prime[threadIdx.x] = result[(blockIdx.x+2)*mesh2 + threadIdx.x];
	__syncthreads();

	if (blockIdx.x == 0)   // fpol
	{
//		prime[threadIdx.x] = result[(blockIdx.x+2)*mesh2 + threadIdx.x];
//		__syncthreads();
		if (threadIdx.x == 0)
		{
			profile[mesh1] = 0;
			for (i=1;i<mesh2;i++)
			{
				profile[mesh1-i] = profile[mesh1-i+1] + del*(prime[mesh1-i]+prime[mesh1-i+1])/2;
//				printf("mesh1-i %d ff %f ff' %f del %f \n ",mesh1-i,profile[mesh1-i],prime[mesh1-i],del);
			}
		}
		__syncthreads();
		result[blockIdx.x*mesh2 + threadIdx.x] = sqrt(btr*btr + 2*normfactor*profile[threadIdx.x]);
//		printf("thread %d FF' %f F %f \n",threadIdx.x,prime[threadIdx.x],profile[threadIdx.x]);
//		printf("F thread %d flux %f \n",threadIdx.x,result[blockIdx.x*mesh2 + threadIdx.x]);
	}
	else if (blockIdx.x == 1)  // p
	{
//		prime[threadIdx.x] = result[(blockIdx.x+2)*mesh2 + threadIdx.x];
//		__syncthreads();
		if (threadIdx.x == 0)
		{
			profile[mesh1] = 0;
			for (i=1;i<mesh2;i++)
			{
				profile[mesh1-i] = profile[mesh1-i+1] + del*(prime[mesh1-i]+prime[mesh1-i+1])/2;
//				printf("mesh1-i %d p %f p' %f \n ",mesh1-i,profile[mesh1-i],prime[mesh1-i]);
			}
		}
		__syncthreads();
		result[blockIdx.x*mesh2 + threadIdx.x] = normfactor*profile[threadIdx.x]; 
//		printf("thread %d P' %f P %f \n",threadIdx.x,prime[threadIdx.x],profile[threadIdx.x]);
//		printf("P thread %d flux %f \n",threadIdx.x,result[blockIdx.x*mesh2 + threadIdx.x]);
//		printf("P' thread %d flux %f \n",threadIdx.x,result[(blockIdx.x+2)*mesh2 + threadIdx.x]);
	}

}

/*
__global__ void profile_compute(float* test_right, float* result, float* xpoint, float* boundflux, float* btmr, float* current_representation)
{
	float normflux,normfactor,btr,x[20],sum,y[21];
	int kppcur,kffcur,pbdry,ffbdry,fitdz,num,i;

	kppcur = (int)current_representation[1];
	kffcur = (int)current_representation[2];
	pbdry = (int)current_representation[3]; 
	ffbdry = (int)current_representation[4];
	num = (int)current_representation[25];

	normflux = float(threadIdx.x)/mesh1;
	btr = 1.69550002*btmr[0];
	normfactor = xpoint[4] - boundflux[0]; 
	sum = 0;
	y[0] = 1;

	for (i=0; i<num; i++)
	{
		x[i] = float(test_right[EF_NFCOIL + i]);
		y[i+1] = y[i]*normflux; 		
	}

	if (blockIdx.x == 0)   // ff
	{
		for (i=0; i<kffcur; i++)
		{
			sum = sum + x[i]*((1-y[i+1])/(i+1) - (ffbdry)*(1-y[kffcur+1])/(kffcur+1));                       //  (y[i]-(1-pbdry)*y[kppcur]);
		}
		result[blockIdx.x*mesh2 + threadIdx.x] = sqrt(btr*btr + 4*pi*normfactor*sum/(darea*10000000));
		//result[blockIdx.x*mesh2 + threadIdx.x] = sqrt(btr*btr+4*pi*2*normfactor*((ff1 - 0.5*ff1) - (ff1*normflux - 0.5*ff1*normflux*normflux))/(darea*10000000)); 
	}
	else if (blockIdx.x == 1)  // p
	{
		for (i=0; i<kppcur; i++)
		{
			sum = sum + x[i]*((1-y[i+1])/(i+1) - (pbdry)*(1-y[kppcur+1])/(kppcur+1));                       //  (y[i]-(1-pbdry)*y[kppcur]);
		}
		result[blockIdx.x*mesh2 + threadIdx.x] = normfactor*sum/darea; 
		//result[blockIdx.x*mesh2 + threadIdx.x] = normfactor*((p1 + 0.5*p2 -(p1+p2)/3) - (p1*normflux + 0.5*p2*normflux*normflux - normflux*normflux*normflux*(p1+p2)/3))/darea;
	}
	else if (blockIdx.x == 2)   // ffprime
	{
		for (i=0; i<kffcur; i++)
		{
			sum = sum + x[kppcur+i]*(y[i]-(ffbdry)*y[kffcur]);
		}
//		result[blockIdx.x*mesh2 + threadIdx.x] = 4*pi*(ff1 - ff1*normflux)/(darea*10000000);
		result[blockIdx.x*mesh2 + threadIdx.x] = 4*pi*sum/(darea*10000000);
	}
	else   // pprime
	{
		for (i=0; i<kppcur; i++)
		{
			sum = sum + x[i]*(y[i]-(pbdry)*y[kppcur]);
		}
//		result[blockIdx.x*mesh2 + threadIdx.x] = (p1 + p2*normflux - (p1+p2)*normflux*normflux)/darea;
		result[blockIdx.x*mesh2 + threadIdx.x] = sum/darea;
		
//		printf("thread %d flux %f \n",threadIdx.x,y[1]);
	}
}
*/

__global__ void search_surface(float* cu_result,float* location,float* tani,float* coti,float* xpoint,int* num,float* xy_surface)
{
//	__shared__ int N[32];
	int ra,rb,za,zb,r_axi,z_axi,i,N1,N2;
	float si1,si2,flux,dflux,aflux,factor;

	z_axi = int(location[4]);
	r_axi = int(location[5])+68;
	dflux = xpoint[11]/mesh1;
	aflux = xpoint[4];

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

			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);
		
			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);
				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
				xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;
//				n_surface[i] = i;
				
				__syncthreads();
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
		
			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);

			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);

				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
				xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;

				__syncthreads();
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
			
			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);

			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);

				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        	xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*delZ;
				
				__syncthreads();
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
	
			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);

			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);

				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        	xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;

				__syncthreads();
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
                
			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);

			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);

				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra-factor*delR;
                        	xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ;

				__syncthreads();
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
		
			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);

			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);

				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        	xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;

				__syncthreads();
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
			
			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);

			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);

				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+(rb-ra)*factor*delR;
                        	xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za-factor*delZ;

				__syncthreads();
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

			N1 = int((aflux-si1)/dflux);
			N2 = int((aflux-si2)/dflux);

			for(i=N1+1;i<N2+1;i++)
			{
				flux = aflux - i*dflux;
				factor = (flux-si1)/(si2-si1);

				xy_surface[i*2*numbdry+blockIdx.x] = igridr+delR+delR*ra+factor*delR;
				xy_surface[(i*2+1)*numbdry+blockIdx.x] = igridz+delZ+delZ*za+factor*(zb-za)*delZ; 
				
				__syncthreads();
			}
			__syncthreads();
		}
		__syncthreads();
	}
}

__global__ void q_profile(float* cu_result,float* xpoint,float* xy_surface,float* btor,float* cu_qpsi)
{
	int ra,za;
	float r,z,factor1,factor2,bp,br,bz,bt,dl,a,b;
	float flux1,flux2,flux3,flux4,flux_a,flux_b,flux_c,flux_d;
	__shared__ float sum[numbdry];

	r = xy_surface[(blockIdx.x+1)*2*numbdry+threadIdx.x];
	z = xy_surface[((blockIdx.x+1)*2+1)*numbdry+threadIdx.x];

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
	bt = fabs(1.69550002*btor[0]/r);	

	if (threadIdx.x < (numbdry-1))
	{
		a = xy_surface[(blockIdx.x+1)*2*numbdry+threadIdx.x+1];
		b = xy_surface[((blockIdx.x+1)*2+1)*numbdry+threadIdx.x+1];
		dl = sqrt((r-a)*(r-a)+(z-b)*(z-b));
	}
	else if (threadIdx.x == (numbdry-1))
	{
		a = xy_surface[(blockIdx.x+1)*2*numbdry];
		b = xy_surface[((blockIdx.x+1)*2+1)*numbdry];
		dl = sqrt((r-a)*(r-a)+(z-b)*(z-b)); 
	}
	sum[threadIdx.x] = bt*dl/(2*pi*r*bp);
	__syncthreads();

	for(unsigned int s=64;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<numbdry)
                {
                        sum[threadIdx.x] += sum[k];
                }
		__syncthreads();
	}
	if (threadIdx.x == 1)
	{
		cu_qpsi[blockIdx.x+1] = sum[0];
		__syncthreads();
		cu_qpsi[0] = cu_qpsi[1];
	}
}

__global__ void diagnostic_compute(float* cu_gfc,float* test_right,float* result, float* current_representation)
{
	int num;
	__shared__ float temp_left[64];
	__shared__ float temp_right[64];

	unsigned int index = threadIdx.x * EF_DIAG + blockIdx.x;
	num = EF_NFCOIL+(int)current_representation[25];
	temp_left[threadIdx.x] = cu_gfc[index];
	temp_right[threadIdx.x] = test_right[threadIdx.x]; 

	temp_left[threadIdx.x] = temp_left[threadIdx.x] * temp_right[threadIdx.x];
	
	for(unsigned int s=32;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<num)
                {
                        temp_left[threadIdx.x] += temp_left[k];
                }
		__syncthreads();
	}
	result[blockIdx.x] = temp_left[0];

//	if (threadIdx.x == 1)
//	//	{
//	//		printf("%d %f\n",blockIdx.x, result[blockIdx.x]);
//	//	}
}

__global__ void chi_compute(float* result,float* diagnostic,float* fwt,float* leas_flux)
{
	__shared__ float chi[EF_DIAG];

	float w,m,k;

	m = diagnostic[threadIdx.x];
	w = fwt[threadIdx.x];
	k = result[threadIdx.x];

	__syncthreads();
	
	if (w != 0)
	{
		chi[threadIdx.x] = (m/w - k)*(m/w - k)*(w*w);
	}
	else
	{
		chi[threadIdx.x] = 0;
	}

//	printf("%d %f %f %f\n",threadIdx.x, chi[threadIdx.x],k,w);
	
	__syncthreads();	
	for(unsigned int s=64;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<EF_DIAG)
               	{
			chi[threadIdx.x] += chi[k];
		}
		__syncthreads();
	}
	if (threadIdx.x == 0)
	{	
		printf("chi = %f ", chi[threadIdx.x]);
		leas_flux[590006] = chi[threadIdx.x];	
	}
}	
