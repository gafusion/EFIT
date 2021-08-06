void tri_solver_init_fb()
{
	int i,j;
    	double tri_down[mesh-1];
   	double tri_up[mesh-1];
    	double tri_main[mesh*mesh];
    	double act_trace[mesh*mesh];
	double currentInver[mesh];
	for (i=0;i<(mesh-1);i++)
    	{
        	tri_down[i]=-1-delR/(2*(igridr+delR*(i+2)));
        	tri_up[i]=-1+delR/(2*(igridr+delR*(i+1)));
//		currentInver[i]=(igridr+delR*i)*delR*delR;
		currentInver[i]=(igridr+delR*(i+1))*delR*delR*0.0000001*4*pi;		
    	}
//    	tri_up[mesh-1]=0;
        currentInver[mesh-1]=(igridr+delR*mesh)*delR*delR*0.0000001*4*pi;
    	double c=RTOZ;
    	double cos_temp;
    	double constant_main;
    	for(i=0;i<mesh;i++)
    	{
        	cos_temp=(i+1)*pi/(mesh+1.0);
        	constant_main=2+c*(2-2*cos(cos_temp));   
        	*(tri_main+mesh*i)=constant_main;
        	for(j=1;j<mesh;j++)
        	{
            		*(act_trace+mesh*i+j)=-tri_down[j-1]/(*(tri_main+mesh*i+j-1));
            		*(tri_main+mesh*i+j)= constant_main+*(act_trace+mesh*i+j)*tri_up[j-1];    
        	}	
    	}

    	cudaMalloc((void**)&cu_act_trace_fb,((mesh*mesh)*sizeof(double)));
    	cudaMemcpy(cu_act_trace_fb,act_trace,(mesh*mesh)*sizeof(double),cudaMemcpyHostToDevice);
    	cudaMalloc((void**)&cu_tri_main_fb,((mesh*mesh)*sizeof(double)));
    	cudaMemcpy(cu_tri_main_fb,tri_main,(mesh*mesh)*sizeof(double),cudaMemcpyHostToDevice);
   	double ele_bfft[(mesh+1)*(mesh+1)];
    	double ele[mesh*mesh];
    	double fft_mat[(mesh+1)*(mesh+1)];
    	for(i=0;i<=(mesh*mesh-1);i++)
    	{
       		*(ele+i)=i+1;
    	}
    	for(i=0;i<=(mesh*mesh+2*mesh);i++)
    	{
        	*(ele_bfft+i)=0;
        	fft_mat[i]=0;
    	}
    	for(i=0;i<mesh;i++)
        	for(j=0;j<=mesh-1;j++)
        	{
            		fft_mat[i*(mesh+1)+j]=sin((i+1)*(j+1)*pi/(mesh+1.0))*sqrt(2.0/(mesh+1.0));
           	 	*(ele_bfft+i*(mesh+1)+j)=ele[i*mesh+j];
        	}
	double radix[mesh];
	for(i=0;i<mesh;i++)
	{
		radix[i]=igridr+delR+delR*i;
	}
    	cudaMemcpyToSymbol(cu_radix_fb,radix,(mesh)*sizeof(double));
    	cudaMemcpyToSymbol(cu_currentInver_fb,currentInver,mesh*sizeof(double));
    	cudaMalloc((void**)&cu_mat_fb, (mesh+1)*(mesh+1)*sizeof(double));
    	cudaMemcpy(cu_mat_fb,fft_mat,(mesh+1)*(mesh+1)*sizeof(double) , cudaMemcpyHostToDevice);

	double tani[numbdry];
	double coti[numbdry];
	for (i=0;i<numbdry;i++)
	{
		tani[i] = tan(i*2*pi/numbdry);
		coti[i] = 1/tani[i];
//		printf("% 10.9E % 10.9E\n",tani[i],coti[i]);
	}
	cudaMalloc((void**)&cu_tani_fb, numbdry*sizeof(double));
        cudaMemcpy(cu_tani_fb,tani,numbdry*sizeof(double) , cudaMemcpyHostToDevice);
	cudaMalloc((void**)&cu_coti_fb, numbdry*sizeof(double));
        cudaMemcpy(cu_coti_fb,coti,numbdry*sizeof(double) , cudaMemcpyHostToDevice);
}

void boundData_fb(double* a, double* b)
{
	for(int i=0; i<ELI_NUM*4*(mesh+1);i++)
		b[i] = 0;
/*
	for(int i=0;i<63;i++)
	{
		for(int j=0;j<ELI_NUM;j++)
		{
			b[63*4*i+j] = a[ELI_NUM*(i+1) + j];
			b[63*4*(63+i*2)+j] = a[(i+1)*65*ELI_NUM+j];
			b[63*4*(63+i*2+1)+j] = a[((i+1)*65+64)*ELI_NUM+j];
			b[63*4*(63*3+i)+j] = a[(65*64+1+i)*ELI_NUM+j];						
		}

	}
*/
	for(int i=0;i<ELI_NUM;i++)
	{
		for(int j=0;j<(mesh*4);j++)
		{
			b[(mesh+1)*4*i+j] = a[mesh*4*i+j]; 
//			b[(mesh+1)*4*i+mesh+2*j] = a[nwnh*i+(j+1)*(mesh+2)];
//			b[(mesh+1)*4*i+mesh+2*j+1] = a[nwnh*i+(j+1)*(mesh+2)+(mesh+1)];;
//			b[(mesh+1)*4*i+mesh*3+j] =a[nwnh*i+(j+1)+(mesh+2)*(mesh+1)];
		}

	}
}

int read_gfunction_files_fb()
{
	int length3 = mesh*4*ELI_NUM;
	double* hgridresponse;
	hgridresponse= (double *)malloc(length3*sizeof(double));
 	FILE* Number4;
	Number4=fopen("/u/huangyao/efund_postprocess257_double/gridpc/ef_gridresponse.dat","r");
	if(Number4==0)
	{
		printf("error: can not open ef_gridresponse.dat\n");
		return 1;
	}
	if(!fread(hgridresponse,sizeof(double),length3,Number4))
	{
		printf("File 4 is empty\n");
		return 1;
	}	
	fclose(Number4);

	double* boundresponse;
	boundresponse = (double *)malloc(ELI_NUM*(mesh+1)*4*sizeof(double));

	boundData_fb(hgridresponse,boundresponse);

    	cudaMalloc((void**)&cu_boundresponse_fb,(ELI_NUM*(mesh+1)*4*sizeof(double)));
    	cudaMemcpy(cu_boundresponse_fb,boundresponse,ELI_NUM*(mesh+1)*4*sizeof(double),cudaMemcpyHostToDevice);
	free(hgridresponse);
	free(boundresponse);

	double* hicresponse;
	hicresponse = (double *)malloc(EF_NFCOIL*nwnh*sizeof(double));
	FILE* PFCR;
	PFCR=fopen("/u/huangyao/efund_postprocess257_double/gridfc&ec/ef_pfresponse.dat","r");
	if(PFCR==0)
	{
		printf("error: can not open ef_pfresponse.dat\n");
		return 1;
	}
	if(!fread(hicresponse,sizeof(double),EF_NFCOIL*nwnh,PFCR))
	{
		printf("File 4 is empty\n");
		return 1;
	}	
	fclose(PFCR);
	cudaMalloc((void**)&cu_icresponse_fb,(EF_NFCOIL*nwnh*sizeof(double)));
    	cudaMemcpy(cu_icresponse_fb,hicresponse,EF_NFCOIL*nwnh*sizeof(double),cudaMemcpyHostToDevice);
	free(hicresponse);

	double* ecgridresponse;
	ecgridresponse= (double *)malloc(EF_ECCOIL*nwnh*sizeof(double));
 	FILE* ECCR;
	ECCR=fopen("/u/huangyao/efund_postprocess257_double/gridfc&ec/ef_ecresponse.dat","r");
	if(ECCR==0)
	{
		printf("error: can not open ef_ecresponse.dat\n");
		return 1;
	}
	if(!fread(ecgridresponse,sizeof(double),EF_ECCOIL*nwnh,ECCR))
	{
		printf("File 8 is empty\n");
		return 1;
	}	
	fclose(ECCR);
	cudaMalloc((void**)&cu_ecgridresponse_fb,(nwnh*EF_ECCOIL*sizeof(double)));
    	cudaMemcpy(cu_ecgridresponse_fb,ecgridresponse,nwnh*EF_ECCOIL*sizeof(double),cudaMemcpyHostToDevice);
	free(ecgridresponse);
	
	return 0;
}

__global__ void data_eli_init_fb(double* cu_result, double* effect_flux, int* cu_num_se, double* ffprime, double* pprime, double* r_bdry, double* z_bdry, int* cu_numbdry_gafile, double* cu_psiaxis_gafile, double *cu_psibdry_gafile)
{
	int offset,a,b,c,d,e,f,g;
	double fluxaxis,fluxbound,normfactor,interfactor,x;

	__shared__ double ffprime_buffer[mesh2];
	__shared__ double pprime_buffer[mesh2];

	__shared__ double buffer[mesh][2];
	__shared__ double normedflux[mesh];

	fluxaxis = cu_psiaxis_gafile[0];
	fluxbound = cu_psibdry_gafile[0];

	ffprime_buffer[threadIdx.x] = ffprime[threadIdx.x];
	pprime_buffer[threadIdx.x] = pprime[threadIdx.x];

	if (threadIdx.x < 2)
	{
		ffprime_buffer[threadIdx.x+mesh] = ffprime[threadIdx.x+mesh];
		pprime_buffer[threadIdx.x+mesh] = pprime[threadIdx.x+mesh];
	}

	unsigned int index = (blockIdx.x+1)*mesh2+threadIdx.x+1;
	buffer[threadIdx.x][1] = cu_result[index];
	normfactor=1/(fluxbound-fluxaxis);

	c = (z_bdry[cu_numbdry_gafile[0]]-igridz-delZ)/delZ;
	d = (z_bdry[cu_numbdry_gafile[0]+2]-igridz-delZ)/delZ;
	e = (r_bdry[cu_numbdry_gafile[0]+1]-igridr-delR)/delR;
	f = (r_bdry[cu_numbdry_gafile[0]+3]-igridr-delR)/delR;

	__syncthreads();

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
		normedflux[b]=(buffer[b][0]-fluxaxis)*normfactor;
		if(normedflux[b]>=0&&normedflux[b]<=1&&a>(c-1)&&a<(d+1)&&blockIdx.x>(e-1)&&blockIdx.x<(f+1))
		{
			g = int(normedflux[b]*mesh1);
			interfactor = normedflux[b]*mesh1 - g;

			x = cu_radix_fb[blockIdx.x]*((1-interfactor)*pprime_buffer[g] + interfactor*pprime_buffer[g+1]) + ((1-interfactor)*ffprime_buffer[g] + interfactor*ffprime_buffer[g+1])/cu_radix_fb[blockIdx.x];
//			printf("blockIdx.x %d threadIdx.x %d interfactor %f x %f \n",blockIdx.x, threadIdx.x, interfactor, x);
		}
		else
		{
			x = 0;
		}
		//begin to store the transposed data into global memory
		
		effect_flux[offset+b] = x;	
//		printf("data eli %d = %f \n",offset+b,x);	
	}
}

__global__ void boundflux_fb(double* cu_boundresponse, double* effect_flux,double* leas_flux)
{
	double current;
	double result=0;
	int index;
	for(int i=0;i<32;i++)
	{
		current = effect_flux[blockIdx.x*32+i];
		index = (blockIdx.x*32+i)*1024+threadIdx.x;
		result += cu_boundresponse[index]*current;
	}
	leas_flux[blockIdx.x*1024+threadIdx.x] = result;
}

__global__ void post_boundflux_fb(double* leas_flux,double* cuBoundResult)
{
	double result=0;
	int index;
	if(blockIdx.x<34)
	{
		for(int i=0;i<32;i++)
		{
			index = (blockIdx.x*32+i)*1024+threadIdx.x;
			result += leas_flux[index];
		}
		cuBoundResult[blockIdx.x*1024+threadIdx.x] = result;	
	}
	else
	{
		for(int i=0;i<28;i++)
		{
			index = (blockIdx.x*32+i)*1024+threadIdx.x;
			result += leas_flux[index];
		}
		cuBoundResult[blockIdx.x*1024+threadIdx.x] = result;	
	}
}

__global__ void final_boundflux_fb(double* leas_flux,double* cuBoundResult)
{
	double result=0;
	int index;
	for(int i=0;i<35;i++)
	{
		index = i*1024+threadIdx.x;
		result += leas_flux[index];
	}
	cuBoundResult[blockIdx.x*1024+threadIdx.x] = result;
}

__global__ void data_regroup_fb(double* effect_flux,double* leas_flux,int* cu_num_se)
{
	int offset,a;
	__shared__ double data_buffer[mesh1];
	data_buffer[threadIdx.x] = 0;
	__syncthreads();
	if(blockIdx.x<mesh)
	{
		if(threadIdx.x<cu_num_block[blockIdx.x])
		{
			offset=cu_blo_accu[blockIdx.x];
			a=cu_num_se[offset+threadIdx.x];
			data_buffer[a] =effect_flux[offset + threadIdx.x];	
		}
	}
	__syncthreads();
	leas_flux[blockIdx.x*mesh1+threadIdx.x] = data_buffer[threadIdx.x];
}

__global__ void gs_solver_init_fb(double* cuBoundResult, double* cu_bfft)
{
	__shared__ double sbfft[mesh];
	sbfft[threadIdx.x] = cu_bfft[blockIdx.x*mesh1+threadIdx.x]/(delR*delZ);
	__syncthreads();
	if(blockIdx.x == 0)
	{
		sbfft[threadIdx.x] = (1+delR/(2*(igridr+delR)))*cuBoundResult[threadIdx.x];		
	}
	else if(blockIdx.x == mesh-1)
		sbfft[threadIdx.x] = (1-delR/(2*(egridr-delR)))*cuBoundResult[mesh*3+threadIdx.x];
	else
		sbfft[threadIdx.x] *= cu_currentInver_fb[blockIdx.x];
	if(threadIdx.x == 0)
	{
		sbfft[threadIdx.x] += RTOZ*cuBoundResult[mesh+blockIdx.x*2];
	}
	if(threadIdx.x == mesh-1)
	{
		sbfft[threadIdx.x] += RTOZ*cuBoundResult[mesh+blockIdx.x*2+1];			
	}
	__syncthreads();
	cu_bfft[blockIdx.x*(mesh+1) + threadIdx.x] = sbfft[threadIdx.x];
	
}

__global__ void fft_inver_fb(double* A, double* B, double* C)
{

   	// Block index
   	double Csub = 0;
        int c,k;

   	for (int a=(mesh+1)*BLOCK_SIZE*blockIdx.y, b= BLOCK_SIZE*blockIdx.x; a<=(mesh+1)*BLOCK_SIZE*blockIdx.y+mesh;a+=BLOCK_SIZE,b+=BLOCK_SIZE*(mesh+1)) 
	{
         	__shared__ double As[BLOCK_SIZE][BLOCK_SIZE+1];
		__shared__ double Bs[BLOCK_SIZE][BLOCK_SIZE+1];

         	As[threadIdx.y][threadIdx.x] = A[a + (mesh+1) * threadIdx.y + threadIdx.x];
         	Bs[threadIdx.y][threadIdx.x] = B[b + (mesh+1) * threadIdx.y + threadIdx.x];

         	__syncthreads();
         	for (k = 0; k < BLOCK_SIZE; ++k)
            		Csub += As[threadIdx.x][k] * Bs[k][threadIdx.y];

         	__syncthreads();
   	}
        c = (mesh+1) * BLOCK_SIZE*blockIdx.x + BLOCK_SIZE*blockIdx.y;    
   	C[c + (mesh+1) * threadIdx.y + threadIdx.x] = Csub;
}

__global__ void tri_solver_fb(double* cu_result, double* cu_act_trace, double* cu_tri_main,double* cu_bfft)
{
    	__shared__ double sharedata[mesh][3];
//    	__shared__ double act_trace[mesh];
//    	__shared__ double tri_main[mesh];
//    	__shared__ double result[mesh+1];
//    	__shared__ double tri_up[mesh];
    	unsigned int index=blockIdx.x*(mesh+1)+threadIdx.x;
    	sharedata[threadIdx.x][0]=cu_result[index];
    	unsigned int index2=blockIdx.x*mesh+threadIdx.x;
    	sharedata[threadIdx.x][1]=cu_act_trace[index2];
    	sharedata[threadIdx.x][2]=cu_tri_main[index2];
//    	tri_up[threadIdx.x]=cu_tri_up[threadIdx.x];
//    	sharedata[mesh][0]=0;    
    	__syncthreads();
    	
	int i,k;
    	double tempa,tempb;
//	if(blockIdx.x<100)

    	for(i=1;i<=p_pm;i++)
    	{
        	k=1<<(i-1);
        	if(threadIdx.x>=k)
        	{
            		tempa=sharedata[threadIdx.x][0]+sharedata[threadIdx.x-k][0]*sharedata[threadIdx.x][1];
            		tempb=sharedata[threadIdx.x][1]*sharedata[threadIdx.x-k][1];            
        	}
        	__syncthreads();
        	if(threadIdx.x>=k)
        	{
//            		tempa=sharedata[threadIdx.x][0]+sharedata[threadIdx.x-k][0]*sharedata[threadIdx.x][1];
            		sharedata[threadIdx.x][0]=tempa;
            		sharedata[threadIdx.x][1]=tempb;
                }
        	__syncthreads();      
    	}
        
   	sharedata[threadIdx.x][0]=sharedata[threadIdx.x][0]/sharedata[threadIdx.x][2];
//	__syncthreads();
    	sharedata[threadIdx.x][2]=-(-1+delR/(2*(igridr+delR*(threadIdx.x+1))))/sharedata[threadIdx.x][2];
//    	sharedata[threadIdx.x][2]=-cu_tri_up[threadIdx.x]/sharedata[threadIdx.x][2];

    	__syncthreads();

    	for(i=1;i<=p_pm;i++)
    	{
        	k=2<<(i-1)>>1;
         	if(threadIdx.x<=mesh-1-k)
        	{
            		tempa=sharedata[threadIdx.x][0]+sharedata[threadIdx.x+k][0]*sharedata[threadIdx.x][2];
            		tempb=sharedata[threadIdx.x][2]*sharedata[threadIdx.x+k][2];            
        	}
		__syncthreads();
         	if(threadIdx.x<=mesh-1-k)
        	{
            		sharedata[threadIdx.x][0]=tempa;
            		sharedata[threadIdx.x][2]=tempb;            
        	}
		__syncthreads();

    	}

    	cu_bfft[index]=sharedata[threadIdx.x][0];    
}

__global__ void	ecFlux_fb(double* test_diagnostic2, double* cu_ecgridresponse, double* cu_ecflux)
{
	__shared__ double sresponse[129][6];
	unsigned int index = blockIdx.x * 774 + threadIdx.y*6 + threadIdx.x;
	sresponse[threadIdx.y][threadIdx.x] = cu_ecgridresponse[index];
	__syncthreads();
	sresponse[threadIdx.y][threadIdx.x] *= test_diagnostic2[threadIdx.x]; 
	__syncthreads();
	double temp_sum = 0;
	if(threadIdx.x == 0)
	{
		temp_sum = sresponse[threadIdx.y][0]+sresponse[threadIdx.y][1]+sresponse[threadIdx.y][2]+sresponse[threadIdx.y][3]+sresponse[threadIdx.y][4]+sresponse[threadIdx.y][5];	
		cu_ecflux[blockIdx.x*129+threadIdx.y]=temp_sum;
	}
}

__global__ void data_process_fb_right(double* cu_result, double* ec_flux, double* boundflux, double* r_bdry, double* z_bdry, double* flux_bdry_right,double* fwtbdry)
{
	int rID,zID,rIDe,zIDe;
        double flux,flux1,flux2,flux3,flux4,m1,m2,m3,m4,r,z,r0,z0; 

	__shared__ double flux_buffer[1];

	r = r_bdry[threadIdx.x];
	z = z_bdry[threadIdx.x];

	rID = int((r-igridr-delR)/delR);
	zID = int((z-igridz-delZ)/delZ);

	rIDe = int((r-igridr)/delR);
	zIDe = int((z-igridz)/delZ);

	r0 = (rID+1)*delR+igridr;
	z0 = (zID+1)*delZ+igridz;

	flux1 = cu_result[rID*mesh1+zID] + ec_flux[rIDe*mesh2+zIDe];
	flux2 = cu_result[(rID+1)*mesh1+zID] + ec_flux[(rIDe+1)*mesh2+zIDe];
	flux3 = cu_result[rID*mesh1+zID+1] + ec_flux[rIDe*mesh2+zIDe+1];
	flux4 = cu_result[(rID+1)*mesh1+zID+1] + ec_flux[(rIDe+1)*mesh2+zIDe+1];

	m1 = (r0+delR-r)*(z0+delZ-z);
	m2 = (r-r0)*(z0+delZ-z);
	m3 = (r0+delR-r)*(z-z0);
	m4 = (r-r0)*(z-z0);

	flux = (m1*flux1+m2*flux2+m3*flux3+m4*flux4)/darea;

	if (threadIdx.x == 0)
	{
		flux_buffer[0] = flux;
	}
	
	__syncthreads();

	if (threadIdx.x == 0)
	{
		flux_bdry_right[threadIdx.x] = (boundflux[0] - flux)*fwtbdry[threadIdx.x];
	}
	else
	{
//		flux_bdry_right[threadIdx.x] = boundflux[0] - flux;
		flux_bdry_right[threadIdx.x] = (flux_buffer[0] - flux)*fwtbdry[threadIdx.x];
	}
}

__global__ void data_process_fb_left(double* cu_icresponse,double* r_bdry,double* z_bdry,double* bdry_response,int* numbdry_gafile,double* fwtbdry)
{
	int rID,zID,n1,n2,n3,n4,num;
	double m,m1,m2,m3,m4,r,z,r0,z0; 

	__shared__ double flux_buffer[1];

	r = r_bdry[threadIdx.x];
	z = z_bdry[threadIdx.x];
	num = numbdry_gafile[0]-1;

	rID = int((r-igridr)/delR);
	zID = int((z-igridz)/delZ);

	r0 = rID*delR+igridr;
	z0 = zID*delZ+igridz;

	m1 = (r0+delR-r)*(z0+delZ-z);
	m2 = (r-r0)*(z0+delZ-z);
	m3 = (r0+delR-r)*(z-z0);
	m4 = (r-r0)*(z-z0);

	n1 = rID*mesh2+zID;
	n2 = (rID+1)*mesh2+zID;
	n3 = rID*mesh2+zID+1;
	n4 = (rID+1)*mesh2+zID+1;

	m = (m1*cu_icresponse[n1*EF_NFCOIL+blockIdx.x] + m2*cu_icresponse[n2*EF_NFCOIL+blockIdx.x] + m3*cu_icresponse[n3*EF_NFCOIL+blockIdx.x] + m4*cu_icresponse[n4*EF_NFCOIL+blockIdx.x])/darea;

	if (threadIdx.x == 0)
	{
		bdry_response[blockIdx.x*num] = m*fwtbdry[threadIdx.x];
		flux_buffer[0] = m;
	}

	__syncthreads();

	if (threadIdx.x > 0)
	{
		bdry_response[blockIdx.x*num + threadIdx.x] = (m - flux_buffer[0])*fwtbdry[threadIdx.x];
	}
}

__global__ void changePos_1_fb(double* cu_dgfc, double* d_A, int* numbdry_gafile)
{
	int index;
       	index = blockIdx.x * (numbdry_gafile[0]-1) + threadIdx.x;
	d_A[index] = cu_dgfc[index];
}

__global__ void changePos_2_fb(double* d_U,double* d_S,double* d_VT,int* numbdry_gafile,double* cu_U_matrix,double* cu_S_vector,double* cu_VT_matrix)
{
	printf("changePos_2 start ! \n");
	int index1,index2,index3;
       	
	index1 = blockIdx.x * (numbdry_gafile[0]-1) + threadIdx.x;
	index2 = blockIdx.x * (EF_NFCOIL) + threadIdx.x;
	index3 = threadIdx.x;

	cu_U_matrix[index1] = d_U[index1];

	if(blockIdx.x<(EF_NFCOIL) && threadIdx.x<(EF_NFCOIL))
	{
		cu_VT_matrix[index2] = d_VT[index2];
	}
	else if (blockIdx.x == (EF_NFCOIL) && threadIdx.x<(EF_NFCOIL))
	{
		cu_S_vector[index3] = d_S[index3];
	}

//	printf("U matrix %d  %f  \n",index1,cu_U_matrix[index1]);
	printf("changePos_2 done ! \n");
}

__global__ void responsematrix_inverse_U_fb(double* m, double* U, float* current_representation, double* result)
{
	int index1, index2, num;
	__shared__ double temp[250];

	num = EF_DIAG+(int)current_representation[26];

	index1 = blockIdx.x * num + threadIdx.x;
	index2 = threadIdx.x;

	temp[index2] = U[index1] * m[index2];
//	if (blockIdx.x == 0)
//		printf("%d U'*m %f U' %f m %f \n",index1,temp[index2],U[index1],m[index2]);
	__syncthreads();

	for(unsigned int s=64;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<num) 
		{
			temp[threadIdx.x] += temp[k];
		}
		__syncthreads();
	}
	
	if (threadIdx.x == 0)
	{
		result[blockIdx.x] = temp[0];
//		printf("U'*m %d  %f  \n",blockIdx.x,temp[0]);
	}
//	printf("inverse1 done! \n");
}

__global__ void responsematrix_inverse_S_fb(double* result_temp, double* S, float* current_representation)
{
//	int num;
	double s_temp; //s_max,s_truncate;
//	__shared__ double temp[60];
//	__shared__ double UT[121];
	
//	printf("inverse 1 start!! \n");

//	num = (EF_NFCOIL+(int)current_representation[25]);

//	index1 = blockIdx.x * num + threadIdx.x;
//	index2 = threadIdx.x;

	s_temp = S[threadIdx.x];
//	s_max = S[0];
//	s_truncate = s_max * truncate_factor;
//	printf("S %d = %f \n",threadIdx.x,s_temp);
//	if (s_temp == 0 || s_temp < s_truncate)
//	{
//		result_temp[threadIdx.x] = 0;
//	}
//	else
//	{
	result_temp[threadIdx.x] = result_temp[threadIdx.x]/s_temp; 
//	}
//		printf("ST*U'*M %d = %f \n",threadIdx.x,result_temp[threadIdx.x]);
//	printf("inverse2 done! \n");
}

__global__ void responsematrix_inverse_VT_fb(double* result_temp, double* VT, double* test_right, float* current_representation)
{
	int index1, index2, num;
	__shared__ double temp[60];
//	__shared__ double UT[121];
	
//	printf("inverse 1 start!! \n");

	num = (EF_NFCOIL+(int)current_representation[25]);

	index1 = blockIdx.x * num + threadIdx.x;
	index2 = threadIdx.x;

	temp[index2] = VT[index1] * result_temp[index2];
	__syncthreads();

//	if (blockIdx.x == 0)
	//	printf("%d V* ST* U'*m %f VT %f ST* U'*m %f \n",index1,temp[index2],VT[index1],result_temp[index2]);
//		printf(" %f \n",temp[index2]);

	for(unsigned int s=16;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<num) 
		{
			temp[threadIdx.x] += temp[k];
		}
		__syncthreads();
	}
	__syncthreads();

	if (threadIdx.x == 0)
	{
		test_right[blockIdx.x] = temp[0];
//		printf("test_right %d = %f \n",blockIdx.x,test_right[blockIdx.x]);
	}

//	printf("inverse3 done! \n");
}

__global__ void data_process_fb(double* cu_result, double* ec_flux, double* ic_flux, double* r_bdry, double* z_bdry, double* flux_bdry_right, double* fwtbdry)
{
	int rID,zID,rIDe,zIDe;
        double flux,flux1,flux2,flux3,flux4,m1,m2,m3,m4,r,z,r0,z0; 

	__shared__ double flux_buffer[1];

	r = r_bdry[threadIdx.x];
	z = z_bdry[threadIdx.x];

	rID = int((r-igridr-delR)/delR);
	zID = int((z-igridz-delZ)/delZ);

	rIDe = int((r-igridr)/delR);
	zIDe = int((z-igridz)/delZ);

	r0 = (rID+1)*delR+igridr;
	z0 = (zID+1)*delZ+igridz;

	flux1 = cu_result[rID*mesh1+zID] + ec_flux[rIDe*mesh2+zIDe];
	flux2 = cu_result[(rID+1)*mesh1+zID] + ec_flux[(rIDe+1)*mesh2+zIDe];
	flux3 = cu_result[rID*mesh1+zID+1] + ec_flux[rIDe*mesh2+zIDe+1];
	flux4 = cu_result[(rID+1)*mesh1+zID+1] + ec_flux[(rIDe+1)*mesh2+zIDe+1];

	m1 = (r0+delR-r)*(z0+delZ-z);
	m2 = (r-r0)*(z0+delZ-z);
	m3 = (r0+delR-r)*(z-z0);
	m4 = (r-r0)*(z-z0);

	flux = (m1*flux1+m2*flux2+m3*flux3+m4*flux4)/darea;

	if (threadIdx.x == 0)
	{
		flux_buffer[0] = flux;

		flux1 = ic_flux[rIDe*mesh2+zIDe];
		flux2 = ic_flux[(rIDe+1)*mesh2+zIDe];
		flux3 = ic_flux[rIDe*mesh2+zIDe+1];
		flux4 = ic_flux[(rIDe+1)*mesh2+zIDe+1];

		flux = (m1*flux1+m2*flux2+m3*flux3+m4*flux4)/darea;
	}
	
	__syncthreads();

	if (threadIdx.x == 0)
	{
		flux_bdry_right[threadIdx.x] = flux*fwtbdry[threadIdx.x];
	}
	else
	{
//		flux_bdry_right[threadIdx.x] = boundflux[0] - flux;
		flux_bdry_right[threadIdx.x] = (flux_buffer[0] - flux)*fwtbdry[threadIdx.x];
	}
}

__global__ void	pfFlux_fb(double* right_side, double* cu_icresponse, double* cu_icflux)
{
	__shared__ double sresponse[48][18];
	unsigned int index = blockIdx.x * 864 + threadIdx.y*18 + threadIdx.x;
//	sresponse[threadIdx.y][threadIdx.x] = ;
//	__syncthreads();
	sresponse[threadIdx.y][threadIdx.x] = right_side[threadIdx.x] * cu_icresponse[index]; 
	__syncthreads();

//	if(blockIdx.x == 0 && threadIdx.y == 0)
//	{
//		printf("pfflux %d %f \n",threadIdx.x,right_side[threadIdx.x]);
//	}

//	double temp_sum = 0;

//	if(threadIdx.x<6)
//	{
//		temp_sum +=sresponse[threadIdx.y][threadIdx.x]+sresponse[threadIdx.y][threadIdx.x+6]+sresponse[threadIdx.y][threadIdx.x+12];
//		__syncthreads();
//		sresponse[threadIdx.y][threadIdx.x] = temp_sum;		   
//	}
//	__syncthreads();

	if(threadIdx.x == 0)
	{
		for (int i=1;i<EF_NFCOIL;i++)
		{
			sresponse[threadIdx.y][0] = sresponse[threadIdx.y][0] + sresponse[threadIdx.y][i];
		}
		__syncthreads();
		cu_icflux[blockIdx.x*48+threadIdx.y] = sresponse[threadIdx.y][0];
//		printf("pf->grid psi  %d = %f \n", blockIdx.x*48+threadIdx.y, sresponse[threadIdx.y][0]);
	}

}

__global__ void	Fluxadd_fb(double* cu_icflux,double* cu_ecflux,double* cu_ipflux,double* cu_result)
{
	if(blockIdx.x>0&&blockIdx.x<256&&threadIdx.x>0&&threadIdx.x<256)
	{
		cu_result[(blockIdx.x-1)*256+threadIdx.x-1] = cu_icflux[blockIdx.x*257+threadIdx.x]+cu_ecflux[blockIdx.x*257+threadIdx.x]+cu_ipflux[(blockIdx.x-1)*256+threadIdx.x-1];
//		cu_result[(blockIdx.x-1)*256+threadIdx.x-1] = -1* cu_result[(blockIdx.x-1)*256+threadIdx.x-1];
	}
}

__global__ void expansion_finder_fb(double *cu_result,double* leas_flux)
{
	double br,bz,c,d,invr;
	__shared__ double sDiff[4][63];
	__shared__ double position[63];
	position[threadIdx.x] = double(threadIdx.x);
	sDiff[3][threadIdx.x] = 9999;
	__syncthreads();
	if (blockIdx.x<122)
	{
		invr=1/cu_radix_fb[blockIdx.x+23];
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
		invr=1/cu_radix_fb[blockIdx.x+23-122];
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

__global__ void norm1_fb(double* leas_flux)
{
	double c, d;
	if(blockIdx.x == 0)
	{
		__shared__ double flux[122];
		__shared__ double fflux[2][122];
		flux[threadIdx.x] = leas_flux[threadIdx.x+3000];        // si
		fflux[0][threadIdx.x] = leas_flux[threadIdx.x+3000+122];    //  z
		fflux[1][threadIdx.x] = double(threadIdx.x);  		// r   
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
		__shared__ double bp[122];
		__shared__ double bpp[2][122];
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
		__shared__ double bp[122];
		__shared__ double bpp[2][122];
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

__global__ void xpoint_brbz_fb(double* leas_flux, double* cu_result)
{
	int rID, zID;
	double br,bz,invr;
	__shared__ double tempbr[5];
	__shared__ double tempbz[5];
	if(blockIdx.z == 0)
	{
		rID = int(leas_flux[1])+23-1+blockIdx.x;
		zID = int(leas_flux[0])-1+blockIdx.y;
		invr = 1/cu_radix_fb[rID];

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
		invr = 1/cu_radix_fb[rID];

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
		invr = 1/cu_radix_fb[rID];

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

__global__ void post_xpoint_brbz_fb(double* leas_flux)
{
	__shared__ double br[4];
	__shared__ double bz[4];

	br[threadIdx.x] = leas_flux[4000+100*blockIdx.x+threadIdx.x];
	bz[threadIdx.x] = leas_flux[5000+100*blockIdx.x+threadIdx.x];

	double temp;

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


__global__ void xpoint_pos_fb(double* leas_flux, double* cu_result, double* xpoint, double* cu_xlocation)
{
	double ebr, ebz, ebrr, ebrz, ebzr, ebzz, dr, dz, er, ez, xr, xz, xflux;
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
		xflux += cu_radix_fb[rID+23]*(-dz*ebr+ebz*dr);
//		printf("\ndown flux new %f \n", xflux);
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
		xflux += cu_radix_fb[rID+23]*(-dz*ebr+ebz*dr);
		xpoint[2] = xflux;
		xpoint[7] = xr ;	
		xpoint[8] = xz ;	
	//	printf("\nup flux new %f \n", xflux);
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
		xflux += cu_radix_fb[rID+68]*(-dz*ebr+ebz*dr);
		xpoint[4] = xflux;
		xpoint[9] = er ;	
		xpoint[10] = xz ;	
//		printf("\nmid flux new %f \n", xflux);
		cu_xlocation[4] = double(zID);
		cu_xlocation[5] = double(rID);
	}
}

__global__ void limiterFlux_fb(float* limiter, double* cu_result, double* xpoint,double* cuLimFlux)
{
	int rID,zID;
	double flux1,flux2,flux3,flux4;
	double m1,m2,m3,m4;
	double r,z,zx1,zx2,r0,z0;	
	
	r = double(limiter[threadIdx.x]);
	z = double(limiter[threadIdx.x + LIMITER_NUM]);

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

__global__ void norm3_fb(double* xpoint,double*cuLimFlux)
{
	__shared__ double slimFlux[LIMITER_NUM];
	double c,d;
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

__global__ void data_eli_fb(double* cu_result,double* effect_flux,int* cu_num_se,double* ffprime,double* pprime,double* r_bdry,double* z_bdry,double* xpoint,double* cuLimFlux,double* cu_bound_flux)
{
	int offset,a,b,c,d,e,f,g;
	double fluxaxis,fluxbound,normfactor,interfactor,x;

	__shared__ double ffprime_buffer[mesh2];
	__shared__ double pprime_buffer[mesh2];

	__shared__ double buffer[mesh][2];
	__shared__ double normedflux[mesh];

	ffprime_buffer[threadIdx.x] = ffprime[threadIdx.x];
	pprime_buffer[threadIdx.x] = pprime[threadIdx.x];

	if (threadIdx.x < 2)
	{
		ffprime_buffer[threadIdx.x+mesh] = ffprime[threadIdx.x+mesh];
		pprime_buffer[threadIdx.x+mesh] = pprime[threadIdx.x+mesh];
	}

	if(cuLimFlux[0]>xpoint[1])
		fluxbound = cuLimFlux[0];
	else
		fluxbound = xpoint[1];

	fluxaxis = xpoint[4];
	
	if(threadIdx.x==0&&blockIdx.x == 0)
	{
//		fluxbound = 0.3375518918;
		cu_bound_flux[0] = fluxbound;
		xpoint[11] = fluxaxis-fluxbound;
	}

//	if(threadIdx.x==0&&blockIdx.x == 0)              
//		printf("lim%20.18lf  bdry%20.18lf  ax%20.18lf  diff%20.18lf\n",cuLimFlux[0],cu_bound_flux[0],xpoint[4],xpoint[11]);	

	unsigned int index=blockIdx.x*(mesh+1)+threadIdx.x;
	buffer[threadIdx.x][1]=cu_result[index];
	normfactor=1/(fluxbound-fluxaxis);

	c = (xpoint[6]-igridz-delZ)/delZ;
	d = (xpoint[8]-igridz-delZ)/delZ;
	e = (1.016 - igridr - delR)/delR;
	f = (2.3540 - igridr - delR)/delR;

	__syncthreads();

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
		normedflux[b]=(buffer[b][0]-fluxaxis)*normfactor;
		if(normedflux[b]>=0&&normedflux[b]<=1&&a>(c-1)&&a<(d+1)&&blockIdx.x>(e-1)&&blockIdx.x<(f+1))
		{
			g = int(normedflux[b]*mesh1);
			interfactor = normedflux[b]*mesh1 - g;

			x = cu_radix_fb[blockIdx.x]*((1-interfactor)*pprime_buffer[g] + interfactor*pprime_buffer[g+1]) + ((1-interfactor)*ffprime_buffer[g] + interfactor*ffprime_buffer[g+1])/cu_radix_fb[blockIdx.x];
//			printf("blockIdx.x %d threadIdx.x %d interfactor %f x %f \n",blockIdx.x, threadIdx.x, interfactor, x);
		}
		else
		{
			x = 0;
		}
		//begin to store the transposed data into global memory
		
		effect_flux[offset+b] = x;			
	}
}

__global__ void compute_converge_init_fb(double* cu_result,double* cu_result_compare,int* converge)
{
	unsigned int index1 = blockIdx.x * mesh2 + threadIdx.x;
	unsigned int index2 = blockIdx.x * mesh1 + threadIdx.x;
	cu_result_compare[index2] = cu_result[index1];

	if(threadIdx.x == 0 && blockIdx.x == 0)
	{
		converge[0] = 0;	
	}
}

__global__ void compute_converge_fb(double* cu_result,double* cu_result_compare)
{
	__shared__ double buffer[mesh];
	double c,d;
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

__global__ void compute_converge_post_fb(double* cu_result_compare,double* xpoint,double* cu_bound_flux,int* converge,double* leas_flux)
{
	__shared__ double buffer[mesh];
	double c,d;
	double error;

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
	
	error = 0.0000001*(xpoint[4] - cu_bound_flux[0]);

	if(threadIdx.x == 0)
	{
		if(buffer[0] < error)
		{
			converge[0] = 1;
			leas_flux[0] = buffer[0];
			printf("error %f ",buffer[0]/error);
//			printf("%f %f\n",error,leas_flux[907]);
		}
		else
		{
			converge[0] = 0;

//			printf("error %f \n",buffer[0]/error);
//			printf("%d \n",converge[0]);
//			printf("%f %f\n",error,buffer[0]);
		}
	}
}

__global__ void result_regroup_fb(double* cu_result,double* result,double* cu_bound_resut,double* cu_icflux,double* cu_ecflux)
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

__global__ void error_amend(double* effect_flux, double* leas_flux)
{
	__shared__ double plasmaCurrent[32];
	unsigned int index = blockIdx.x * 32 + threadIdx.x;
	if(blockIdx.x<ELI_CORE)
	{
		plasmaCurrent[threadIdx.x] = effect_flux[index];
		for(unsigned int s=16;s>0;s>>=1)
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<32) 
			{
				plasmaCurrent[threadIdx.x] += plasmaCurrent[k];
			}
		}
	}
	else 
	{
		if(threadIdx.x<ELI_REST)		
			plasmaCurrent[threadIdx.x] = effect_flux[index];
	
		for(unsigned int s=16;s>0;s>>=1)
		{
			int k=threadIdx.x+s;
			if (threadIdx.x<s&&k<ELI_REST) 
			{
				plasmaCurrent[threadIdx.x] += plasmaCurrent[k];
			}
		}

	}
	leas_flux[blockIdx.x] = plasmaCurrent[0];
} 


__global__ void error_amend2(double* leas_flux, double* init_diag_dat)
{
	__shared__ double plasmaCurrent[ELI_CORE];

	plasmaCurrent[threadIdx.x] = leas_flux[threadIdx.x];
	if (threadIdx.x+512 < ELI_CORE)
		plasmaCurrent[threadIdx.x+512] = leas_flux[threadIdx.x+512];

	for(unsigned int s=512;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<ELI_CORE) 
		{
			plasmaCurrent[threadIdx.x] += plasmaCurrent[k];
		}
		__syncthreads();
	}
	__syncthreads();
	
	init_diag_dat[EF_DIAG-1] = plasmaCurrent[0];
//	printf("%15.7e\n",leas_flux[0]);
} 

__global__ void searchbdry_fb(double* cu_result,double* leas_flux,double* xpoint,double* cu_bound_flux,double* tani,double* coti,double* xybdry,int* maxgridbdry)
{
	int ra,rb,za,zb,r_axi,z_axi;
	double si1,si2,dif,bdry,factor;
	__shared__ double xbdry[128];
	__shared__ double ybdry[128];
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

__global__ void postsearchbdry_fb(double* xybdry,double* leas_flux,double* xpoint,double* cu_bound_flux)
{
	__shared__ double rzbdry[numbdry*2];
	int z_axi,r_axi,num_u,num_b;
	double ra,za,rxu,zxu,rxb,zxb,agu,agb,rtempu,ztempu,agtempu,rtempb,ztempb,agtempb;

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

__global__ void search_surface_fb(double* cu_result,double* location,double* tani,double* coti,double* xpoint,int* num,double* xy_surface)
{
//	__shared__ int N[32];
	int ra,rb,za,zb,r_axi,z_axi,i,N1,N2;
	double si1,si2,flux,dflux,aflux,factor;

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

__global__ void q_profile_fb(double* cu_result,double* xpoint,double* xy_surface,double* btor,double* cu_qpsi)
{
	int ra,za;
	double r,z,factor1,factor2,bp,br,bz,bt,dl,a,b;
	double flux1,flux2,flux3,flux4,flux_a,flux_b,flux_c,flux_d;
	__shared__ double sum[numbdry];

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

__global__ void forbetapli_fb(double* cu_result,double* leas_flux,double* profile_result,double* cu_bound_flux,double* xpoint)
{
	double invr,fluxaxi,fluxbound,br,bz,normflux,factor,interfactor;
	int a,b,c,d,e;
	__shared__ double sDiff[7][mesh];
	__shared__ double profile[mesh2];

	invr = 1/cu_radix_fb[blockIdx.x];
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
		interfactor = double(normflux*256 - e); 

//		printf("e %d interfactor %f p %f \n ",e,interfactor,(interfactor*profile[e+1]+(1-interfactor)*profile[e]));

//		sDiff[3][threadIdx.x] = factor*2*pi*double(cu_radix[blockIdx.x])*(br*br + bz*bz);  

		sDiff[3][threadIdx.x] = factor*2*pi*cu_radix_fb[blockIdx.x]*(br*br + bz*bz); 
		sDiff[4][threadIdx.x] = factor*2*pi*(cu_radix_fb[blockIdx.x]) * (interfactor*profile[e+1]+(1-interfactor)*profile[e]);
		sDiff[5][threadIdx.x] = factor*2*pi*cu_radix_fb[blockIdx.x];
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

__global__ void postbetapli_fb(double* leas_flux,double* xybdry,double* cu_compute_diagnotics)
{
	__shared__ double rzbdry[(numbdry+3)*2];
	__shared__ double sumbdry[numbdry+2];
	__shared__ double sum[mesh][4];
	double a,b;
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

//		printf("volume %f betap %f cross_area %f li %f sumbdry %f energy %f ",leas_flux[590000],leas_flux[590001],leas_flux[590002],leas_flux[590003],sumbdry[0],leas_flux[590005]);

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
