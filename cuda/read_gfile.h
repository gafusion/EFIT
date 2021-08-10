#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int read_gafile()
{
	FILE *gfile=fopen("g133221.03000","r+");    

	islice = 1;
	cudaMalloc((void**)&cu_ishot,(islice*sizeof(float)));
	ishot = (float *)malloc(islice*sizeof(float));
    	
	itime = (float *)malloc(islice*sizeof(float));
	cudaMalloc((void**)&cu_itime,(islice*sizeof(float)));

	char a[30]="";
	char mystring[20];
	char FL[70];
	char aa[10];
	char aaa[10];
	
	fgets(FL,70,gfile);
	
//	printf("FL %s \n",FL);

	char *p;
	p = strtok(FL," ");
//	printf("p0 %s \n",p);
	
	int i,j;
	for(j=0;j<6;j++)
	{
		p=strtok(NULL," ");
//		printf("line1 information #%d %s \n",j+1,p);
		if (j == 1)
		{
			strncpy(aa,p+1,7);
			ishot[0] = float(atof(aa));
//			printf("shot = %f\n",ishot[0]);
		}
		if (j == 2)
		{
			strncpy(aaa,p,4);
			itime[0] = float(atof(aaa));
//			printf("time = %f",itime[0]);
		}
	}
	 
	grid_gafile = atoi(p);
//	printf("grid %d \n",grid_gafile);
	
	cudaMemcpy(cu_ishot,ishot,islice*sizeof(float),cudaMemcpyHostToDevice);
    	cudaMemcpy(cu_itime,itime,islice*sizeof(float),cudaMemcpyHostToDevice);

	plasma_parameter_gafile = (double *)malloc(20*sizeof(double));
	ffprime_gafile = (double *)malloc(grid_gafile*sizeof(double));
	pprime_gafile = (double *)malloc(grid_gafile*sizeof(double));
	fpol_gafile = (double *)malloc(grid_gafile*sizeof(double));
	pres_gafile = (double *)malloc(grid_gafile*sizeof(double));
	psi_gafile = (double *)malloc(grid_gafile*grid_gafile*sizeof(double));
	q_gafile = (double *)malloc(grid_gafile*sizeof(double));
	diagnostic1_fb = (double *)malloc(islice*EF_DIAG*sizeof(double));
	current_gafile = (double *)malloc(1*sizeof(double));

///////////////////////////// read plasma parameter/////////////////
	for(i=0;i<20;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				plasma_parameter_gafile[i] = atof(a);
			}
		} 
	}
	current_gafile[0] = plasma_parameter_gafile[10];
//	printf("plasma current % 10.9E \n",current_gafile[0]);
//	printf("psi boundary % 10.9E \n",plasma_parameter_gafile[8]);
//	printf("psi axis % 10.9E \n",plasma_parameter_gafile[7]);

//	float* abflux_gafile;
//	abflux_gafile = (float *)malloc(2*sizeof(float));
//	abflux_gafile[0] = plasma_parameter_gafile[7];
//	abflux_gafile[1] = plasma_parameter_gafile[8];
//	cudaMalloc((void**)&cu_abflux_gafile,(2*sizeof(float)));
//	cudaMemcpy(cu_abflux_gafile, abflux_gafile, 2*sizeof(float), cudaMemcpyHostToDevice);
	
///////////////////////////////read fpol/////////////////////////////////////////////
    	for(i=0;i<grid_gafile;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				fpol_gafile[i] = atof(a);
//				printf("fpol %d % 10.9E \n",i,fpol[i]);
			}
		} 
	}
//////////////////////////////read pres////////////////////////////////////////
    	for(i=0;i<grid_gafile;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				pres_gafile[i] = atof(a);
//				printf("pres %d % 10.9E \n",i,pres[i]);
			}
		} 
	}
////////////////////////////////read ffprime////////////////////////////////////////
    	for(i=0;i<grid_gafile;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				ffprime_gafile[i] = atof(a);
//				printf("ffprime %d % 10.9E \n",i,ffprime[i]);
			}
		} 
	}
//////////////////////////////read pprime////////////////////////////////////////
    	for(i=0;i<grid_gafile;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				pprime_gafile[i] = atof(a);
//				printf("pprime %d % 10.9E \n",i,pprime[i]);
			}
		} 
	}
///////////////////////////////////grid psi////////////////////////////////////
    	for(i=0;i<grid_gafile*grid_gafile;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				psi_gafile[i] = atof(a);
//				printf("psi %d % 10.9E \n",i,psi[i]);
			}
		} 
	}
///////////////////////////////////q profile/////////////////////////////////
	for(i=0;i<grid_gafile;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				q_gafile[i] = atof(a);
//				printf("q %d % 10.9E \n",i,q[i]);
			}
		} 
	}
//////////////////// boundary and limiter number
	if ( fgets (mystring , 17, gfile) != NULL )	
	{
		if(mystring[0]=='\n')
		{
			fgets (mystring , 17, gfile);
		}
//		printf("%s \n",mystring);
		p = strtok(mystring," ");
		numbdry_gafile = atoi(p);
		p = strtok(NULL," ");
		numlimter_gafile = atoi(p);
//		printf("numbdry %d numlimter %d \n",numbdry_gafile,numlimter_gafile);

		bdry_gafile = (double *)malloc(numbdry_gafile*2*sizeof(double));
		limiter_gafile = (double *)malloc(numlimter_gafile*2*sizeof(double));
	}
	
///////////////////// boundary points
	for(i=0;i<numbdry_gafile*2;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				bdry_gafile[i] = atof(a);
//				printf("bdry %d % 10.9E \n",i,bdry[i]);
			}
		} 
	} 
//////////////////// limiter points
	for(i=0;i<numlimter_gafile*2;i++)
	{
		if ( fgets (mystring , 17, gfile) != NULL )	
		{ 
			if(mystring[0]=='\n')
			{
				i=i-1;
			}
			else
			{
				strcpy(a,mystring);
				limiter_gafile[i] = atof(a);
//				printf("limiter %d % 10.9E \n",i,limiter[i]);
			}
		} 
	}

/////////////////read afile///////////////////////////////////////
	FILE *afile=fopen("a133221.03000","r+");    

	host_btor_fb = (double *)malloc(1*sizeof(double));
	cudaMalloc((void**)&test_btor_fb,(1*sizeof(double)));
		
	for (i=0;i<26;i++)
	{
		fgets(FL,70,afile);
		if (i == 4)
		{
			strncpy(a,FL+33,17);
			host_btor_fb[0] = atof(a);
//			printf("BTOR = % 10.9E \n",host_btor_fb[0]);
			cudaMemcpy(test_btor_fb, host_btor_fb, 1*sizeof(double), cudaMemcpyHostToDevice);

			strncpy(a,FL+49,17);
			diagnostic1_fb[EF_DIAG-1] = atof(a);
//			printf("diagnostic1_fb[120] = % 10.9E \n",diagnostic1_fb[EF_DIAG-1]);
	//		cudaMemcpy(test_btor_fb, host_btor_fb, 1*sizeof(double), cudaMemcpyHostToDevice);
		}
	}
//	printf("FL %s \n",FL);

	p = strtok(FL," ");
//	printf("afile line26 information #%d %s \n",1,p);
	for(j=0;j<3;j++)
	{
		p=strtok(NULL," ");
//		printf("afile line26 information #%d %s \n",j+2,p);
	}

	num_ecurrent_gafile = atoi(p);
	ecurrent_gafile = (double *)malloc(num_ecurrent_gafile*sizeof(double));

	for (i=0;i<35;i++)
	{
		fgets(FL,70,afile);
	}

	fgets (mystring , 18, afile);
	strcpy(a,mystring);
	ecurrent_gafile[0] = atof(a);
//	printf("ecurrent %d % 10.9E \n",0,ecurrent_gafile[0]);
	for(j=0;j<num_ecurrent_gafile-1;j++)
	{
		if ( fgets (mystring , 17, afile) != NULL )	
		{
			if(mystring[0]=='\n')
			{
//				j=j-1;
				fgets (mystring , 18, afile);
				strcpy(a,mystring);
				ecurrent_gafile[j+1] = atof(a);
//				printf("ecurrent %d % 10.9E \n",j+1,ecurrent_gafile[j+1]);
			}
			else
			{
				strcpy(a,mystring);
				ecurrent_gafile[j+1] = atof(a);
//				printf("ecurrent %d % 10.9E \n",j+1,ecurrent_gafile[j+1]);
			}
		}
	}

	return 1;
}

int read_gafile_post()
{
	double *psi_gafile_temp;
	double *pprime_gafile_temp;
	double *ffprime_gafile_temp;
	double *psibdry_gafile_temp;
	double *psiaxis_gafile_temp;
	int *temp_numbdry_gafile;
	double* temp_ecurrent_gafile;

 	psi_gafile_temp = (double *)malloc(mesh2*mesh2*sizeof(double));
	pprime_gafile_temp = (double *)malloc(mesh2*sizeof(double));
	ffprime_gafile_temp = (double *)malloc(mesh2*sizeof(double));
	psibdry_gafile_temp = (double *)malloc(1*sizeof(double));
	psiaxis_gafile_temp = (double *)malloc(1*sizeof(double));
	temp_numbdry_gafile = (int *)malloc(1*sizeof(int));
	temp_ecurrent_gafile = (double *)malloc(num_ecurrent_gafile*sizeof(double));

	r_bdry_gafile = (double *)malloc((numbdry_gafile+4)*sizeof(double));
	z_bdry_gafile = (double *)malloc((numbdry_gafile+4)*sizeof(double));
	fwtbdry = (double *)malloc((numbdry_gafile-1)*sizeof(double));

	temp_profile_result_fb = (double *)malloc(4*mesh2*sizeof(double));

	double interfactor = mesh1/(grid_gafile-1);
	int i,j,ii,jj,nrmax,nrmin,nzmax,nzmin;
	double r,z,rr,zz,dr,dz,factor,factor1,factor2,factor3,factor4;
//	printf("interfactor = %f \n",interfactor);

	if (current_gafile[0] > 0)
	{
		factor = -1;
	}
	else
	{
		factor = 1;
	}

	psibdry_gafile_temp[0] = factor*plasma_parameter_gafile[8];
	psiaxis_gafile_temp[0] = factor*plasma_parameter_gafile[7];

	if (interfactor == 1)
	{
		for (i=0;i<mesh2;i++)
		{
			for(j=0;j<mesh2;j++)
			{
				psi_gafile_temp[i*mesh2+j] = factor*psi_gafile[j*grid_gafile+i];
			}
		}
		for (i=0;i<mesh2;i++)
		{
			ffprime_gafile_temp[i] = factor*ffprime_gafile[i]*(darea*10000000)/(4*pi);
			pprime_gafile_temp[i] = factor*pprime_gafile[i]*darea;
			temp_profile_result_fb[i] = factor*fpol_gafile[i]; //fpol profile
			temp_profile_result_fb[i+mesh2] = pres_gafile[i]; //pres profile
			temp_profile_result_fb[i+mesh2*2] = factor*ffprime_gafile[i]; //ffprime profile
			temp_profile_result_fb[i+mesh2*3] = factor*pprime_gafile[i]; //pprime profile
		}
	}
	if (interfactor > 1)
	{
		for (i=0;i<mesh2;i++)
		{
			for(j=0;j<mesh2;j++)
			{
				r = igridr + i*delR;
				z = igridz + i*delZ;

				ii = int(i/interfactor);
				jj = int(j/interfactor);

				rr = igridr + (interfactor*delR)*ii; 
				zz = igridz + (interfactor*delZ)*jj;

				dr = r-rr;
				dz = z-zz;

				factor3 = (dr*dz)/(darea*interfactor*interfactor);
				factor2 = ((interfactor*delR-dr)*dz)/(darea*interfactor*interfactor);
				factor4 = (dr*(interfactor*delZ-dz))/(darea*interfactor*interfactor);
				factor1 = ((interfactor*delR-dr)*(interfactor*delZ-dz))/(darea*interfactor*interfactor);

				psi_gafile_temp[i*mesh2+j] = factor*(factor1*psi_gafile[jj*grid_gafile+ii] + factor2*psi_gafile[(jj+1)*grid_gafile+ii] + factor3*psi_gafile[(jj+1)*grid_gafile+ii+1] + factor4*psi_gafile[jj*grid_gafile+ii+1]);
			}
		}
	}
	if (interfactor < 1)
	{
		for (i=0;i<mesh2;i++)
		{
			for(j=0;j<mesh2;j++)
			{
				psi_gafile_temp[i*mesh2+j] = psi_gafile[int((j/interfactor)*grid_gafile+(i/interfactor))];
			}
		}
	}

	cudaMalloc((void**)&cu_ffprime_gafile,(mesh2*sizeof(double)));
	cudaMalloc((void**)&cu_pprime_gafile,(mesh2*sizeof(double)));
	cudaMalloc((void**)&cu_psi_gafile,(mesh2*mesh2*sizeof(double)));
	cudaMalloc((void**)&cu_psibdry_gafile,(1*sizeof(double)));
	cudaMalloc((void**)&cu_psiaxis_gafile,(1*sizeof(double)));

	cudaMemcpy(cu_ffprime_gafile, ffprime_gafile_temp, mesh2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_pprime_gafile, pprime_gafile_temp, mesh2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_psi_gafile, psi_gafile_temp, mesh2*mesh2*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_psibdry_gafile, psibdry_gafile_temp, 1*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_psiaxis_gafile, psiaxis_gafile_temp, 1*sizeof(double), cudaMemcpyHostToDevice);

	r_bdry_gafile[numbdry_gafile] = bdry_gafile[0];
	r_bdry_gafile[numbdry_gafile+1] = bdry_gafile[0];
	r_bdry_gafile[numbdry_gafile+2] = bdry_gafile[0];
	r_bdry_gafile[numbdry_gafile+3] = bdry_gafile[0];

	z_bdry_gafile[numbdry_gafile] = bdry_gafile[0];
	z_bdry_gafile[numbdry_gafile+1] = bdry_gafile[0];
	z_bdry_gafile[numbdry_gafile+2] = bdry_gafile[0];
	z_bdry_gafile[numbdry_gafile+3] = bdry_gafile[0];

	nrmax = 0;
	nrmin = 0;
	nzmax = 0;
	nzmin = 0;
	
	for (i=0;i<numbdry_gafile;i++)
	{
		r_bdry_gafile[i] = bdry_gafile[i*2];
		z_bdry_gafile[i] = bdry_gafile[i*2+1];

//		printf("r %f, z %f\n",r_bdry_gafile[i],z_bdry_gafile[i]);
		
		if (r_bdry_gafile[i] > r_bdry_gafile[numbdry_gafile+3])
		{
			r_bdry_gafile[numbdry_gafile+3] = r_bdry_gafile[i];
			z_bdry_gafile[numbdry_gafile+3] = z_bdry_gafile[i];
			nrmax = i;
//			printf("nrmax = %d \n",i);

		}
		if (r_bdry_gafile[i] < r_bdry_gafile[numbdry_gafile+1])
		{
			r_bdry_gafile[numbdry_gafile+1] = r_bdry_gafile[i];
			z_bdry_gafile[numbdry_gafile+1] = z_bdry_gafile[i];
			nrmin = i;
//			printf("nrmin = %d \n",i);
		}

		if (z_bdry_gafile[i] > z_bdry_gafile[numbdry_gafile+2])
		{
			z_bdry_gafile[numbdry_gafile+2] = z_bdry_gafile[i];
			r_bdry_gafile[numbdry_gafile+2] = r_bdry_gafile[i];
			nzmax = i;
//			printf("nzmax = %d \n",i);
		}
		if (z_bdry_gafile[i] < z_bdry_gafile[numbdry_gafile])
		{
			z_bdry_gafile[numbdry_gafile] = z_bdry_gafile[i];
			r_bdry_gafile[numbdry_gafile] = r_bdry_gafile[i];
			nzmin = i;
//			printf("nzmin = %d \n",i);
		}
	}
//	printf("nrmax = %d  nrmin = %d  nzmax = %d  nzmin = %d\n",nrmax,nrmin,nzmax,nzmin);
	for(i=0;i<numbdry_gafile-1;i++)
	{
		fwtbdry[i] = 1;
	}
	fwtbdry[nzmin] = 10;
	fwtbdry[nrmin] = 10;
	fwtbdry[nzmax] = 10;
	fwtbdry[nrmax] = 10;

//	printf("r_max %f r_min %f z_max %f z_min %f \n", r_bdry_gafile[numbdry_gafile],r_bdry_gafile[numbdry_gafile+1],z_bdry_gafile[numbdry_gafile],z_bdry_gafile[numbdry_gafile+1]);
//
	cudaMalloc((void**)&cu_r_bdry_gafile,((numbdry_gafile+4)*sizeof(double)));
	cudaMalloc((void**)&cu_z_bdry_gafile,((numbdry_gafile+4)*sizeof(double)));
	cudaMalloc((void**)&cu_fwtbdry,((numbdry_gafile-1)*sizeof(double)));

	cudaMemcpy(cu_r_bdry_gafile, r_bdry_gafile, (numbdry_gafile+4)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_z_bdry_gafile, z_bdry_gafile, (numbdry_gafile+4)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_fwtbdry, fwtbdry, (numbdry_gafile-1)*sizeof(double), cudaMemcpyHostToDevice);

	temp_numbdry_gafile[0] = numbdry_gafile;
	cudaMalloc((void**)&cu_numbdry_gafile,1*sizeof(int));
	cudaMemcpy(cu_numbdry_gafile, temp_numbdry_gafile, 1*sizeof(int), cudaMemcpyHostToDevice);

	for(i=0;i<num_ecurrent_gafile;i++)
	{
		temp_ecurrent_gafile[i] = ecurrent_gafile[i];
	}
	cudaMalloc((void**)&cu_ecurrent_gafile,num_ecurrent_gafile*sizeof(double));
	cudaMemcpy(cu_ecurrent_gafile, temp_ecurrent_gafile, num_ecurrent_gafile*sizeof(double), cudaMemcpyHostToDevice);
	
	return 1;
}
