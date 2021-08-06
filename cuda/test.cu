#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "cutil.h"
#include <cublas_v2.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sched.h>
#include "lapacke.h"
#include <unistd.h>
#include "macro.h"
#include <pthread.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <cusolverDn.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#include <helper_cuda.h>
#include <helper_functions.h>

__constant__ int cu_num_block[mesh];
__constant__ int cu_blo_accu[mesh];
__constant__ float cu_currentInver[mesh]; 
__constant__ float cu_radix[mesh];

__constant__ double cu_currentInver_fb[mesh]; 
__constant__ double cu_radix_fb[mesh];

pthread_mutex_t mut;
pthread_t read_sig_t;
pthread_t pefit_run_t;

float* superb;
int* cu_num_se;
float* cu_limiter;
float* b_green_up;
float* b_green_down;
float* b_green_mid;
float* pf_green_up;
float* pf_green_down;
float* pf_green_mid;
float* ec_green_up;
float* ec_green_down;
float* ec_green_mid;
float* cu_efDiagresponse;
float* cu_ecgridresponse;
float* cu_gfc;
float* cu_ec_g;
float* cu_icresponse;
float* cu_boundresponse;
float* cu_initNormFlux;
float* cu_weight;
float* init_diag_dat;
float* cu_diagnostic1;
float* cu_diagnostic2;
float* cu_pro_diag;
float* test_diagnostic1;
float* test_diagnostic2;
float* test_psiref;
float* test_diagbit;
float* test_weight;
float* test_ec;
float* cu_diagbit1;
float* slice;
int islice;
float* ishot;
float* itime;
float* cu_slice;
float* cu_ishot;
float* cu_itime;
float* cu_psiref;
float* cu_mat;
float* responsematrix;
float* cu_act_trace;
float* cu_tri_main;
float* temp_cu_result;
float* temp_dmate;
float* cu_bfft;
float* cu_result;
float* cu_result_compare;
float* redu;
float* xpoint;
float* leas_flux;
float* fwt;
float* cubuffer;
float* cuLbuffer;
float* cu_dgfc;
float* cu_icflux;
float* cu_ecflux;
float* right_side;
float* test_right;
float* dmata;
float* hmata;
float* hright;
float* cuBoundResult;
float* effect_flux;
float* cuLimFlux;
float* cu_tani;
float* cu_coti;
float* xybdry;
int* limZ;
float* cu_bound_flux;
float* cu_xlocation;
float* time_values;
int* ipiv;
int info;
int small_loop = 0;
int time_key=0;
float* host_buffer;
float* de_buffer;
int num_shapes = 0;
int num_seg = 0;
float* keyTimes;
int* cu_inter_pos;
float* cu_inter_arg;
float* cu_grid_pos;
float* out_buf;
float* cu_out_buf;
int time_offset;
int num_grid = 0;
float next_ctime;
int s_key = 0;
int no_slice;
int slice_end;
int *cu_converge;
int *converge;
float cTime;
float * diagnostic2;
int* maxgridbdry;
float* xysurface;
float* cu_qpsi;
float* cu_compute_diagnotics;
float* cu_converge_error;
float* btor;
float* cu_btor;
float* test_btor;
float* host_btor;
float* diagnostic1;

float* result_psi;
float* result_bound;
float* result_rzxpt;
float* result_brsp;
float* result_axis;
float* result_nullpoint;
float* cu_result_output;
float* result_current;
float* cu_profile_result;
float* result_fpol_profile;
float* result_pres_profile;
float* result_ffprime_profile;
float* result_pprime_profile;
float* result_xybdry;
float* result_xylimtr;
float* cu_pcurrent;
float* result_pcurrent;
float* result_xysurface;
float* result_qpsi;
float* result_parameter;
float* compute_diagnotics;
float* converge_error;

float* current_representation;
float* cu_current_representation;

float* h_A;
float* h_U;
float* h_S;
float* h_VT;
float* d_A;
float* d_U;
float* d_S;
float* d_VT;
float* cu_U_matrix;
float* cu_S_vector;
float* cu_VT_matrix;
float* cu_response_inverse_temp;

int grid_gafile;
double* plasma_parameter_gafile;
double *ffprime_gafile;
double *pprime_gafile;
double *fpol_gafile;
double *pres_gafile;
double* psi_gafile;
double* q_gafile;
int numbdry_gafile;
int numlimter_gafile;
double* bdry_gafile;
double* limiter_gafile;
double* current_gafile;
int num_ecurrent_gafile; 
double* ecurrent_gafile;
double* r_bdry_gafile;
double* z_bdry_gafile;

double* cu_r_bdry_gafile;
double* cu_z_bdry_gafile;
double *cu_ffprime_gafile;
double *cu_pprime_gafile;
double *cu_psi_gafile;
double *cu_psiaxis_gafile;
double *cu_psibdry_gafile;
int* cu_numbdry_gafile;
double* cu_ecurrent_gafile;

double* cu_act_trace_fb;
double* cu_tri_main_fb;
double* cu_mat_fb;
double* cu_tani_fb;
double* cu_coti_fb;
double* cu_boundresponse_fb;
double* cu_icresponse_fb;
double* cu_ecgridresponse_fb;
double* effect_flux_fb;
double* cuLbuffer_fb;
double* cubuffer_fb;
double* cuBoundResult_fb;
double* cu_bfft_fb;
double* cu_ipflux_fb;
double* cu_ecflux_fb;
double* cu_pro_diag_fb;
double* cu_dgfc_fb;
double* superb_fb;
double* cu_U_matrix_fb;
double* cu_S_vector_fb;
double* cu_VT_matrix_fb;
double* cu_response_inverse_temp_fb;
double* test_right_fb;
double* h_A_fb;
double* h_U_fb;
double* h_S_fb;
double* h_VT_fb;
double* d_A_fb;
double* d_U_fb;
double* d_S_fb;
double* d_VT_fb;
double* cu_icflux_fb;
double* cu_result_fb;
double* leas_flux_fb;
double* xpoint_fb;
double* cu_xlocation_fb;
double* cuLimFlux_fb;
double* cu_bound_flux_fb;
double* cu_result_compare_fb;
double* cu_converge_error_fb;
double* fwtbdry;
double* cu_fwtbdry;

double* result_bound_fb;
double* result_nullpoint_fb;
double* temp_profile_result_fb;
double* cu_profile_result_fb;
double* result_fpol_profile_fb;
double* result_pres_profile_fb;
double* result_ffprime_profile_fb;
double* result_pprime_profile_fb;
double* cu_result_output_fb;
double* cu_compute_diagnotics_fb;
double* compute_diagnotics_fb;
double* xybdry_fb;
double* xysurface_fb;
double* host_btor_fb;
double* test_btor_fb;
double* cu_qpsi_fb;
double* result_qpsi_fb;
double* result_xybdry_fb;
double* result_parameter_fb;
double* converge_error_fb;
double* hright_fb;
double* result_psi_fb;
double* diagnostic1_fb;

int* num_knots;
float* knots_psi;
float* knots_pprime;
float* knots_ffprime;
float* grid_mse_bp;
float* pf_mse_bp;
float* ec_mse_bp;
float* num_channels_mse;
int* mpress;
int* mjphi;
float* cu_spline_psi;
float* cu_spline_left;
float* cu_spline_right;
float* cu_pres_psi;
float* cu_press;
float* cu_pres_left;
float* cu_pres_right;
float* cu_press_uncer;
int* cu_num_jphi;
float* cu_jphi_npsi;
float* cu_j_phi;
float* cu_j_phi_uncer;
float* cu_jphi_surface;
float* cu_jphi_left;
float* cu_jphi_right;
int* diag_num_all;
int* cu_diag_num_all;
int* offset_dA_mse;
int* offset_dA_jphi;
int* offset_dA_pres;
int* offset_dA_spline;
float* h_A_left;
float* cu_diag_right;

int work_size = 0;
int *devInfo;
double *work;

dim3 grid((mesh+1)/BLOCK_SIZE,(mesh+1)/BLOCK_SIZE,1);
dim3 threads(BLOCK_SIZE,BLOCK_SIZE,1);
dim3 grid2(mesh,1,1);
dim3 threads2(mesh,1,1);
dim3 grid3(1116,1,1);
dim3 threads3(32,32,1);
dim3 grid4(1,1,1);
dim3 threads4(EF_DIAG,POLY_NUM,1);
dim3 grid5(1,1,1);
dim3 threads5(UNKNOW_ALL,UNKNOW_ALL,1);
dim3 grid6(1116,1,1);
dim3 threads6(32,POLY_NUM,1);
dim3 grid7(1376,1,1);
dim3 threads7(18,48,1);
dim3 grid8(1116,1,1);
dim3 threads8(1024,1,1);
dim3 grid9(2,1,1);
dim3 threads9(15,20,1);
dim3 grid10(23,23,1);
dim3 threads10(32,1,1);
dim3 threads11(128,6);
dim3 grid17(512,1,1);
dim3 threads17(6,129,1);
dim3 gridbox(122,3,1);
#include "norm.h"
#include "xlocate.h"
#include "betali.h"
#include "write_gfile.h"
#include "read_gfile.h"
#include "fixbdry.h"

void num_regu(int* num,int*num_se,int* num_block,int* blo_accu)
{
	int i,j;
	int count=0;
	int blockcount=0;
//	int m;
	for(i=0;i<mesh;i++)
	{
		blockcount=0;
//		m=0;
		for(j=0;j<mesh;j++)
		{
			if(num[i*mesh+j]==1)
			{
				num_se[count]=j;
//				num_se[count+1]=m;
//				m++;
				blockcount++;
				count=count+1;
			}
		}
		num_block[i]=blockcount;
		if(i==0)
			blo_accu[i]=0;
		else
			blo_accu[i]=blo_accu[i-1]+num_block[i-1];
	}
//	for(int i=0;i<mesh;i++)
//		printf("%d\n",blo_accu[i]);

}

void execute_mode_init()
{
	////////// read current representation /////////////////////
	float* current_representation_temp;
	current_representation_temp = (float *)malloc(5000*sizeof(float));
	current_representation = (float *)malloc(30*sizeof(float));

	FILE* Number11;
	Number11=fopen("current_representation.dat","rt");
	if(Number11 == 0)
	{
		printf("error: can not open current_representation.dat\n");
//		return 1;
	}
	for(int i=0; i<1+4+1024*2; i++)
	{
		fscanf(Number11,"%f",current_representation_temp+i);
	}
	fclose(Number11);
	
	if (current_representation_temp[0] == 1 || current_representation_temp[0] == 3)
	{
		for (int i=0; i<30; i++)
		{
			if(i<5)
			{
				current_representation[i] = current_representation_temp[i];
			}
			else if (i<15)
			{
				current_representation[i] = 0.8*current_representation_temp[i];
			}
			else if (i<25)
			{
				current_representation[i] = 0.8*current_representation_temp[i+1024-10];
			}
			else
			{
				current_representation[i] = 0;
			}
		}
		current_representation[25] = current_representation[1] + current_representation[2];
		current_representation[26] = 2 - (current_representation[3] + current_representation[4]);

//	for(int i=0; i<30; i++)
//	{
//		printf("current_representation %d = %f \n",i,current_representation[i]);
//	}

	//////////////read gafile function////////////////////////
		if (current_representation[0] == 3)
		{
			read_gafile();
			read_gafile_post();

			current_representation[25] = 0;
			current_representation[26] = (numbdry_gafile-1)-EF_DIAG;
		}
	//////////////////////////////////////////////////////////
	}

	if (current_representation_temp[0] == 2)
	{
		read_gafile();
		read_gafile_post();

		num_knots = (int *)malloc(1*sizeof(int));
		knots_psi = (float *)malloc(10*sizeof(float));
		knots_pprime = (float *)malloc(10*sizeof(float));
		knots_ffprime = (float *)malloc(10*sizeof(float));

		for (int i=0; i<30; i++)
		{
			current_representation[i] = current_representation_temp[i];
		}

		num_knots[0] = current_representation[1];
		for (int i=0; i<num_knots[0]; i++)
		{
			knots_psi[i] = current_representation[2+i];
			printf("knots_psi %d = %f",i,knots_psi[i]);
			knots_pprime[i] = current_representation[2 + num_knots[0] + i];
			knots_ffprime[i] = current_representation[2 + num_knots[0]*2 + i];
		}
		
		current_representation[26] = 0;
		current_representation[25] = (num_knots[0]-1)*8;

		cudaMalloc((void**)&cu_spline_psi,num_knots[0]*sizeof(float));	
		cudaMemcpy(cu_spline_psi,knots_psi,num_knots[0]*sizeof(float),cudaMemcpyHostToDevice);
	}

	cudaMalloc((void**)&cu_current_representation,(30*sizeof(float)));
    	cudaMemcpy(cu_current_representation,current_representation,30*sizeof(float),cudaMemcpyHostToDevice);
}

void tri_solver_init()
{
	int i,j;
    	float tri_down[mesh-1];
   	float tri_up[mesh-1];
    	float tri_main[mesh*mesh];
    	float act_trace[mesh*mesh];
	float currentInver[mesh];
	for (i=0;i<(mesh-1);i++)
    	{
        	tri_down[i]=-1-delR/(2*(igridr+delR*(i+2)));
        	tri_up[i]=-1+delR/(2*(igridr+delR*(i+1)));
//		currentInver[i]=(igridr+delR*i)*delR*delR;
		currentInver[i]=(igridr+delR*(i+1))*delR*delR*0.0000001*4*pi;		
    	}
//    	tri_up[mesh-1]=0;
        currentInver[mesh-1]=(igridr+delR*mesh)*delR*delR*0.0000001*4*pi;
    	float c=RTOZ;
    	float cos_temp;
    	float constant_main;
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
    	cudaMalloc((void**)&cu_act_trace,((mesh*mesh)*sizeof(float)));
    	cudaMemcpy(cu_act_trace,act_trace,(mesh*mesh)*sizeof(float),cudaMemcpyHostToDevice);
    	cudaMalloc((void**)&cu_tri_main,((mesh*mesh)*sizeof(float)));
    	cudaMemcpy(cu_tri_main,tri_main,(mesh*mesh)*sizeof(float),cudaMemcpyHostToDevice);
   	float ele_bfft[(mesh+1)*(mesh+1)];
    	float ele[mesh*mesh];
    	float fft_mat[(mesh+1)*(mesh+1)];
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
	float radix[mesh];
	for(i=0;i<mesh;i++)
	{
		radix[i]=igridr+delR+delR*i;

	}
    	cudaMemcpyToSymbol(cu_radix,radix,(mesh)*sizeof(float));
    	cudaMemcpyToSymbol(cu_currentInver,currentInver,mesh*sizeof(float));
    	cudaMalloc((void**)&cu_mat, (mesh+1)*(mesh+1)*sizeof(float));
    	cudaMemcpy(cu_mat,fft_mat,(mesh+1)*(mesh+1)*sizeof(float) , cudaMemcpyHostToDevice);

	float tani[numbdry];
	float coti[numbdry];
	for (i=0;i<numbdry;i++)
	{
		tani[i] = tan(i*2*pi/numbdry);
		coti[i] = 1/tani[i];
//		printf("% 10.9E % 10.9E\n",tani[i],coti[i]);
	}
	cudaMalloc((void**)&cu_tani, numbdry*sizeof(float));
        cudaMemcpy(cu_tani,tani,numbdry*sizeof(float) , cudaMemcpyHostToDevice);
	cudaMalloc((void**)&cu_coti, numbdry*sizeof(float));
        cudaMemcpy(cu_coti,coti,numbdry*sizeof(float) , cudaMemcpyHostToDevice);
}

void dataprocess(float* a, float* b, int ri, int ci, int ro )
{
	for(int i=0;i<ro*ci;i++)
		b[i] = 0;
	for(int i=0;i<ci;i++)
	{
		for(int j=0;j<ri;j++)
		{
			b[i*ro+j] = a[i*ri+j];
		}
	}
}

void boundData(float* a, float* b)
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

int arrange_init()
{
	int num[mesh*mesh];
    	FILE *Number;
	Number=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/arrange.dat","r");
	if(Number==0)
	{
		printf("error: can not open arrange.dat\n");
		return 1;
	}
	if(!fread(num,sizeof(int),mesh*mesh,Number))
	{
		printf("File 1 is empty\n");
		return 1;
	}
	fclose(Number);
	int num_se[ELI_NUM];
	int num_block[mesh];
	int blo_accu[mesh];
	num_regu(num,num_se,num_block,blo_accu);
    	cudaMemcpyToSymbol(cu_num_block, num_block,mesh*sizeof(int));
    	cudaMemcpyToSymbol(cu_blo_accu,blo_accu,mesh*sizeof(int));
//    	cudaMemcpyToSymbol(cu_num_se,num_se,ELI_NUM*sizeof(int));
	cudaMalloc((void**)&cu_num_se,ELI_NUM*sizeof(int));
    	cudaMemcpy(cu_num_se,num_se,ELI_NUM*sizeof(int),cudaMemcpyHostToDevice);

    	return 0;
}

int read_gfunction_files()
{
	int length = EF_NUM_GRIDPTS_FULLSET_INVESSEL * (EF_NSILOP+EF_MAGPRI+EF_NROGOW);
        float * efDiagresponse;
	efDiagresponse= (float *)malloc(ELI_NUM*EF_DIAG*sizeof(float));
	float * refDiagresponse;
	refDiagresponse= (float *)malloc(ELI_NUM*EF_DIAG*sizeof(float));
	FILE *Number2;
	Number2=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/ef_diagresponse.dat","r");
	if(Number2==0)
	{
		printf("error: can not open ef_diagresponse.dat\n");
		return 1;
	}
	if(!fread(efDiagresponse,sizeof(float),length,Number2))
	{
		printf("File 2 is empty\n");
		return 1;
	}	
	fclose(Number2);
	dataprocess(efDiagresponse,refDiagresponse,ELI_NUM,EF_DIAG,ELI_NUM);
    	cudaMalloc((void**)&cu_efDiagresponse,(ELI_NUM*EF_DIAG*sizeof(float)));
    	cudaMemcpy(cu_efDiagresponse,refDiagresponse,ELI_NUM*EF_DIAG*sizeof(float),cudaMemcpyHostToDevice);
    	free(refDiagresponse);

	int length2 = EF_NFCOIL * (EF_NSILOP+EF_MAGPRI+EF_NROGOW);
//	float gfc1[EF_NFCOIL * (EF_NSILOP+EF_MAGPRI+EF_NROGOW)];
	float gfc[EF_NFCOIL * (EF_NSILOP+EF_MAGPRI+EF_NROGOW)];	
 	FILE* Number3;
	Number3=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/ef_gc.dat","r");
	if(Number3==0)
	{
		printf("error: can not open ef_gc.dat\n");
		return 1;
	}
	if(!fread(gfc,sizeof(float),length2,Number3))
	{
		printf("File 3 is empty\n");
		return 1;
	}	
	fclose(Number3);
    	cudaMalloc((void**)&cu_gfc, ((EF_NFCOIL+(int)current_representation[25])*EF_DIAG)*sizeof(float));
	cudaMemcpy(cu_gfc,gfc,length2*sizeof(float),cudaMemcpyHostToDevice);
	responsematrix = cu_gfc + length2;

	int length9 = EF_ECCOIL * (EF_NSILOP+EF_MAGPRI+EF_NROGOW);
//	float gec1[EF_ECCOIL * (EF_NSILOP+EF_MAGPRI+EF_NROGOW)];
	float gec[EF_ECCOIL * (EF_NSILOP+EF_MAGPRI+EF_NROGOW)];	
 	FILE* Number9;
	Number9=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/ef_ec_gc.dat","r");
	if(Number9==0)
	{
		printf("error: can not open ef_ec_gc.dat\n");
		return 1;
	}
	if(!fread(gec,sizeof(float),length9,Number9))
	{
		printf("File 9 is empty\n");
		return 1;
	}	
	fclose(Number9);

	cudaMalloc((void**)&cu_ec_g, (EF_ECCOIL*EF_DIAG)*sizeof(float));
	cudaMemcpy(cu_ec_g,gec,length9*sizeof(float),cudaMemcpyHostToDevice);

	int length3 = mesh*4*ELI_NUM;
	float* hgridresponse;
	hgridresponse= (float *)malloc(length3*sizeof(float));
 	FILE* Number4;
	Number4=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/ef_gridresponse.dat","r");
	if(Number4==0)
	{
		printf("error: can not open ef_gridresponse.dat\n");
		return 1;
	}
	if(!fread(hgridresponse,sizeof(float),length3,Number4))
	{
		printf("File 4 is empty\n");
		return 1;
	}	
	fclose(Number4);

	float* boundresponse;
	boundresponse= (float *)malloc(ELI_NUM*(mesh+1)*4*sizeof(float));
	boundData(hgridresponse,boundresponse);
    	cudaMalloc((void**)&cu_boundresponse,(ELI_NUM*(mesh+1)*4*sizeof(float)));
    	cudaMemcpy(cu_boundresponse,boundresponse,ELI_NUM*(mesh+1)*4*sizeof(float),cudaMemcpyHostToDevice);
	free(hgridresponse);
	free(boundresponse);

	float* hicresponse;
	hicresponse = (float *)malloc(EF_NFCOIL*nwnh*sizeof(float));
	FILE* PFCR;
	PFCR=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/ef_pfresponse.dat","r");
	if(PFCR==0)
	{
		printf("error: can not open ef_pfresponse.dat\n");
		return 1;
	}
	if(!fread(hicresponse,sizeof(float),EF_NFCOIL*nwnh,PFCR))
	{
		printf("File 4 is empty\n");
		return 1;
	}	
	fclose(PFCR);
	cudaMalloc((void**)&cu_icresponse,(EF_NFCOIL*nwnh*sizeof(float)));
    	cudaMemcpy(cu_icresponse,hicresponse,EF_NFCOIL*nwnh*sizeof(float),cudaMemcpyHostToDevice);
	free(hicresponse);

	float* ecgridresponse;
	ecgridresponse= (float *)malloc(EF_ECCOIL*nwnh*sizeof(float));
 	FILE* ECCR;
	ECCR=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/ef_ecresponse.dat","r");
	if(ECCR==0)
	{
		printf("error: can not open ef_ecresponse.dat\n");
		return 1;
	}
	if(!fread(ecgridresponse,sizeof(float),EF_ECCOIL*nwnh,ECCR))
	{
		printf("File 8 is empty\n");
		return 1;
	}	
	fclose(ECCR);
	cudaMalloc((void**)&cu_ecgridresponse,(nwnh*EF_ECCOIL*sizeof(float)));
    	cudaMemcpy(cu_ecgridresponse,ecgridresponse,nwnh*EF_ECCOIL*sizeof(float),cudaMemcpyHostToDevice);
	free(ecgridresponse);

/////////green table for mse///////////////
	num_channels_mse = (float *)malloc(1*sizeof(float));
	if (MSE_FIT == 1)
	{
	        float * grid_mse_br;
		grid_mse_br = (float *)malloc(ELI_NUM*17*mesh2*sizeof(float));
        	float * grid_mse_bz;
		grid_mse_bz = (float *)malloc(ELI_NUM*17*mesh2*sizeof(float));

		FILE *mse_br_grid;
		mse_br_grid=fopen("mse_br_grid.dat","r");
		if(mse_br_grid==0)
		{
			printf("error: can not open mse_br_grid.dat\n");
			return 1;
		}
		if(!fread(grid_mse_br,sizeof(float),ELI_NUM*17*mesh2,mse_br_grid))
		{
			printf("mse_br_grid is empty\n");
			return 1;
		}	
		fclose(mse_br_grid);

		FILE *mse_bz_grid;
		mse_bz_grid=fopen("mse_bz_grid.dat","r");
		if(mse_bz_grid==0)
		{
			printf("error: can not open mse_bz_grid.dat\n");
			return 1;
		}
		if(!fread(grid_mse_bz,sizeof(float),ELI_NUM*17*mesh2,mse_bz_grid))
		{
			printf("mse_bz_grid is empty\n");
			return 1;
		}	
		fclose(mse_bz_grid);

	        float * pf_mse_br;
		pf_mse_br = (float *)malloc(EF_NFCOIL*17*mesh2*sizeof(float));
        	float * pf_mse_bz;
		pf_mse_bz = (float *)malloc(EF_NFCOIL*17*mesh2*sizeof(float));

		FILE *mse_br_pf;
		mse_br_pf=fopen("mse_br_pf.dat","r");
		if(mse_br_pf==0)
		{
			printf("error: can not open mse_br_pf.dat\n");
			return 1;
		}
		if(!fread(pf_mse_br,sizeof(float),EF_NFCOIL*17*mesh2,mse_br_pf))
		{
			printf("mse_br_pf is empty\n");
			return 1;
		}	
		fclose(mse_br_pf);

		FILE *mse_bz_pf;
		mse_bz_pf=fopen("mse_bz_pf.dat","r");
		if(mse_bz_pf==0)
		{
			printf("error: can not open mse_bz_pf.dat\n");
			return 1;
		}
		if(!fread(pf_mse_bz,sizeof(float),EF_NFCOIL*17*mesh2,mse_bz_pf))
		{
			printf("mse_bz_pf is empty\n");
			return 1;
		}	
		fclose(mse_bz_pf);

		float * ec_mse_br;
		ec_mse_br = (float *)malloc(EF_ECCOIL*17*mesh2*sizeof(float));
        	float * ec_mse_bz;
		ec_mse_bz = (float *)malloc(EF_ECCOIL*17*mesh2*sizeof(float));

		FILE *mse_br_ec;
		mse_br_ec=fopen("mse_br_ec.dat","r");
		if(mse_br_ec==0)
		{
			printf("error: can not open mse_br_ec.dat\n");
			return 1;
		}
		if(!fread(ec_mse_br,sizeof(float),EF_ECCOIL*17*mesh2,mse_br_ec))
		{
			printf("mse_br_ec is empty\n");
			return 1;
		}	
		fclose(mse_br_ec);

		FILE *mse_bz_ec;
		mse_bz_ec=fopen("mse_bz_ec.dat","r");
		if(mse_bz_ec==0)
		{
			printf("error: can not open mse_bz_ec.dat\n");
			return 1;
		}
		if(!fread(ec_mse_bz,sizeof(float),EF_ECCOIL*17*mesh2,mse_bz_ec))
		{
			printf("mse_bz_ec is empty\n");
			return 1;
		}	
		fclose(mse_bz_ec);

		FILE* Number11;
		Number11=fopen("num_channels_mse.dat","rt");
		if(Number11 == 0)
		{
			printf("error: can not open num_channels_mes.dat\n");
			return 1;
		}
		for(int i=0; i<1; i++)
		{
			fscanf(Number11,"%f",num_channels_mse+i);
		}
		fclose(Number11);
		
		float * mse_rz;
		mse_rz = (float *)malloc(num_channels_mse[0]*2*sizeof(float));		
		FILE* Number12;
		Number12=fopen("mse_rz.dat","rt");
		if(Number12 == 0)
		{
			printf("error: can not open mse_rz.dat\n");
			return 1;
		}
		for(int i=0; i<num_channels_mse[0]*2; i++)
		{
			fscanf(Number12,"%f",mse_rz+i);
		}
		fclose(Number12);

		float r,z,r0,z0,a,b;
		float *factor;
		int index01,index02,index03,index04,rr,zz;
		factor = (float *)malloc(num_channels_mse[0]*4*sizeof(float));
		grid_mse_bp = (float *)malloc(num_channels_mse[0]*ELI_NUM*sizeof(float));
		pf_mse_bp = (float *)malloc(num_channels_mse[0]*EF_NFCOIL*sizeof(float));
		ec_mse_bp = (float *)malloc(num_channels_mse[0]*EF_ECCOIL*sizeof(float));
		
		for (int i=0;i<num_channels_mse[0];i++)
		{
			r = mse_rz[i*2];
			r0 = int((r-igridr)/delR)*delR+igridr;
			rr = int((r-igridr)/delR);

			z = mse_rz[i*2+1];
			z0 = int((z-(-0.1))/delZ)*delZ+(-0.1);
			zz = int((z-(-0.1))/delZ);
			
			factor[i*4] = (r0+delR-r)*(z0+delZ-z)/darea;
			factor[i*4+1] = (r-r0)*(z0+delZ-z)/darea;
			factor[i*4+2] = (r0+delR-r)*(z-z0)/darea;
			factor[i*4+3] = (r-r0)*(z-z0)/darea;

			index01 = zz*mesh2 + rr;
			index02 = zz*mesh2 + rr + 1;
			index03 = (zz + 1)*mesh2 + rr;
			index04 = (zz + 1)*mesh2 + rr + 1;

			for (int j=0;j<ELI_NUM;j++)
			{
				grid_mse_bp[i*ELI_NUM+j] = a*(grid_mse_br[index01*ELI_NUM+j]*factor[i*4]+grid_mse_br[index02*ELI_NUM+j]*factor[i*4+1]+grid_mse_br[index03*ELI_NUM+j]*factor[i*4+2]+grid_mse_br[index04*ELI_NUM+j]*factor[i*4+3]) + b*(grid_mse_bz[index01*ELI_NUM+j]*factor[i*4]+grid_mse_bz[index02*ELI_NUM+j]*factor[i*4+1]+grid_mse_bz[index03*ELI_NUM+j]*factor[i*4+2]+grid_mse_bz[index04*ELI_NUM+j]*factor[i*4+3]);	
			}
			for (int j=0;j<EF_NFCOIL;j++)
			{
				pf_mse_bp[i*ELI_NUM+j] = a*(pf_mse_br[index01*ELI_NUM+j]*factor[i*4]+pf_mse_br[index02*ELI_NUM+j]*factor[i*4+1]+pf_mse_br[index03*ELI_NUM+j]*factor[i*4+2]+pf_mse_br[index04*ELI_NUM+j]*factor[i*4+3]) + b*(pf_mse_bz[index01*ELI_NUM+j]*factor[i*4]+pf_mse_bz[index02*ELI_NUM+j]*factor[i*4+1]+pf_mse_bz[index03*ELI_NUM+j]*factor[i*4+2]+pf_mse_bz[index04*ELI_NUM+j]*factor[i*4+3]);
			}
			for (int j=0;j<EF_ECCOIL;j++)
			{
				ec_mse_bp[i*ELI_NUM+j] = a*(ec_mse_br[index01*ELI_NUM+j]*factor[i*4]+ec_mse_br[index02*ELI_NUM+j]*factor[i*4+1]+ec_mse_br[index03*ELI_NUM+j]*factor[i*4+2]+ec_mse_br[index04*ELI_NUM+j]*factor[i*4+3]) + b*(ec_mse_bz[index01*ELI_NUM+j]*factor[i*4]+ec_mse_bz[index02*ELI_NUM+j]*factor[i*4+1]+ec_mse_bz[index03*ELI_NUM+j]*factor[i*4+2]+ec_mse_bz[index04*ELI_NUM+j]*factor[i*4+3]);
			}
		}	
	}

///////////////////////////////////////////
	
	return 0;
}

int read_setup_files()
{
	float * normFlux;
	normFlux= (float *)malloc((mesh+2)*(mesh+2)*sizeof(float));
 	FILE* Number5;
	Number5=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/psi_norm.dat","r");
	if(Number5==0)
	{
		printf("error: can not open xpsi.dat\n");
		return 1;
	}
	if(!fread(normFlux,sizeof(float),(mesh+2)*(mesh+2),Number5))
	{
		printf("File normFlux is empty\n");
		return 1;
	}	
	fclose(Number5);

	float * initNormFlux;
	initNormFlux= (float *)malloc(mesh*mesh*sizeof(float));
	for(int i=0;i<mesh;i++)
		for(int j=0;j<mesh;j++)
			initNormFlux[i*mesh+j] = normFlux[(i+1)*(mesh+2)+j+1];

    	cudaMalloc((void**)&cu_initNormFlux,(mesh*mesh*sizeof(float)));
   	cudaMemcpy(cu_initNormFlux,initNormFlux,mesh*mesh*sizeof(float),cudaMemcpyHostToDevice);
	free(initNormFlux);
	free(normFlux);

	float limiter[LIMITER_NUM*2];
    	FILE *Number8;
	Number8=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/limiter.dat","r");
	if(Number8==0)
	{
		printf("error: can not open limiter.dat\n");
		return 1;
	}
	for(int i=0; i< LIMITER_NUM*2; i++)
	{
		fscanf(Number8,"%f",limiter+i);
	}
	fclose(Number8);
	cudaMalloc((void**)&cu_limiter,(LIMITER_NUM*2*sizeof(float)));
    	cudaMemcpy(cu_limiter,limiter,LIMITER_NUM*2*sizeof(float),cudaMemcpyHostToDevice);

	return 0;
}

int test_diag_read(int argc, char* argv[])
{
/**********************************************************
 *       read the diagnostic data                       *
 *********************************************************/
	slice = (float *)malloc(1*sizeof(float));
	FILE* Number11;
	Number11=fopen("sliceno.dat","rt");
	if(Number11 == 0)
	{
		printf("error: can not open sliceno.dat\n");
		return 1;
	}
	for(int i=0; i<1; i++)
	{
		fscanf(Number11,"%f",slice+i);
	}
	fclose(Number11);
	cudaMalloc((void**)&cu_slice,(1*sizeof(float)));
    	cudaMemcpy(cu_slice,slice,1*sizeof(float),cudaMemcpyHostToDevice);
	
	islice = int(slice[0]);

//	float * ishot;
	ishot = (float *)malloc(islice*sizeof(float));
	
	Number11=fopen("shotno.dat","rt");
	if(Number11 == 0)
	{
		printf("error: can not open shotno.dat\n");
		return 1;
	}
	for(int i=0; i<islice; i++)
	{
		fscanf(Number11,"%f",ishot+i);
	}
	fclose(Number11);
	cudaMalloc((void**)&cu_ishot,(islice*sizeof(float)));
    	cudaMemcpy(cu_ishot,ishot,islice*sizeof(float),cudaMemcpyHostToDevice);

//	float * itime;
	itime = (float *)malloc(islice*sizeof(float));
	
	Number11=fopen("time.dat","rt");
	if(Number11 == 0)
	{
		printf("error: can not open time.dat\n");
		return 1;
	}
	for(int i=0; i<islice; i++)
	{
		fscanf(Number11,"%f",itime+i);
	}
	fclose(Number11);
	cudaMalloc((void**)&cu_itime,(islice*sizeof(float)));
    	cudaMemcpy(cu_itime,itime,islice*sizeof(float),cudaMemcpyHostToDevice);

	float * psiref;
	psiref = (float *)malloc(islice*sizeof(float));
	
	Number11=fopen("siref.dat","rt");
	if(Number11 == 0)
	{
		printf("error: can not open siref.dat\n");
		return 1;
	}
	for(int i=0; i<islice; i++)
	{
		fscanf(Number11,"%f",psiref+i);
	}
	fclose(Number11);
	cudaMalloc((void**)&cu_psiref,(islice*sizeof(float)));
    	cudaMemcpy(cu_psiref,psiref,islice*sizeof(float),cudaMemcpyHostToDevice);


	btor = (float *)malloc(islice*sizeof(float));
	
	Number11=fopen("btor.dat","rt");
	if(Number11 == 0)
	{
		printf("error: can not open siref.dat\n");
		return 1;
	}
	for(int i=0; i<islice; i++)
	{
		fscanf(Number11,"%f",btor+i);
	}
	fclose(Number11);
	cudaMalloc((void**)&cu_btor,(islice*sizeof(float)));
    	cudaMemcpy(cu_btor,btor,islice*sizeof(float),cudaMemcpyHostToDevice);


//	float * diagnostic1;
	diagnostic1 = (float *)malloc(islice*EF_DIAG*sizeof(float));
 	
	Number11=fopen("magnetic.dat","rt");
	if(Number11==0)
	{
		printf("error: can not open magnetic.dat\n");
		return 1;
	}
	for(int i = 0;i<islice*EF_DIAG;i++)
	{
		fscanf(Number11,"%f",diagnostic1+i);
	}
	fclose(Number11);
	cudaMalloc((void**)&cu_diagnostic1,(islice*EF_DIAG*sizeof(float)));
    	cudaMemcpy(cu_diagnostic1,diagnostic1,islice*EF_DIAG*sizeof(float),cudaMemcpyHostToDevice);

//	test_diagnostic1 = cu_diagnostic1;
	if (current_representation[0] == 1) 
	{ 
		cudaMalloc((void**)&cu_pro_diag,((EF_DIAG+(int)current_representation[26])*sizeof(float)));
	}
	if (current_representation[0] == 2)
	{
		cudaMalloc((void**)&cu_pro_diag,EF_DIAG*sizeof(float));
	}
//	if (current_representation[0] == 3)
//	{
//		cudaMalloc((void**)&cu_pro_diag_fb,numbdry_gafile*sizeof(double));
//	}
	
	float * diagbit1;
	diagbit1 = (float *)malloc(islice*EF_DIAG*sizeof(float));
 	FILE* Number111;
	Number111=fopen("diagbit.dat","rt");
	if(Number111==0)
	{
		printf("error: can not open diagbit1.dat\n");
		return 1;
	}
	for(int i = 0;i<islice*EF_DIAG;i++)
	{
		fscanf(Number111,"%f",diagbit1+i);
	}
	fclose(Number111);
	cudaMalloc((void**)&cu_diagbit1,(islice*EF_DIAG*sizeof(float)));
    	cudaMemcpy(cu_diagbit1,diagbit1,islice*EF_DIAG*sizeof(float),cudaMemcpyHostToDevice);

	diagnostic2 = (float *)malloc(islice*EF_ECCOIL*sizeof(float));
 	FILE* Number12;
	Number12=fopen("ecurrent.dat","rt");
	if(Number12==0)
	{
		printf("error: can not open ecurrent.dat\n");
		return 1;
	}
	for(int i = 0;i<islice*EF_ECCOIL;i++)
	{
		fscanf(Number12,"%f",diagnostic2+i);
	}
	fclose(Number12);
	cudaMalloc((void**)&cu_diagnostic2,(islice*EF_ECCOIL*sizeof(float)));
    	cudaMemcpy(cu_diagnostic2,diagnostic2,islice*EF_ECCOIL*sizeof(float),cudaMemcpyHostToDevice);	

//	cudaMalloc((void**)&cu_pro_diag2,(EF_DIAG*sizeof(float)));
//	test_diagnostic2 = cu_diagnostic2;

	float weight[EF_DIAG*islice];
    	FILE *Number10;
	Number10=fopen("weight.dat","r");
//	Number10=fopen("weight.dat","r");
	if(Number10==0)
	{
		printf("error: can not open weight.dat\n");
		return 1;
	}
	for(int i = 0;i<islice*EF_DIAG;i++)
		fscanf(Number10,"%f",weight+i);
	fclose(Number10);
	cudaMalloc((void**)&cu_weight,(islice*EF_DIAG*sizeof(float)));
    	cudaMemcpy(cu_weight,weight,islice*EF_DIAG*sizeof(float),cudaMemcpyHostToDevice);

	result_xylimtr = (float *)malloc(2*numlim*sizeof(float));
	FILE *Number19;
	Number19=fopen("/u/huangyao/P-EFIT/New_Folder/pefit_257_jt_new/bin/xylim.dat","r");
	if(Number19==0)
	{
		printf("error: can not open xylim.dat\n");
		return 1;
	}
	for(int i=0; i<2*numlim; i++)
		fscanf(Number19,"%f", result_xylimtr+i);
	fclose(Number19);

	if (MSE_FIT == 1)
	{
	}
	if (current_representation[0] == 2)
	{
		mpress = (int *)malloc(1*sizeof(int));
		if (PRES_FIT == 1)
		{
			Number11=fopen("pressure_input.dat","rt");
			if(Number11==0)
			{
				printf("error: can not open pressure.dat\n");
				return 1;
			}

			fscanf(Number11,"%d",mpress);
		
			float *pressure;
			pressure = (float *)malloc(mpress[0]*sizeof(float));
			float *pres_psi;
			pres_psi = (float *)malloc(mpress[0]*sizeof(float));
			float *pres_uncer;
			pres_uncer = (float *)malloc(mpress[0]*sizeof(float));

			for(int i = 0;i<mpress[0];i++)
			{
				fscanf(Number11,"%f",pres_psi+i);
				printf("pres %d = %f \n",i+1,pres_psi[i]);
			}
			for(int i = 0;i<mpress[0];i++)
			{
				fscanf(Number11,"%f",pressure+i);
				printf("pressure %d = %f \n",i+1,pressure[i]);
			}
			for(int i = 0;i<mpress[0];i++)
			{
				fscanf(Number11,"%f",pres_uncer+i);
				printf("pres_uncer %d = %f \n",i+1,pres_uncer[i]);
			}

			fclose(Number11);

			cudaMalloc((void**)&cu_press,mpress[0]*sizeof(float));
    			cudaMemcpy(cu_press,pressure,mpress[0]*sizeof(float),cudaMemcpyHostToDevice);

			cudaMalloc((void**)&cu_pres_psi,mpress[0]*sizeof(float));
    			cudaMemcpy(cu_pres_psi,pres_psi,mpress[0]*sizeof(float),cudaMemcpyHostToDevice);
			
			cudaMalloc((void**)&cu_press_uncer,mpress[0]*sizeof(float));
    			cudaMemcpy(cu_press_uncer,pres_uncer,mpress[0]*sizeof(float),cudaMemcpyHostToDevice);
		}
		mjphi = (int *)malloc(1*sizeof(int));
		if (JPHI_FIT == 1)
		{
			Number11=fopen("jphi_input.dat","rt");
			if(Number11==0)
			{
				printf("error: can not open jphi.dat\n");
				return 1;
			}

			fscanf(Number11,"%d",mjphi);
		
			float *jphi_npsi;
			jphi_npsi = (float *)malloc(mjphi[0]*sizeof(float));
			float *j_phi;
			j_phi = (float *)malloc(mjphi[0]*sizeof(float));
			float *j_phi_uncer;
			j_phi_uncer = (float *)malloc(mjphi[0]*sizeof(float));

			for(int i = 0;i<mjphi[0];i++)
			{
				fscanf(Number11,"%f",jphi_npsi+i);
				printf("jphi_npsi %d = %f \n",i+1,jphi_npsi[i]);
			}
			for(int i = 0;i<mjphi[0];i++)
			{
				fscanf(Number11,"%f",j_phi+i);
				printf("j_phi %d = %f \n",i+1,j_phi[i]);
			}
			for(int i = 0;i<mjphi[0];i++)
			{
				fscanf(Number11,"%f",j_phi_uncer+i);
				printf("j_phi_uncer %d = %f \n",i+1,j_phi_uncer[i]);
			}

			fclose(Number11);

			cudaMalloc((void**)&cu_num_jphi,1*sizeof(int));
    			cudaMemcpy(cu_num_jphi,mjphi,1*sizeof(int),cudaMemcpyHostToDevice);

			cudaMalloc((void**)&cu_jphi_npsi,mjphi[0]*sizeof(float));
    			cudaMemcpy(cu_jphi_npsi,jphi_npsi,mjphi[0]*sizeof(float),cudaMemcpyHostToDevice);

			cudaMalloc((void**)&cu_j_phi,mjphi[0]*sizeof(float));
    			cudaMemcpy(cu_j_phi,j_phi,mjphi[0]*sizeof(float),cudaMemcpyHostToDevice);

			cudaMalloc((void**)&cu_j_phi_uncer,mjphi[0]*sizeof(float));
    			cudaMemcpy(cu_j_phi_uncer,j_phi_uncer,mjphi[0]*sizeof(float),cudaMemcpyHostToDevice);
		}
	}

	return 0;
}

int allocate_inter_vars()
{
	if (current_representation[0] == 1 || current_representation[0] == 2)
	{	
		int* offset_da = (int *)malloc(1*sizeof(int));

		diag_num_all = (int *)malloc(1*sizeof(int));

		cudaMalloc((void**)&cu_diag_num_all, 1*sizeof(int));
		
		cudaMalloc((void**)&temp_cu_result, (mesh+1)*(mesh+1)*sizeof(float));

    		cudaMalloc((void**)&temp_dmate, 20000*sizeof(float));

   	 	cudaMalloc((void**)&cu_bfft, (mesh+1)*(mesh+1)*sizeof(float));

  	  	cudaMalloc((void**)&cu_result, 500*(mesh+1)*(mesh+1)*sizeof(float));

		cudaMalloc((void**)&cu_result_compare, ((mesh+1)*(mesh+1)+5000)*sizeof(float));

 	   	cudaMalloc((void**)&xpoint, 12*500*sizeof(float));

 	   	cudaMalloc((void**)&leas_flux, (mesh*EF_DIAG*(EF_NFCOIL+(int)current_representation[25])*10)*sizeof(float));

  	  	cudaMalloc((void**)&fwt, (EF_DIAG)*sizeof(float));

 	   	cudaMalloc((void**)&cubuffer, (40000)*sizeof(float));

  	  	cudaMalloc((void**)&cuLbuffer, (1500000)*sizeof(float));

		if (current_representation[0] == 1)
		{
			cudaMalloc((void**)&cu_dgfc, (((int)current_representation[25]+EF_NFCOIL)*(EF_DIAG+(int)current_representation[26]))*sizeof(float));

			/////////////////////////// Mapped memory for SVD////////////////////////////

			cudaHostAlloc((void**)&h_A,((int)current_representation[25]+EF_NFCOIL)*((int)current_representation[26]+EF_DIAG)*sizeof(float),cudaHostAllocMapped);

			cudaHostAlloc((void**)&h_S,((int)current_representation[25]+EF_NFCOIL)*sizeof(float),cudaHostAllocMapped);

			cudaHostAlloc((void**)&h_U,((int)current_representation[26]+EF_DIAG)*((int)current_representation[26]+EF_DIAG)*sizeof(float),cudaHostAllocMapped);
	
			cudaHostAlloc((void**)&h_VT,((int)current_representation[25]+EF_NFCOIL)*((int)current_representation[25]+EF_NFCOIL)*sizeof(float),cudaHostAllocMapped);

			cudaHostGetDevicePointer((void**)&d_A,(void*)h_A,0);

			cudaHostGetDevicePointer((void**)&d_S,(void*)h_S,0);

			cudaHostGetDevicePointer((void**)&d_U,(void*)h_U,0);

			cudaHostGetDevicePointer((void**)&d_VT,(void*)h_VT,0);

			//////////////////////////////Device memory for SVD/////////////////////////////
	
			cudaMalloc((void**)&cu_U_matrix,((int)current_representation[26]+EF_DIAG)*((int)current_representation[26]+EF_DIAG)*sizeof(float));

   	 		cudaMalloc((void**)&cu_S_vector,((int)current_representation[25]+EF_NFCOIL)*sizeof(float));

  	 	 	cudaMalloc((void**)&cu_VT_matrix,((int)current_representation[25]+EF_NFCOIL)*((int)current_representation[25]+EF_NFCOIL)*sizeof(float));

			cudaMalloc((void**)&cu_response_inverse_temp,((int)current_representation[26]+EF_DIAG)*sizeof(float));
	
			////////////////////////////////////////////////////////////
	
		}

		if (current_representation[0] == 2)
		{
			/////////////////////////// Mapped memory for SVD////////////////////////////

			printf("allocate_inter_vars 01 is ok\n");

			diag_num_all[0] = EF_DIAG + MSE_FIT*num_channels_mse[0] + PRES_FIT*mpress[0] + JPHI_FIT*mjphi[0] + (num_knots[0]-2)*6;

			printf("diag_num_all = %d \n",diag_num_all[0]);

			cudaMemcpy(cu_diag_num_all,diag_num_all,1*sizeof(int),cudaMemcpyHostToDevice);

			cudaMalloc((void**)&cu_dgfc, (((int)current_representation[25]+EF_NFCOIL)*EF_DIAG)*sizeof(float));

			cudaMalloc((void**)&h_A_left, (((int)current_representation[25]+EF_NFCOIL)*diag_num_all[0])*sizeof(float));

			cudaHostAlloc((void**)&h_A,((int)current_representation[25]+EF_NFCOIL)*diag_num_all[0]*sizeof(float),cudaHostAllocDefault);

			cudaHostAlloc((void**)&h_S,((int)current_representation[25]+EF_NFCOIL)*sizeof(float),cudaHostAllocMapped);

			cudaHostAlloc((void**)&h_U,diag_num_all[0]*diag_num_all[0]*sizeof(float),cudaHostAllocMapped);
	
			cudaHostAlloc((void**)&h_VT,((int)current_representation[25]+EF_NFCOIL)*((int)current_representation[25]+EF_NFCOIL)*sizeof(float),cudaHostAllocMapped);

			cudaHostGetDevicePointer((void**)&d_A,(void*)h_A,0);

			cudaHostGetDevicePointer((void**)&d_S,(void*)h_S,0);

			cudaHostGetDevicePointer((void**)&d_U,(void*)h_U,0);

			cudaHostGetDevicePointer((void**)&d_VT,(void*)h_VT,0);

			//////////////////////////////Device memory for SVD/////////////////////////////
	
			cudaMalloc((void**)&cu_U_matrix,diag_num_all[0]*diag_num_all[0]*sizeof(float));

   	 		cudaMalloc((void**)&cu_S_vector,((int)current_representation[25]+EF_NFCOIL)*sizeof(float));

  	 	 	cudaMalloc((void**)&cu_VT_matrix,((int)current_representation[25]+EF_NFCOIL)*((int)current_representation[25]+EF_NFCOIL)*sizeof(float));

			cudaMalloc((void**)&cu_response_inverse_temp,diag_num_all[0]*sizeof(float));
	
			////////////////////////////////////////////////////////////

			cudaMalloc((void**)&cu_spline_right,(num_knots[0]-2)*6*sizeof(float));

			cudaMalloc((void**)&cu_spline_left,(num_knots[0]-2)*6*((int)current_representation[25]+EF_NFCOIL)*sizeof(float));

			printf("allocate_inter_vars 02 is ok\n");

			if (PRES_FIT == 1)
			{
				cudaMalloc((void**)&cu_pres_right,mpress[0]*sizeof(float));

				cudaMalloc((void**)&cu_pres_left,mpress[0]*((int)current_representation[25]+EF_NFCOIL)*sizeof(float));

				offset_da[0] = EF_DIAG + MSE_FIT*num_channels_mse[0] + JPHI_FIT*mjphi[0];

				cudaMalloc((void**)&offset_dA_pres,1*sizeof(int));

				cudaMemcpy(offset_dA_pres,offset_da,1*sizeof(int),cudaMemcpyHostToDevice);
			}
			if (JPHI_FIT == 1)
			{
				cudaMalloc((void**)&cu_jphi_surface,mjphi[0]*numbdry*2*sizeof(float));

				cudaMalloc((void**)&cu_jphi_right,mjphi[0]*sizeof(float));

				cudaMalloc((void**)&cu_jphi_left,mjphi[0]*((int)current_representation[25]+EF_NFCOIL)*sizeof(float));

				offset_da[0] = EF_DIAG + MSE_FIT*num_channels_mse[0];

				cudaMalloc((void**)&offset_dA_jphi,1*sizeof(int));

				cudaMemcpy(offset_dA_jphi,offset_da,1*sizeof(int),cudaMemcpyHostToDevice);
			}
			if (MSE_FIT == 1)
			{
				offset_da[0] = EF_DIAG;

				cudaMalloc((void**)&offset_dA_mse,1*sizeof(int));

				cudaMemcpy(offset_dA_mse,offset_da,1*sizeof(int),cudaMemcpyHostToDevice);
			}
			offset_da[0] = EF_DIAG + MSE_FIT*num_channels_mse[0] + PRES_FIT*mpress[0] + JPHI_FIT*mjphi[0];

			cudaMalloc((void**)&offset_dA_spline,1*sizeof(int));

			cudaMemcpy(offset_dA_spline,offset_da,1*sizeof(int),cudaMemcpyHostToDevice);

			cudaMalloc((void**)&cu_diag_right,diag_num_all[0]*sizeof(float));
		}

		cudaMalloc((void**)&cu_icflux,nwnh*sizeof(float));

		cudaMalloc((void**)&cu_ecflux,nwnh*sizeof(float));

		cudaMalloc((void**)&test_ec,(EF_DIAG)*sizeof(float));

  	 	cudaMalloc((void**)&test_right,(((int)current_representation[25]+EF_NFCOIL)*sizeof(float)));

 	   	cudaMalloc((void**)&cuBoundResult,((mesh+1)*4*sizeof(float)));

 	   	cudaMalloc((void**)&effect_flux,((int)current_representation[25]+1)*ELI_NUM*sizeof(float));

		cudaMalloc((void**)&cuLimFlux,(LIMITER_NUM*sizeof(float)));

		cudaMalloc((void**)&limZ,(2*sizeof(int)));

		cudaMalloc((void**)&cu_bound_flux,(500*sizeof(float)));

		cudaMalloc((void**)&cu_xlocation,(8*sizeof(float)));
	
		time_values = (float*)malloc((5000*sizeof(float)));

		cudaMalloc((void**)&cu_out_buf,(100*sizeof(float)));

		cudaHostAlloc((void**)&out_buf,100*sizeof(float), 0);

		cudaMalloc((void**)&de_buffer,(100*sizeof(float)));

		cudaMalloc((void**)&xybdry,(numbdry+3)*2*sizeof(float));

		cudaMalloc((void**)&maxgridbdry,numbdry*sizeof(int));

		cudaMalloc((void**)&xysurface,mesh2*numbdry*2*sizeof(float));

		cudaMalloc((void**)&cu_converge, 1*sizeof(int));

		cudaMalloc((void**)&cu_converge_error, 1*sizeof(float));

		cudaHostAlloc((void**)&converge,1*sizeof(int),cudaHostAllocMapped);

		cudaHostGetDevicePointer((void**)&cu_converge,(void*)converge,0);

		result_psi = (float *)malloc(mesh2*mesh2*sizeof(float));

		result_bound = (float *)malloc(8*sizeof(float));

		result_rzxpt = (float *)malloc(8*sizeof(float));

		result_nullpoint = (float *)malloc(12*sizeof(float));

		result_brsp = (float *)malloc(23*sizeof(float));

		result_current = (float *)malloc(1*sizeof(float));

		result_fpol_profile = (float *)malloc(mesh2*sizeof(float));

		result_pres_profile = (float *)malloc(mesh2*sizeof(float));

		result_ffprime_profile = (float *)malloc(mesh2*sizeof(float));

		result_pprime_profile = (float *)malloc(mesh2*sizeof(float));

		result_xybdry = (float *)malloc(2*(numbdry+3)*sizeof(float));

		result_xysurface = (float *)malloc(mesh2*2*numbdry*sizeof(float));

		result_qpsi = (float *)malloc(mesh2*sizeof(float));

		result_parameter = (float *)malloc(10*sizeof(float));

		cudaMalloc((void**)&cu_result_output,(mesh2*mesh2*sizeof(float)));

		cudaMalloc((void**)&cu_profile_result,(4*mesh2*sizeof(float))); 

		cudaMalloc((void**)&cu_pcurrent,(mesh2*mesh2*sizeof(float))); 

		cudaMalloc((void**)&cu_qpsi,(mesh2*sizeof(float)));

		result_pcurrent = (float *)malloc(mesh2*mesh2*sizeof(float));

		cudaMalloc((void**)&cu_compute_diagnotics,(EF_DIAG*sizeof(float))); 

		compute_diagnotics = (float *)malloc(EF_DIAG*sizeof(float));

		converge_error = (float *)malloc(1*sizeof(float));

		ipiv = (int*)malloc((EF_NFCOIL+(int)current_representation[25])*sizeof(float));

		superb = (float*)malloc(((EF_NFCOIL+(int)current_representation[25])-1)*sizeof(float));

		hright = (float*)malloc((EF_NFCOIL+(int)current_representation[25])*sizeof(float));
	}
	
	if (current_representation[0] == 3)
	{
		cudaMalloc((void**)&cu_dgfc_fb, (EF_NFCOIL+(int)current_representation[25])*(EF_DIAG+(int)current_representation[26])*sizeof(double));

		/////////////////////////// Mapped memory for SVD////////////////////////////

		cudaHostAlloc((void**)&h_A_fb,(EF_NFCOIL+(int)current_representation[25])*(EF_DIAG+(int)current_representation[26])*sizeof(double),cudaHostAllocDefault);

		///////////////////////////////////////////////////////////////////////////////////

		cudaHostAlloc((void**)&h_S_fb,(EF_DIAG+(int)current_representation[26])*sizeof(double),cudaHostAllocWriteCombined);

		cudaHostAlloc((void**)&h_U_fb,(EF_DIAG+(int)current_representation[26])*(EF_DIAG+(int)current_representation[26])*sizeof(double),cudaHostAllocDefault);
	
		cudaHostAlloc((void**)&h_VT_fb,(EF_NFCOIL+(int)current_representation[25])*(EF_NFCOIL+(int)current_representation[25])*sizeof(double),cudaHostAllocDefault);

		//////////////////////////////Device memory for SVD/////////////////////////////
	
		cudaMalloc((void**)&cu_U_matrix_fb,(EF_DIAG+(int)current_representation[26])*(EF_DIAG+(int)current_representation[26])*sizeof(double));

   	 	cudaMalloc((void**)&cu_S_vector_fb,(EF_DIAG+(int)current_representation[26])*sizeof(double));

  	  	cudaMalloc((void**)&cu_VT_matrix_fb,(EF_NFCOIL+(int)current_representation[25])*(EF_NFCOIL+(int)current_representation[25])*sizeof(double));

		cudaMalloc((void**)&cu_response_inverse_temp_fb,(EF_DIAG+(int)current_representation[26])*sizeof(double));
	
		////////////////////////////////////////////////////////////

		cudaMalloc((void**)&effect_flux_fb,((int)current_representation[25]+1)*ELI_NUM*sizeof(double));

		cudaMalloc((void**)&test_right_fb,(((int)current_representation[25]+EF_NFCOIL)*sizeof(double)));

		superb_fb = (double*)malloc(((EF_NFCOIL+(int)current_representation[25])-1)*sizeof(double));

		cudaMalloc((void**)&cu_ecflux_fb,nwnh*sizeof(double));

		cudaMalloc((void**)&cu_ipflux_fb,nwnh*sizeof(double));

		cudaMalloc((void**)&cu_bfft_fb, (mesh+1)*(mesh+1)*sizeof(double));

		cudaMalloc((void**)&cuBoundResult_fb,((mesh+1)*4*sizeof(double)));

		cudaMalloc((void**)&cubuffer_fb, (40000)*sizeof(double));

		cudaMalloc((void**)&cuLbuffer_fb, (1500000)*sizeof(double));

		cudaMalloc((void**)&cu_icflux_fb,nwnh*sizeof(double));

		cudaMalloc((void**)&cu_result_fb, 500*(mesh+1)*(mesh+1)*sizeof(double));

		cudaMalloc((void**)&leas_flux_fb, (mesh*EF_DIAG*(EF_NFCOIL+(int)current_representation[25])*10)*sizeof(double));

		cudaMalloc((void**)&xpoint_fb, 12*500*sizeof(double));

		cudaMalloc((void**)&cu_xlocation_fb,(8*sizeof(double)));

		cudaMalloc((void**)&cuLimFlux_fb,(LIMITER_NUM*sizeof(double)));

		cudaMalloc((void**)&cu_bound_flux_fb,(500*sizeof(float)));

		cudaMalloc((void**)&cu_result_compare_fb, ((mesh+1)*(mesh+1)+5000)*sizeof(double));

		cudaMalloc((void**)&cu_converge, 1*sizeof(int));

		cudaHostAlloc((void**)&converge,1*sizeof(int),cudaHostAllocMapped);

		cudaHostGetDevicePointer((void**)&cu_converge,(void*)converge,0);

		cudaMalloc((void**)&cu_converge_error_fb, 1*sizeof(double));

		result_bound_fb = (double *)malloc(8*sizeof(double));

		result_nullpoint_fb = (double *)malloc(12*sizeof(double));

		cudaMalloc((void**)&cu_profile_result_fb,(4*mesh2*sizeof(double))); 

		result_fpol_profile_fb = (double *)malloc(mesh2*sizeof(double));

		result_pres_profile_fb = (double *)malloc(mesh2*sizeof(double));

		result_ffprime_profile_fb = (double *)malloc(mesh2*sizeof(double));

		result_pprime_profile_fb = (double *)malloc(mesh2*sizeof(double));

		cudaMalloc((void**)&cu_result_output_fb,(mesh2*mesh2*sizeof(double)));

		cudaMalloc((void**)&cu_compute_diagnotics_fb,(EF_DIAG*sizeof(double)));

		compute_diagnotics_fb = (double *)malloc(EF_DIAG*sizeof(double));

		cudaMalloc((void**)&maxgridbdry,numbdry*sizeof(int));

		cudaMalloc((void**)&xybdry_fb,(numbdry+3)*2*sizeof(double));

		cudaMalloc((void**)&xysurface_fb,mesh2*numbdry*2*sizeof(double));

		cudaMalloc((void**)&cu_qpsi_fb,(mesh2*sizeof(double)));

		result_qpsi_fb = (double *)malloc(mesh2*sizeof(double));

		result_xybdry_fb = (double *)malloc(2*(numbdry+3)*sizeof(double));

		result_parameter_fb = (double *)malloc(10*sizeof(double));

		converge_error_fb = (double *)malloc(1*sizeof(double));

		hright_fb = (double *)malloc((EF_NFCOIL+(int)current_representation[25])*sizeof(double));

		compute_diagnotics_fb = (double *)malloc(EF_DIAG*sizeof(double));

		result_psi_fb = (double *)malloc(mesh2*mesh2*sizeof(double));

		cudaMalloc((void**)&cu_pro_diag_fb,numbdry_gafile*sizeof(double));

	}

	return 0;

}

__global__ void tridaig_mul(float* effect_flux, float* cu_efDiagresponse, float* leas_flux, float* current_representation)
{
	int k,i,num;
	num = (int)current_representation[25];
	__shared__ float s_flux[32][32];
	__shared__ float s_gdata[32][32];
	unsigned int index = ELI_NUM * threadIdx.y + blockIdx.x * 32 + threadIdx.x;
	unsigned int index2;
	if(threadIdx.y<num)
	{	 
		s_flux[threadIdx.y][threadIdx.x] = effect_flux[index];
	}
	__syncthreads();
	
	float tempsub = 0;
	//following code has been loop unrolled
		
	i=0;		//1st iteration
	index2=(i*32+threadIdx.y)*ELI_NUM+blockIdx.x*32+ threadIdx.x; 
	// the input cu_efDiagresponse is transposed, and the lines up and down read 24*32 elements into shared memory.
	s_gdata[threadIdx.y][threadIdx.x] = cu_efDiagresponse[index2];
	__syncthreads();
	// the "if" section bellow does a sub-matrix multiplication, it will creats a 4*24 sub-matrix.
	if(threadIdx.x<32&&threadIdx.y<num)
	{	
		for(k=0;k<32;k++)			
			tempsub += s_flux[threadIdx.y][k] * s_gdata[threadIdx.x][k];
		leas_flux[blockIdx.x*EF_DIAG*num+threadIdx.y*EF_DIAG+i*32+threadIdx.x] = tempsub;
	}
	__syncthreads();

	i++;        //2nd iteration
	index2=(i*32+threadIdx.y)*ELI_NUM+blockIdx.x*32+ threadIdx.x;
//	if(threadIdx.y==0)		
	s_gdata[threadIdx.y][threadIdx.x] = cu_efDiagresponse[index2];
	__syncthreads();
	if(threadIdx.x<32&&threadIdx.y<num)
	{	
		tempsub =0;
		for(k=0;k<32;k++)			
			tempsub += s_flux[threadIdx.y][k] * s_gdata[threadIdx.x][k];
		leas_flux[blockIdx.x*EF_DIAG*num+threadIdx.y*EF_DIAG+i*32+threadIdx.x] = tempsub;
	}
	__syncthreads();

	i++;		//3rd iteration
	index2=(i*32+threadIdx.y)*ELI_NUM+blockIdx.x*32+ threadIdx.x;
	s_gdata[threadIdx.y][threadIdx.x] = cu_efDiagresponse[index2];
	__syncthreads();
	if(threadIdx.x<32&&threadIdx.y<num)
	{	
		tempsub=0;
		for(k=0;k<32;k++)			
			tempsub += s_flux[threadIdx.y][k] * s_gdata[threadIdx.x][k];
		leas_flux[blockIdx.x*EF_DIAG*num+threadIdx.y*EF_DIAG+i*32+threadIdx.x] = tempsub;
	}
	__syncthreads();

	i++;		//last iteration		
	if(threadIdx.y<25)
	{
		index2=(i*32+threadIdx.y)*ELI_NUM+blockIdx.x*32+ threadIdx.x;
	//	if(threadIdx.y==0)		
		s_gdata[threadIdx.y][threadIdx.x] = cu_efDiagresponse[index2];
	}
	__syncthreads();

	if(threadIdx.x<25&&threadIdx.y<num)
	{	
		tempsub = 0;
		for(k=0;k<32;k++)			
			tempsub += s_flux[threadIdx.y][k] * s_gdata[threadIdx.x][k];
		leas_flux[blockIdx.x*EF_DIAG*num+(threadIdx.y)*EF_DIAG+i*32+threadIdx.x] = tempsub;
//		printf("least flux %f \n",tempsub);
	}			
}

__global__ void mul_postprocess(float* leas_flux, float* responsematrix, float* current_representation)
{
	float k,k2;
	int num;
	num = (int)current_representation[25];
//	__shared__ float temp_data[POLY_NUM][EF_DIAG];
//	__shared__ float temp_buffer[POLY_NUM][EF_DIAG];
//	temp_data[threadIdx.y][threadIdx.x] = leas_flux[threadIdx.y*EF_DIAG+threadIdx.x];

//	k = leas_flux[threadIdx.y*EF_DIAG+threadIdx.x];
	k = leas_flux[blockIdx.x*EF_DIAG+threadIdx.x];

	for(int i=0;i<(ELI_CORE-1);i++)
	{
		k2= leas_flux[(i+1)*EF_DIAG*num + blockIdx.x*EF_DIAG + threadIdx.x];
		k += k2; 
	}
	unsigned int index = blockIdx.x * EF_DIAG +threadIdx.x ;
	__syncthreads();	
	responsematrix[index] = k;
//	printf("least flux %f \n",k);
}

__global__ void	data_process_ec(float* cu_ec_g,float* test_diagnostic2,float* test_ec)
{
	__shared__ float ec[EF_ECCOIL];
	__shared__ float gdata[EF_ECCOIL];

	unsigned int index = blockIdx.x * 6 + threadIdx.x;

	ec[threadIdx.x] = test_diagnostic2[threadIdx.x];
	gdata[threadIdx.x] = cu_ec_g[index];

	ec[threadIdx.x] = ec[threadIdx.x] * gdata[threadIdx.x];

	int k;
//	float temp_sum = 0;

	for(int s=4;s>0;s>>=1)
	{
		k = threadIdx.x + s;
		if(threadIdx.x<s&&k<6)
			ec[threadIdx.x] += ec[k];
	}

	test_ec[blockIdx.x] = ec[0];
}

__global__ void	data_process_right(float* fwt,float* right_side,float* diagbit,float* test_ec,float* cu_weight,float* cu_pro_diag,float* test_psiref)
{
	float k,m,w,t,x,q,ref;
	k = right_side[blockIdx.x];
	w = cu_weight[blockIdx.x];
	t = diagbit[blockIdx.x];
	q = test_ec[blockIdx.x];
	ref = test_psiref[0];

	if (blockIdx.x<1)
	{
		x = 0.001*(k+ref);
		if (fabs(x)<10*t)
		{
			x = 10*t;
		}
	}
	else if (blockIdx.x<EF_NSILOP)
	{
		x = uncertain*k;
		if (fabs(x)<10*t)
		{
			x = 10*t;
		}
	}
	else if (blockIdx.x<(EF_MAGPRI+EF_NSILOP))
	{
		x = uncertain*k;
		if (fabs(x)<10*t)
		{
			x = 10*t;
		}
	}
	else 
	{
		x = uncertain*k;
		if (fabs(x)<10*t)
		{
			x = 10*t;
		}
	}

	k=k-q;
	m=0;
	if(blockIdx.x<EF_NSILOP)
	{
		k = k + ref;	
		if (w!=0)
		{
			m = 1/sqrt(x*x)*w;
		}
		cu_pro_diag[blockIdx.x] = k*m;
		fwt[blockIdx.x] = m;	
//		printf("%f %f %f\n", x, m, k );
	}
	else if(blockIdx.x<(EF_MAGPRI+EF_NSILOP))
	{
		if (w!=0)
		{
			m = 1/sqrt(x*x)*w;
		}
		cu_pro_diag[blockIdx.x] = k*m;
		fwt[blockIdx.x] = m;
//		printf("%f %f %f\n", x, m, k );		
	}
	else 
	{
		if (w!=0)
		{
			m = 1/sqrt(x*x)*w;
		}
		cu_pro_diag[blockIdx.x] = k*m;
		fwt[blockIdx.x] = m;						
	}
//	printf("mm %d %f\n",blockIdx.x,cu_pro_diag[blockIdx.x]);
}

__global__ void	data_process_left(float* fwt,float* cu_gfc,float* cu_dgfc,float* current_representation)
{
	int num = (int)current_representation[26];
	__shared__ float tempA[EF_DIAG];
	tempA[threadIdx.x] = cu_gfc[blockIdx.x*EF_DIAG + threadIdx.x];
//	tempA[1][threadIdx.x] = fwt[threadIdx.x];
	tempA[threadIdx.x] = tempA[threadIdx.x] * fwt[threadIdx.x];
	cu_dgfc[blockIdx.x*(EF_DIAG+num) + threadIdx.x] = tempA[threadIdx.x];
}

__global__ void rightHandSideInit(float* cu_dgfc,float* cu_diag_dat, float* right_side, float* current_representation)
{
	float temp_sum,fsum;
	int n;
	n = (int)current_representation[26];

	__shared__ float temp_left[EF_DIAG+2];
	__shared__ float temp_right[EF_DIAG+2];

	unsigned int index = blockIdx.x * (EF_DIAG+n) + threadIdx.x;
	temp_left[threadIdx.x] = cu_dgfc[index];
	temp_right[threadIdx.x] = cu_diag_dat[threadIdx.x]; 
	temp_left[threadIdx.x] = temp_left[threadIdx.x] * temp_right[threadIdx.x];
	__syncthreads();

	temp_sum = 0;
	fsum = 0;
	if(threadIdx.x<11)
	{
		for(int i=0;i<11;i++)
		{
			temp_sum += temp_left[threadIdx.x*11 + i];
		}
		temp_left[threadIdx.x] = temp_sum;		   
	}
	else if(threadIdx.x==11)
	{
		temp_sum = temp_left[threadIdx.x*11]+temp_left[threadIdx.x*11+1];
		temp_left[threadIdx.x] = temp_sum;
	}
	if(threadIdx.x == 0)
	{
		for(int i=0;i<12;i++)
		{
			fsum += temp_left[i];
		}
		right_side[blockIdx.x] = fsum;
	}	
}	

__global__ void leftHandSideInit(float* cu_dgfc, float* leas_flux, float* current_representation)
{
	int num,n;
	float temp=0,a,b;

	num = EF_NFCOIL + (int)current_representation[25];
	n = (int)current_representation[26];
	
	a = cu_dgfc[blockIdx.x*(EF_DIAG+n)+threadIdx.x];
 	b = cu_dgfc[blockIdx.y*(EF_DIAG+n)+threadIdx.x];
	temp += a*b;
	a = cu_dgfc[blockIdx.x*(EF_DIAG+n)+threadIdx.x+32];
 	b = cu_dgfc[blockIdx.y*(EF_DIAG+n)+threadIdx.x+32];
	temp += a*b;
	a = cu_dgfc[blockIdx.x*(EF_DIAG+n)+threadIdx.x+64];
 	b = cu_dgfc[blockIdx.y*(EF_DIAG+n)+threadIdx.x+64];
	temp += a*b;
	if(threadIdx.x<25+n)
	{
		a = cu_dgfc[blockIdx.x*(EF_DIAG+n)+threadIdx.x+96];
 		b = cu_dgfc[blockIdx.y*(EF_DIAG+n)+threadIdx.x+96];
		temp += a*b;
	}
	leas_flux[(blockIdx.y*num+blockIdx.x)*32+threadIdx.x] = temp;
/*
	sdgfc[threadIdx.x+32][0] = cu_dgfc[blockIdx.x*EF_DIAG+threadIdx.x+32];
	sdgfc[threadIdx.x][1] = cu_dgfc[blockIdx.y*EF_DIAG+threadIdx.x];
	sdgfc[threadIdx.x][1] = cu_dgfc[blockIdx.y*EF_DIAG+threadIdx.x];
*/		
}

__global__ void leftInit(float* temp_dmate, float* dmate)
{	
	__shared__ float sdgfc[32];
	unsigned int index = blockIdx.x*32 + threadIdx.x;
	sdgfc[threadIdx.x] = temp_dmate[index];
	for(unsigned int s=16;s>0;s>>=1)
	{
		int k=threadIdx.x+s;
		if (threadIdx.x<s&&k<32) 
		{
			sdgfc[threadIdx.x] += sdgfc[k];
		}
	}
	if(threadIdx.x == 0)
		dmate[blockIdx.x] = sdgfc[0];
/*
	sdgfc[threadIdx.x+32][0] = cu_dgfc[blockIdx.x*EF_DIAG+threadIdx.x+32];
	sdgfc[threadIdx.x][1] = cu_dgfc[blockIdx.y*EF_DIAG+threadIdx.x];
	sdgfc[threadIdx.x][1] = cu_dgfc[blockIdx.y*EF_DIAG+threadIdx.x];
*/		
}

__global__ void changePos(float* right_side, float* test_right)
{
	test_right[threadIdx.x] = right_side[threadIdx.x];
}

__global__ void changePos_1(float* cu_dgfc, float* d_A, float* current_representation)
{
	int index;
       	index = blockIdx.x * (EF_DIAG+(int)current_representation[26]) + threadIdx.x;
	d_A[index] = cu_dgfc[index];
}

__global__ void changePos_2(float* d_U,float* d_S,float* d_VT,float* current_representation,float* cu_U_matrix,float* cu_S_vector,float* cu_VT_matrix)
{
	int index1,index2,index3;
       	
	index1 = blockIdx.x * (EF_DIAG+(int)current_representation[26]) + threadIdx.x;
	index2 = blockIdx.x * (EF_NFCOIL+(int)current_representation[25]) + threadIdx.x;
	index3 = threadIdx.x;

	cu_U_matrix[index1] = d_U[index1];

	if(blockIdx.x<(EF_NFCOIL+(int)current_representation[25]) && threadIdx.x<(EF_NFCOIL+(int)current_representation[25]))
	{
		cu_VT_matrix[index2] = d_VT[index2];
	}
	else if (blockIdx.x == (EF_NFCOIL+(int)current_representation[25]) && threadIdx.x<(EF_NFCOIL+(int)current_representation[25]))
	{
		cu_S_vector[index3] = d_S[index3];
	}

//	printf("U matrix %d  %f  \n",index1,cu_U_matrix[index1]);
//	printf("changePos_2 done ! \n");
}

__global__ void changePos_3(float* d_U,float* d_S,float* d_VT,float* current_representation,float* cu_U_matrix,float* cu_S_vector,float* cu_VT_matrix,int* diag_num)
{
	int index1,index2,index3;
       	
	index1 = blockIdx.x * diag_num[0] + threadIdx.x;
	index2 = blockIdx.x * (EF_NFCOIL+(int)current_representation[25]) + threadIdx.x;
	index3 = threadIdx.x;

	cu_U_matrix[index1] = d_U[index1];

	if(blockIdx.x<(EF_NFCOIL+(int)current_representation[25]) && threadIdx.x<(EF_NFCOIL+(int)current_representation[25]))
	{
		cu_VT_matrix[index2] = d_VT[index2];
	}
	else if (blockIdx.x == (EF_NFCOIL+(int)current_representation[25]) && threadIdx.x<(EF_NFCOIL+(int)current_representation[25]))
	{
		cu_S_vector[index3] = d_S[index3];
	}

//	printf("U matrix %d  %f  \n",index1,cu_U_matrix[index1]);
//	printf("changePos_2 done ! \n");
}

__global__ void responsematrix_inverse_U(float* m, float* U, float* current_representation, float* result)
{
	int index1, index2, num;
	__shared__ float temp[256];
//	__shared__ float UT[121];
	
//	printf("inverse 1 start!! \n");

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

__global__ void responsematrix_inverse_U_spline(float* m, float* U, float* current_representation, float* result, int* num_diag)
{
	int index1, index2, num;
	__shared__ float temp[256];
//	__shared__ float UT[121];
	
//	printf("inverse 1 start!! \n");

	num = num_diag[0];

	index1 = blockIdx.x * num + threadIdx.x;
	index2 = threadIdx.x;

	temp[index2] = U[index1] * m[index2];
//	if (blockIdx.x == 0)
//		printf("%d U'*m %f U' %f m %f \n",index1,temp[index2],U[index1],m[index2]);
	__syncthreads();

	for(unsigned int s=128;s>0;s>>=1)
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

__global__ void responsematrix_inverse_S(float* result_temp, float* S, float* current_representation)
{
//	int num;
	float s_temp,s_max,s_truncate;
//	__shared__ float temp[60];
//	__shared__ float UT[121];
	
//	printf("inverse 1 start!! \n");

//	num = (EF_NFCOIL+(int)current_representation[25]);

//	index1 = blockIdx.x * num + threadIdx.x;
//	index2 = threadIdx.x;

	s_temp = S[threadIdx.x];
	s_max = S[0];
	s_truncate = s_max * truncate_factor;
//	printf("S %d = %f \n",threadIdx.x,s_temp);
	if (s_temp == 0 || s_temp < s_truncate)
	{
		result_temp[threadIdx.x] = 0;
	}
	else
	{
		result_temp[threadIdx.x] = result_temp[threadIdx.x]/s_temp; 
	}
//		printf("ST*U'*M %d = %f \n",threadIdx.x,result_temp[threadIdx.x]);
}

__global__ void responsematrix_inverse_VT(float* result_temp, float* VT, float* test_right, float* current_representation)
{
	int index1, index2, num;
	__shared__ float temp[64];
//	__shared__ float UT[121];
	
//	printf("inverse 1 start!! \n");

	num = (EF_NFCOIL+(int)current_representation[25]);

	index1 = blockIdx.x * num + threadIdx.x;
	index2 = threadIdx.x;

	temp[index2] = VT[index1] * result_temp[index2];
	__syncthreads();

//	if (blockIdx.x == 0)
	//	printf("%d V* ST* U'*m %f VT %f ST* U'*m %f \n",index1,temp[index2],VT[index1],result_temp[index2]);
//		printf(" %f \n",temp[index2]);

	for(unsigned int s=32;s>0;s>>=1)
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
	}
	if (blockIdx.x ==0 && threadIdx.x == 0)
	{
		for (int i=0;i<num;i++)
		{
			printf("test_right %d = %15.7e \n",i+1,test_right[i]);
		}
	}

//	printf("inverse1 done! \n");
}

__global__ void	currentRefresh(float* right_side,float* effect_flux, float* current_representation)
{
	int i,num;
	float sum;
	__shared__ float polyCoef[32][32];
	unsigned int index = threadIdx.y*ELI_NUM + blockIdx.x * 32 + threadIdx.x;
	polyCoef[threadIdx.y][threadIdx.x] = effect_flux[index]*right_side[EF_NFCOIL+threadIdx.y];
	__syncthreads();
	if(threadIdx.y == 0)
	{
		sum = 0;
		num = (int)current_representation[25];
		
		for (i=0;i<num;i++)
		{
			sum = sum + polyCoef[i][threadIdx.x];
		}
		effect_flux[blockIdx.x*32+threadIdx.x] = sum;	
	}
}

int shot_set()
{
	float* cpos = (float*)malloc(35*20*2*sizeof(float));
	int* inter_pos = (int*)malloc(2*35*20*sizeof(float));
	float* inter_arg = (float*)malloc(4*35*20*sizeof(float));
	float* grid_pos = (float*)malloc(35*2*2*sizeof(float));
	keyTimes = (float*)malloc(35*sizeof(float));
//	init_read_rfm(&num_shapes, &num_seg, cpos, keyTimes, \
			grid_pos,inter_pos,inter_arg);

	cudaMalloc((void**)&cu_inter_pos,(2*num_shapes*num_seg*sizeof(float)));

	cudaMalloc((void**)&cu_inter_arg,(4*num_shapes*num_seg*sizeof(float)));

	cudaMalloc((void**)&cu_grid_pos,(30*2*2*sizeof(float)));

	cudaMemcpy(cu_inter_pos, inter_pos, 2*num_shapes*num_seg*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_inter_arg, inter_arg, 4*num_shapes*num_seg*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_grid_pos, grid_pos, num_shapes*2*2*sizeof(float), cudaMemcpyHostToDevice);
	printf("\ncpos: %d %d\n",num_shapes , num_seg);
	free(cpos);
//	printf("still OK\n");
	free(inter_pos);
	free(inter_arg);
	free(grid_pos);
	printf("shot set succed\n");
	return 0;
}

__global__ void error_amend3(float* effect_flux, float* leas_flux)
{
	float k = leas_flux[0];
//	printf("%15.7e\n",k);
	unsigned int index = blockIdx.x * 32 + threadIdx.x;
	if(blockIdx.x<ELI_CORE)
	{
		effect_flux[index] *= k;
	}
	else 
	{
		if(threadIdx.x<ELI_REST)		
			effect_flux[index] *= k;
	}
}

__global__ void boundflux(float* cu_boundresponse, float* effect_flux,float* leas_flux)
{
	float current;
	float result=0;
	int index;
	for(int i=0;i<32;i++)
	{
		current = effect_flux[blockIdx.x*32+i];
		index = (blockIdx.x*32+i)*1024+threadIdx.x;
		result += cu_boundresponse[index]*current;
	}
	leas_flux[blockIdx.x*1024+threadIdx.x] = result;
}

__global__ void post_boundflux(float* leas_flux,float* cuBoundResult)
{
	float result=0;
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

__global__ void final_boundflux(float* leas_flux,float* cuBoundResult)
{
	float result=0;
	int index;
	for(int i=0;i<35;i++)
	{
		index = i*1024+threadIdx.x;
		result += leas_flux[index];
	}
	cuBoundResult[blockIdx.x*1024+threadIdx.x] = result;
}

__global__ void data_regroup(float* effect_flux,float* leas_flux,int* cu_num_se)
{
	int offset,a;
	__shared__ float data_buffer[mesh1];
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

__global__ void gs_solver_init(float* cuBoundResult, float* cu_bfft)
{
	__shared__ float sbfft[mesh];
	sbfft[threadIdx.x] = cu_bfft[blockIdx.x*mesh1+threadIdx.x]/(delR*delZ);
	__syncthreads();
	if(blockIdx.x == 0)
	{
		sbfft[threadIdx.x] = (1+delR/(2*(igridr+delR)))*cuBoundResult[threadIdx.x];		
	}
	else if(blockIdx.x == mesh-1)
		sbfft[threadIdx.x] = (1-delR/(2*(egridr-delR)))*cuBoundResult[mesh*3+threadIdx.x];
	else
		sbfft[threadIdx.x] *= cu_currentInver[blockIdx.x];
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

__global__ void fft_inver(float* A, float* B, float* C)
{

   	// Block index
   	float Csub = 0;
        int c,k;

   	for (int a=(mesh+1)*BLOCK_SIZE*blockIdx.y, b= BLOCK_SIZE*blockIdx.x; a<=(mesh+1)*BLOCK_SIZE*blockIdx.y+mesh;a+=BLOCK_SIZE,b+=BLOCK_SIZE*(mesh+1)) 
	{
         	__shared__ float As[BLOCK_SIZE][BLOCK_SIZE+1];
		__shared__ float Bs[BLOCK_SIZE][BLOCK_SIZE+1];

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

__global__ void tri_solver(float* cu_result, float* cu_act_trace, float* cu_tri_main,float* cu_bfft)
{
    	__shared__ float sharedata[mesh][3];
//    	__shared__ float act_trace[mesh];
//    	__shared__ float tri_main[mesh];
//    	__shared__ float result[mesh+1];
//    	__shared__ float tri_up[mesh];
    	unsigned int index=blockIdx.x*(mesh+1)+threadIdx.x;
    	sharedata[threadIdx.x][0]=cu_result[index];
    	unsigned int index2=blockIdx.x*mesh+threadIdx.x;
    	sharedata[threadIdx.x][1]=cu_act_trace[index2];
    	sharedata[threadIdx.x][2]=cu_tri_main[index2];
//    	tri_up[threadIdx.x]=cu_tri_up[threadIdx.x];
//    	sharedata[mesh][0]=0;    
    	__syncthreads();
    	
	int i,k;
    	float tempa,tempb;
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

__global__ void	pfFlux(float* right_side, float* cu_icresponse, float* cu_icflux)
{
	__shared__ float sresponse[48][18];
	unsigned int index = blockIdx.x * 864 + threadIdx.y*18 + threadIdx.x;
//	sresponse[threadIdx.y][threadIdx.x] = ;
//	__syncthreads();
	sresponse[threadIdx.y][threadIdx.x] = right_side[threadIdx.x] * cu_icresponse[index]; 
	__syncthreads();

//	if(blockIdx.x == 0 && threadIdx.y == 0)
//	{
//		printf("pfflux %d %f \n",threadIdx.x,right_side[threadIdx.x]);
//	}

//	float temp_sum = 0;

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

__global__ void	ecFlux(float* test_diagnostic2, float* cu_ecgridresponse, float* cu_ecflux)
{
	__shared__ float sresponse[129][6];
	unsigned int index = blockIdx.x * 774 + threadIdx.y*6 + threadIdx.x;
	sresponse[threadIdx.y][threadIdx.x] = cu_ecgridresponse[index];
	__syncthreads();
	sresponse[threadIdx.y][threadIdx.x] *= test_diagnostic2[threadIdx.x]; 
	__syncthreads();
	float temp_sum = 0;
	if(threadIdx.x == 0)
	{
		temp_sum = sresponse[threadIdx.y][0]+sresponse[threadIdx.y][1]+sresponse[threadIdx.y][2]+sresponse[threadIdx.y][3]+sresponse[threadIdx.y][4]+sresponse[threadIdx.y][5];	
		cu_ecflux[blockIdx.x*129+threadIdx.y]=temp_sum;
	}
}

__global__ void	Fluxadd(float* cu_icflux,float* cu_ecflux,float* cu_result)
{
	if(blockIdx.x>0&&blockIdx.x<256&&threadIdx.x>0&&threadIdx.x<256)
	{
		cu_result[(blockIdx.x-1)*256+threadIdx.x-1] += cu_icflux[blockIdx.x*257+threadIdx.x]+cu_ecflux[blockIdx.x*257+threadIdx.x];
//		cu_result[(blockIdx.x-1)*256+threadIdx.x-1] = -1* cu_result[(blockIdx.x-1)*256+threadIdx.x-1];
	}
}

__global__ void truncate(float* dmate)
{
	float k = truncation;
	unsigned int index = blockIdx.x * UNKNOW_ALL + blockIdx.x;
	dmate[index] += k;
}

int result_output(int m)
{
	int factor;

	if (current_representation[0] == 1 || current_representation[0] == 2)
	{
		cudaMemcpy(result_bound,cu_bound_flux+m, 1*sizeof(float),cudaMemcpyDeviceToHost);  // boundary psi
		cudaMemcpy(result_nullpoint, xpoint+12*m, 12*sizeof(float),cudaMemcpyDeviceToHost);  // (x & axi) (psi & locations)
		
		if (current_representation[0] == 1)
		{
			prime_profile_compute<<<2,mesh2>>>(test_right,cu_profile_result,cu_current_representation);   // ff' p' profile
		}
		if (current_representation[0] == 2)
		{
			prime_profile_compute_spline<<<2,mesh2>>>(test_right,cu_profile_result,cu_current_representation,cu_spline_psi);   // ff' p' profile
		}

		profile_compute<<<2,mesh2>>>(cu_profile_result,xpoint+12*m,cu_bound_flux+m,test_btor,cu_current_representation);   // pres f profile
		cudaMemcpy(result_fpol_profile, cu_profile_result, mesh2*sizeof(float),cudaMemcpyDeviceToHost); //fpol profile
		cudaMemcpy(result_pres_profile, cu_profile_result+mesh2, mesh2*sizeof(float),cudaMemcpyDeviceToHost); //pres profile
		cudaMemcpy(result_ffprime_profile, cu_profile_result+mesh2*2, mesh2*sizeof(float),cudaMemcpyDeviceToHost); //ffprime profile
		cudaMemcpy(result_pprime_profile, cu_profile_result+mesh2*3, mesh2*sizeof(float),cudaMemcpyDeviceToHost); //pprime profile

		result_regroup<<<mesh2,mesh2>>>(cu_result+m*mesh1*mesh1,cu_result_output,cuBoundResult,cu_icflux,cu_ecflux); 
		cudaMemcpy(result_psi,cu_result_output, mesh2*mesh2*sizeof(float),cudaMemcpyDeviceToHost);

		diagnostic_compute<<<EF_DIAG,EF_NFCOIL+(int)current_representation[25]>>>(cu_gfc,test_right,cu_compute_diagnotics,cu_current_representation);  // compute diagnostic signals
		cudaMemcpy(compute_diagnotics,cu_compute_diagnotics, EF_DIAG*sizeof(float),cudaMemcpyDeviceToHost);

		chi_compute<<<1,EF_DIAG>>>(cu_compute_diagnotics,cu_pro_diag,fwt,leas_flux);

		searchbdry<<<numbdry,128>>>(cu_result+m*mesh1*mesh1,cu_xlocation,xpoint+12*m,cu_bound_flux+m,cu_tani,cu_coti,xybdry,maxgridbdry);	 // bdy location compute
		postsearchbdry<<<1,numbdry+3>>>(xybdry,cu_xlocation,xpoint+12*m,cu_bound_flux+m);

		search_surface<<<numbdry,128>>>(cu_result+m*mesh1*mesh1,cu_xlocation,cu_tani,cu_coti,xpoint+12*m,maxgridbdry,xysurface);
		q_profile<<<mesh1,numbdry>>>(cu_result+m*mesh1*mesh1,xpoint+12*m,xysurface,test_btor,cu_qpsi); // q profile
		cudaMemcpy(result_qpsi,cu_qpsi, mesh2*sizeof(float),cudaMemcpyDeviceToHost);
//		result_qpsi[0] = ((1.69550002*fabs(host_btor[0]))/(result_nullpoint[9]*result_nullpoint[9]))/((hright[18]*result_nullpoint[9]+hright[20]/result_nullpoint[9])*2*pi/(10000000*darea));

		forbetapli<<<mesh,mesh>>>(cu_result+m*mesh1*mesh1,leas_flux,cu_profile_result,cu_bound_flux+m,xpoint+12*m);
		postbetapli<<<1,mesh>>>(leas_flux,xybdry,cu_compute_diagnotics);	// beta li
	
		cudaMemcpy(result_xybdry, xybdry, (numbdry+3)*2*sizeof(float),cudaMemcpyDeviceToHost);  
		cudaMemcpy(result_parameter, leas_flux+590000, 8*sizeof(float),cudaMemcpyDeviceToHost);
		cudaMemcpy(converge_error,cu_converge_error, 1*sizeof(float),cudaMemcpyDeviceToHost);
		cudaMemcpy(hright, test_right, (EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);

		if (compute_diagnotics[120] > 0)
		{
			factor = -1;
		}
		else
		{
			factor = 1;
		}
	}

	if (current_representation[0] == 3)
	{
		cudaMemcpy(result_bound_fb,cu_bound_flux_fb+m, 1*sizeof(double),cudaMemcpyDeviceToHost);  // boundary psi
		cudaMemcpy(result_nullpoint_fb, xpoint_fb+12*m, 12*sizeof(double),cudaMemcpyDeviceToHost);  // (x & axi) (psi & locations)

		cudaMemcpy(cu_profile_result_fb, temp_profile_result_fb, 4*mesh2*sizeof(double), cudaMemcpyHostToDevice);

		cudaMemcpy(result_fpol_profile_fb, cu_profile_result_fb, mesh2*sizeof(double),cudaMemcpyDeviceToHost); //fpol profile
		cudaMemcpy(result_pres_profile_fb, cu_profile_result_fb+mesh2, mesh2*sizeof(double),cudaMemcpyDeviceToHost); //pres profile
		cudaMemcpy(result_ffprime_profile_fb, cu_profile_result_fb+mesh2*2, mesh2*sizeof(double),cudaMemcpyDeviceToHost); //ffprime profile
		cudaMemcpy(result_pprime_profile_fb, cu_profile_result_fb+mesh2*3, mesh2*sizeof(double),cudaMemcpyDeviceToHost); //pprime profile

		result_regroup_fb<<<mesh2,mesh2>>>(cu_result_fb+m*mesh1*mesh1,cu_result_output_fb,cuBoundResult_fb,cu_icflux_fb,cu_ecflux_fb); 
		cudaMemcpy(result_psi_fb,cu_result_output_fb, mesh2*mesh2*sizeof(double),cudaMemcpyDeviceToHost);

//		error_amend<<<ELI_CORE,32>>>(effect_flux_fb,leas_flux_fb);
//		error_amend2<<<1,512>>>(leas_flux_fb,cu_compute_diagnotics_fb);
//		cudaMemcpy(compute_diagnotics_fb,cu_compute_diagnotics_fb, EF_DIAG*sizeof(double),cudaMemcpyDeviceToHost);

		cudaMemcpy(cu_compute_diagnotics_fb+(EF_DIAG-1), current_gafile, 1*sizeof(double), cudaMemcpyHostToDevice);

		searchbdry_fb<<<numbdry,128>>>(cu_result_fb+m*mesh1*mesh1,cu_xlocation_fb,xpoint_fb+12*m,cu_bound_flux_fb+m,cu_tani_fb,cu_coti_fb,xybdry_fb,maxgridbdry);	 // bdy location compute
		postsearchbdry_fb<<<1,numbdry+3>>>(xybdry_fb,cu_xlocation_fb,xpoint_fb+12*m,cu_bound_flux_fb+m);

		search_surface_fb<<<numbdry,128>>>(cu_result_fb+m*mesh1*mesh1,cu_xlocation_fb,cu_tani_fb,cu_coti_fb,xpoint_fb+12*m,maxgridbdry,xysurface_fb);
		q_profile_fb<<<mesh1,numbdry>>>(cu_result_fb+m*mesh1*mesh1,xpoint_fb+12*m,xysurface_fb,test_btor_fb,cu_qpsi_fb); // q profile
		cudaMemcpy(result_qpsi_fb,cu_qpsi_fb, mesh2*sizeof(double),cudaMemcpyDeviceToHost);
//		result_qpsi[0] = ((1.69550002*fabs(host_btor[0]))/(result_nullpoint[9]*result_nullpoint[9]))/((hright[18]*result_nullpoint[9]+hright[20]/result_nullpoint[9])*2*pi/(10000000*darea));

		forbetapli_fb<<<mesh,mesh>>>(cu_result_fb+m*mesh1*mesh1,leas_flux_fb,cu_profile_result_fb,cu_bound_flux_fb+m,xpoint_fb+12*m);
		postbetapli_fb<<<1,mesh>>>(leas_flux_fb,xybdry_fb,cu_compute_diagnotics_fb);	// beta li
	
		cudaMemcpy(result_xybdry_fb, xybdry_fb, (numbdry+3)*2*sizeof(double),cudaMemcpyDeviceToHost);  
		cudaMemcpy(result_parameter_fb, leas_flux_fb+590000, 8*sizeof(double),cudaMemcpyDeviceToHost);
		cudaMemcpy(converge_error_fb,cu_converge_error_fb, 1*sizeof(double),cudaMemcpyDeviceToHost);
		cudaMemcpy(hright_fb, test_right_fb, (EF_NFCOIL+(int)current_representation[25])*sizeof(double),cudaMemcpyDeviceToHost);

		if (current_gafile[0] > 0)
		{
			factor = -1;
		}
		else
		{
			factor = 1;
		}
	}

	

///////////////////////////////////////////write gfile/////////////////////////////////////////////////////////////////////////////
	FILE *gfile;
	char seg[30];
	FILE *afile;
	char sea[30];

	if (current_representation[0] == 1 || current_representation[0] == 2)
	{
		sprintf(seg,"g%d.%05d_pefit",int(ishot[m]),int(itime[m]));
		gfile = fopen(seg,"w");

		sprintf(sea,"a%d.%05d_pefit",int(ishot[m]),int(itime[m]));
		afile = fopen(sea,"w");

		fprintf(gfile,"  EFITD    12/10/2013    #%d  %dms           3 257 257\n",int(ishot[m]),int(itime[m]));
		// line 1  P-EFIT date shot time grid
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",0.170000000E+01, 0.320000000E+01, 0.169550002E+01, 0.840000000E+00, 0.000000000E+00);
		// line 2  drim zdim rcentr rleft zmid
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",result_nullpoint[9],result_nullpoint[10],result_nullpoint[4]*factor,result_bound[0]*factor,host_btor[0]);
		// line 3  rmaxis zmaxis simag sidry bcentr
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",compute_diagnotics[120],result_nullpoint[4]*factor,0.0,result_nullpoint[9],0.0);
		// line 4  curent simag xdum rmaxis xdum
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",result_nullpoint[10],0.0,result_bound[0]*factor,0.0,0.0);
		// line 5  zmaxis xdum sibry xdum xdum
	
		write_gfile_function(257, result_fpol_profile, gfile,factor);  // line 6~18  fpol 
		write_gfile_function(257, result_pres_profile, gfile,1);  // line 19~31  pres 
		write_gfile_function(257, result_ffprime_profile, gfile,factor);  // line 32~44 ffprime 
		write_gfile_function(257, result_pprime_profile, gfile,factor);  // line 45~57   pprime

		write_gfile_function(257*257, result_psi, gfile, factor);     // line 58~902 si
		write_gfile_function(257, result_qpsi, gfile,1);              // line 903~915 q(si)
	
		fprintf(gfile,"   %d   %d\n",int(numbdry+3),int(numlim));
		write_gfile_function((numbdry+3)*2, result_xybdry, gfile,1);   // rbbbs zbbbs
		write_gfile_function(numlim*2, result_xylimtr, gfile,1);   // xlim ylim
 
		fprintf(gfile,"    %d% 10.9E    %d\n", 0, 1.7, 0);
		fprintf(gfile,"    %d\n", 0);

//		write_gfile_function(mesh2*2*numbdry,result_xysurface,gfile,1);   // surface locations for testing
//		fprintf(gfile,"%d %d %d %d\n", 65, 65, int(ishot[no_slice]),int(itime[no_slice]));
//		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E\n", 0.84, 2.54, -1.6, 1.6);

//		write_gfile_function_double(EF_NFCOIL, hright, gfile);   // brsp
//		write_gfile_function(EF_ECCOIL, diagnostic2+m*EF_ECCOIL, gfile);  // ecurrent
//		write_gfile_function(65*65, result_pcurrent, gfile);  //pcurrent	
	}

	if (current_representation[0] == 3)
	{
		sprintf(seg,"g%d.%05d_pefit",int(ishot[m]),int(itime[m]));
		gfile = fopen(seg,"w");

		sprintf(sea,"a%d.%05d_pefit",int(ishot[m]),int(itime[m]));
		afile = fopen(sea,"w");

		fprintf(gfile,"  EFITD    12/10/2013    #%d  %dms           3 257 257\n",int(ishot[m]),int(itime[m]));
		// line 1  P-EFIT date shot time grid
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",0.170000000E+01, 0.320000000E+01, 0.169550002E+01, 0.840000000E+00, 0.000000000E+00);
		// line 2  drim zdim rcentr rleft zmid
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",result_nullpoint_fb[9],result_nullpoint_fb[10],result_nullpoint_fb[4]*factor,result_bound_fb[0]*factor,host_btor_fb[0]);
		// line 3  rmaxis zmaxis simag sidry bcentr
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",current_gafile[0],result_nullpoint_fb[4]*factor,0.0,result_nullpoint_fb[9],0.0);
		// line 4  curent simag xdum rmaxis xdum
		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",result_nullpoint_fb[10],0.0,result_bound_fb[0]*factor,0.0,0.0);
		// line 5  zmaxis xdum sibry xdum xdum
	
		write_gfile_function_double(257, result_fpol_profile_fb, gfile,factor);  // line 6~18  fpol 
		write_gfile_function_double(257, result_pres_profile_fb, gfile,1);  // line 19~31  pres 
		write_gfile_function_double(257, result_ffprime_profile_fb, gfile,factor);  // line 32~44 ffprime 
		write_gfile_function_double(257, result_pprime_profile_fb, gfile,factor);  // line 45~57   pprime

		write_gfile_function_double(257*257, result_psi_fb, gfile, factor);     // line 58~902 si
		write_gfile_function_double(257, result_qpsi_fb, gfile,1);              // line 903~915 q(si)
	
		fprintf(gfile,"   %d   %d\n",int(numbdry+3),int(numlimter_gafile));
		write_gfile_function_double((numbdry+3)*2, result_xybdry_fb, gfile,1);   // rbbbs zbbbs
		write_gfile_function_double(numlimter_gafile*2, limiter_gafile, gfile,1);   // xlim ylim
 
		fprintf(gfile,"    %d% 10.9E    %d\n", 0, 1.7, 0);
		fprintf(gfile,"    %d\n", 0);

//		write_gfile_function(mesh2*2*numbdry,result_xysurface,gfile,1);   // surface locations for testing
//		fprintf(gfile,"%d %d %d %d\n", 65, 65, int(ishot[no_slice]),int(itime[no_slice]));
//		fprintf(gfile,"% 10.9E% 10.9E% 10.9E% 10.9E\n", 0.84, 2.54, -1.6, 1.6);

//		write_gfile_function_double(EF_NFCOIL, hright, gfile);   // brsp
//		write_gfile_function(EF_ECCOIL, diagnostic2+m*EF_ECCOIL, gfile);  // ecurrent
//		write_gfile_function(65*65, result_pcurrent, gfile);  //pcurrent
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////write afile/////////////////////////////////////////////////////////////////////////////

	if (current_representation[0] == 1 || current_representation[0] == 2)
	{
		fprintf(afile," 23-Jan-15 12/10/2013 \n");                                        // line 1
		fprintf(afile," %d               1\n", int(ishot[m]));			// line 2
		fprintf(afile,"  %10.9E\n", float(itime[m]));				//line 3	
		fprintf(afile,"*%.3f             1               0 SNT   3   2 CLC    39   40\n", itime[m]);	//line 4
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", result_parameter[6],0.169550002E+03,host_btor[0],diagnostic1[120]); //line 5
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", compute_diagnotics[120],0.0,0.0,0.0); //line 6
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,result_parameter[0]*1000000);	//line 7
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0);	//line 8
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", result_parameter[1],result_parameter[3],0.0,0.0); //line 9
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.8*result_qpsi[243]+0.2*result_qpsi[244],0.0); //line 10
		fprintf(afile," % 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0); //line 11
		fprintf(afile," % 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0); //line 12
		fprintf(afile," % 10.9E% 10.9E\n", 0.0,0.0); //line 13
		fprintf(afile," % 10.9E% 10.9E\n", 0.0,0.0); //line 14
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 15
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,result_qpsi[256],0.0,0.0); //line 16
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,result_bound[0],result_parameter[2]*10000,result_parameter[5]); //line 17
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", converge_error[0]/result_nullpoint[11],0.0,result_qpsi[0],0.0); //line 18
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 19
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 20
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 21
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 22
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", result_nullpoint[9]*100,result_nullpoint[10]*100,result_nullpoint[4],0.0); //line 23
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 24
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 25
		fprintf(afile,"    %d   %d   %d    %d\n", (int)EF_NSILOP,(int)EF_MAGPRI,(int)EF_NFCOIL,(int)EF_ECCOIL); //line 26
		write_afile_function(EF_NSILOP+EF_MAGPRI, compute_diagnotics, afile);   // magnetic diagnostic
		write_afile_function(EF_NFCOIL, hright, afile);   // brsp
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); 
		fprintf(afile," % 10.9E% 10.9E\n", 0.0,0.0);
	}

	if (current_representation[0] == 3)
	{
		fprintf(afile," 23-Jan-15 12/10/2013 \n");                                        // line 1
		fprintf(afile," %d               1\n", int(ishot[m]));			// line 2
		fprintf(afile,"  %10.9E\n", float(itime[m]));				//line 3	
		fprintf(afile,"*%.3f             1               0 SNT   3   2 CLC    39   40\n", itime[m]);	//line 4
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", result_parameter_fb[6],0.169550002E+03,host_btor_fb[0],diagnostic1_fb[120]); //line 5
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", current_gafile[0],0.0,0.0,0.0); //line 6
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,result_parameter_fb[0]*1000000);	//line 7
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0);	//line 8
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", result_parameter_fb[1],result_parameter_fb[3],0.0,0.0); //line 9
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.8*result_qpsi_fb[243]+0.2*result_qpsi_fb[244],0.0); //line 10
		fprintf(afile," % 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0); //line 11
		fprintf(afile," % 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0); //line 12
		fprintf(afile," % 10.9E% 10.9E\n", 0.0,0.0); //line 13
		fprintf(afile," % 10.9E% 10.9E\n", 0.0,0.0); //line 14
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 15
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,result_qpsi_fb[256],0.0,0.0); //line 16
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,result_bound_fb[0],result_parameter_fb[2]*10000,result_parameter_fb[5]); //line 17
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", converge_error_fb[0]/result_nullpoint_fb[11],0.0,result_qpsi_fb[0],0.0); //line 18
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 19
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 20
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 21
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 22
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", result_nullpoint_fb[9]*100,result_nullpoint_fb[10]*100,result_nullpoint_fb[4],0.0); //line 23
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 24
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); //line 25
		fprintf(afile,"    %d   %d   %d    %d\n", (int)EF_NSILOP,(int)EF_MAGPRI,(int)EF_NFCOIL,(int)EF_ECCOIL); //line 26
		write_afile_function_double(EF_NSILOP+EF_MAGPRI, compute_diagnotics_fb, afile);   // magnetic diagnostic
		write_afile_function_double(EF_NFCOIL, hright_fb, afile);   // brsp
		fprintf(afile," % 10.9E% 10.9E% 10.9E% 10.9E\n", 0.0,0.0,0.0,0.0); 
		fprintf(afile," % 10.9E% 10.9E\n", 0.0,0.0);
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}

int init(int argc, char* argv[])
{
	execute_mode_init();
	printf("execute_mode_init is ok! \n");

	if (current_representation[0] == 1 || current_representation[0] == 2)
	{
		tri_solver_init();
    		read_gfunction_files();
		test_diag_read(argc, argv);
	}
	if (current_representation[0] == 3)
	{
		tri_solver_init_fb();
		read_gfunction_files_fb();
//		test_diag_read(argc, argv);
	}
	printf("current_representation is ok! \n");

	read_setup_files();
	printf("read_setup_files is ok! \n");

    	allocate_inter_vars();
	printf("allocate_inter_vars is ok! \n");

    	arrange_init();
	
	return 0;
}

int pefit_init(int argc, char* argv[],int m)
{
//	cudaFuncSetCacheConfig(tridaig_mul, cudaFuncCachePreferShared);
//    	float alpha2=1.0;
//    	float beta2=0.0;
//    	cublasHandle_t handle;
//   	cublasCreate(&handle);
//	cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, EF_NFCOIL+(int)current_representation[25], EF_NFCOIL+(int)current_representation[25],  EF_DIAG, &alpha2, cu_dgfc, EF_DIAG, cu_dgfc, EF_DIAG, &beta2,dmata , EF_NFCOIL+(int)current_representation[25]);
	if (current_representation[0] == 1 || current_representation[0] == 2)
	{
		if (current_representation[0] == 1)
		{
			data_eli_init<<<mesh,mesh>>>(cu_initNormFlux, effect_flux, cu_num_se, cu_current_representation);
		}
		if (current_representation[0] == 2)
		{
			data_eli_init_spline<<<mesh,mesh>>>(cu_psi_gafile, effect_flux, cu_num_se, cu_current_representation, cu_psiaxis_gafile,cu_psibdry_gafile);
			printf("data_eli_init_spline is ok\n");
		}
/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result1= (float *)malloc(ELI_NUM * (int)current_representation[25]*sizeof(float));
			cudaMemcpy(result1, effect_flux, ELI_NUM * (int)current_representation[25]*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result1;
			Result1=fopen("effect_flux.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<ELI_NUM * (int)current_representation[25];i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result1,"%15.7e\n",result1[i]);
			fclose(Result1);
*////////////////////////////////////////////////////////////////////////////////////////////////////////

		compute_converge_init<<<mesh,mesh>>>(cu_result+m*mesh1*mesh1,cu_result_compare,cu_converge);
		cudaDeviceSynchronize();
		tridaig_mul<<<grid3,threads3>>>(effect_flux, cu_efDiagresponse, leas_flux, cu_current_representation);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result11 = (float *)malloc(1116 * EF_DIAG *(int)current_representation[25]*sizeof(float));
			cudaMemcpy(result11, leas_flux, 1116 * EF_DIAG *(int)current_representation[25]*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result11;
			Result11=fopen("leas_flux.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<1116 * EF_DIAG *(int)current_representation[25];i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result11,"%15.7e\n",result11[i]);
			fclose(Result11);
*////////////////////////////////////////////////////////////////////////////////////////////////////////

		mul_postprocess<<<(int)current_representation[25],EF_DIAG>>>(leas_flux, responsematrix, cu_current_representation);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result2= (float *)malloc(EF_DIAG * (int)current_representation[25]*sizeof(float));
			cudaMemcpy(result2, responsematrix, EF_DIAG * (int)current_representation[25]*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result2;
			Result2=fopen("responsematrix.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<EF_DIAG * (int)current_representation[25];i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result2,"%15.7e\n",result2[i]);
			fclose(Result2);
*////////////////////////////////////////////////////////////////////////////////////////////////////////

		data_process_ec<<<EF_DIAG,EF_ECCOIL>>>(cu_ec_g,test_diagnostic2,test_ec);
		data_process_right<<<EF_DIAG,1>>>(fwt,test_diagnostic1,test_diagbit,test_ec,test_weight,cu_pro_diag,test_psiref);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result3= (float *)malloc((EF_DIAG)*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result3, cu_gfc, (EF_DIAG)*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result3;
			Result3=fopen("responsematrix_1.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<(EF_DIAG)*(EF_NFCOIL+(int)current_representation[25]);i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result3,"% 10.9E\n",result3[i]);
			fclose(Result3);
*////////////////////////////////////////////////////////////////////////////////////////////////////////

		data_process_left<<<EF_NFCOIL+(int)current_representation[25],EF_DIAG>>>(fwt,cu_gfc,cu_dgfc,cu_current_representation);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result4= (float *)malloc((EF_DIAG+2)*23*sizeof(float));
			cudaMemcpy(result4, cu_dgfc, (EF_DIAG+2)*23*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result4;
			Result4=fopen("responsematrix_2.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<(EF_DIAG+2)*23;i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result4,"% 10.9E\n",result4[i]);
			fclose(Result4);
*////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (current_representation[26]>0 && current_representation[0] == 1)
		{
			jtrepresnt<<<((int)current_representation[1]+1)*(1-(int)current_representation[3])+((int)current_representation[2]+1)*(1-(int)current_representation[4]),1>>>(cu_dgfc,cu_pro_diag,cu_current_representation);
		}
		if (current_representation[0] == 2)
		{
			if (PRES_FIT == 1)  
			{
				pprime_flux_function_integral_constrain<<<1,mpress[0]>>>(cu_current_representation,cu_pres_psi,cu_pres_left,cu_spline_psi);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result22= (float *)malloc(mpress[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result22, cu_pres_left, mpress[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result22;
			Result22=fopen("pres_left.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mpress[0]*(EF_NFCOIL+(int)current_representation[25]);i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result22,"%15.7e\n",result22[i]);
			fclose(Result22);
////////////////////////////////////////////////////////////////////////////////////////////////////////				
			}
			if (JPHI_FIT == 1)
			{
//				jphi_constrain_init_1<<<numbdry,100>>>(cu_initNormFlux,cu_xlocation,cu_tani,cu_coti,xpoint+12*m,cu_num_jphi,cu_jphi_npsi,cu_jphi_surface);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result12= (float *)malloc(mjphi[0]*numbdry*2*sizeof(float));
			cudaMemcpy(result12, cu_jphi_surface, mjphi[0]*numbdry*2*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result12;
			Result12=fopen("jphi_surface.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mjphi[0]*numbdry*2;i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result12,"%15.7e\n",result12[i]);
			fclose(Result12);
////////////////////////////////////////////////////////////////////////////////////////////////////////
//				jphi_constrain_2<<<mjphi[0],numbdry>>>(cu_jphi_npsi,cu_jphi_surface,cu_current_representation,cu_spline_psi,cu_jphi_left);
			}
			spline_constrain<<<num_knots[0]-2,dim3(2,2,1)>>>(cu_current_representation,cu_spline_psi,cu_spline_left,cu_spline_right);
			printf("spline_constrain is ok\n");
		}
		if (MSE_FIT == 1)
		{
			
		}

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result5= (float *)malloc(EF_DIAG*34*sizeof(float));
			cudaMemcpy(result5, cu_dgfc, EF_DIAG*34*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result5;
			Result5=fopen("responsematrix_3.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<EF_DIAG*34;i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result5,"% 10.9E\n",result5[i]);
			fclose(Result5);
*////////////////////////////////////////////////////////////////////////////////////////////////////////

//	rightHandSideInit<<<EF_NFCOIL+(int)current_representation[25],EF_DIAG+EF_NFCOIL+(int)current_representation[26]>>>(cu_dgfc,cu_pro_diag,right_side,cu_current_representation);
//	leftHandSideInit<<<dim3(EF_NFCOIL+current_representation[25],EF_NFCOIL+current_representation[25]),threads10>>>(cu_dgfc,temp_dmate,cu_current_representation);
//	leftInit<<<(EF_NFCOIL+current_representation[25])*(EF_NFCOIL+current_representation[25]),32>>>(temp_dmate,dmata);
//	truncate<<<UNKNOW_ALL,1>>>(dmata);
		
		if (current_representation[0] == 1)
		{
			changePos_1<<<(EF_NFCOIL+(int)current_representation[25]),(EF_DIAG+(int)current_representation[26])>>>(cu_dgfc,d_A,cu_current_representation);
		}
		if (current_representation[0] == 2)
		{
			change_Pos_2<<<(EF_NFCOIL+(int)current_representation[25]),EF_DIAG>>>(cu_dgfc,h_A_left,cu_current_representation,cu_diag_num_all,cu_pro_diag,cu_diag_right);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result41 = (float *)malloc(diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result41, h_A_left, diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result41;
			Result41=fopen("responsematrix_new_01.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25]);i++)
			{
				fprintf(Result41,"%15.7e\n",result41[i]);
//				printf("Result4 = % 10.9E\n",result4[i]);
			}
			fclose(Result41);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (MSE_FIT == 1)
			{
//				change_Pos_mse<<<,>>>();
			}
			if (JPHI_FIT == 1)
			{
				change_Pos_jphi<<<(EF_NFCOIL+(int)current_representation[25]),mjphi[0]>>>(cu_jphi_left,h_A_left,cu_j_phi,cu_j_phi_uncer,cu_diag_right,cu_current_representation,cu_diag_num_all,offset_dA_jphi);

////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result42 = (float *)malloc(diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result42, h_A_left, diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result42;
			Result42=fopen("responsematrix_new_02.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25]);i++)
			{
				fprintf(Result42,"%15.7e\n",result42[i]);
//				printf("Result4 = % 10.9E\n",result4[i]);
			}
			fclose(Result42);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			}
			if (PRES_FIT == 1)  
			{
				change_Pos_pres_init<<<(EF_NFCOIL+(int)current_representation[25]),mpress[0]>>>(cu_pres_left,h_A_left,cu_press,cu_press_uncer,cu_diag_right,cu_current_representation,cu_diag_num_all,offset_dA_pres);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result43 = (float *)malloc(diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result43, h_A_left, diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result43;
			Result43=fopen("responsematrix_new_03.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25]);i++)
			{
				fprintf(Result43,"%15.7e\n",result43[i]);
//				printf("Result4 = % 10.9E\n",result4[i]);
			}
			fclose(Result43);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
			change_Pos_spline<<<(EF_NFCOIL+(int)current_representation[25]),(num_knots[0]-2)*6>>>(cu_spline_left,h_A_left,cu_current_representation,cu_diag_num_all,offset_dA_spline,cu_diag_right);

////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result5 = (float *)malloc(diag_num_all[0]*sizeof(float));
			cudaMemcpy(result5, cu_diag_right, diag_num_all[0]*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result5;
			Result5 = fopen("cu_diag_right.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<diag_num_all[0];i++)
			{
				fprintf(Result5,"%15.7e\n",result5[i]);
//				printf("Result4 = % 10.9E\n",result4[i]);
			}
			fclose(Result5);
////////////////////////////////////////////////////////////////////////////////////////////////////////
		}

		cudaDeviceSynchronize();
////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result4 = (float *)malloc(diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result4, h_A_left, diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result4;
			Result4=fopen("responsematrix_new.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25]);i++)
			{
				fprintf(Result4,"%15.7e\n",result4[i]);
//				printf("Result4 = % 10.9E\n",result4[i]);
			}
			fclose(Result4);
////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		if (current_representation[0] == 1)
			info = LAPACKE_sgesvd(LAPACK_COL_MAJOR, 'A', 'A', (EF_DIAG+(int)current_representation[26]), (EF_NFCOIL+(int)current_representation[25]), h_A, (EF_DIAG+(int)current_representation[26]), h_S, h_U, (EF_DIAG+(int)current_representation[26]), h_VT, (EF_NFCOIL+(int)current_representation[25]), superb);

		if (current_representation[0] == 2)
		{
			cudaMemcpy(h_A, h_A_left, diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			info = LAPACKE_sgesvd(LAPACK_COL_MAJOR, 'A', 'A', diag_num_all[0], (EF_NFCOIL+(int)current_representation[25]), h_A, diag_num_all[0], h_S, h_U, diag_num_all[0], h_VT, (EF_NFCOIL+(int)current_representation[25]), superb);
		}

//	LDLTDCMP(UNKNOW_ALL,hmata);
//	LDLTBKSB(UNKNOW_ALL, hmata, hright);
	
	for(int i=0;i<EF_NFCOIL+current_representation[25];i++)
		printf("singular value (%d) = % 10.9E \n",i+1,h_S[i]);
	printf("\n\n");

/*//////////////////////////////////
	FILE *Result6;
	Result6=fopen("u_matrix.txt","w");
	for(int i=0;i<127*127;i++)
	{
		fprintf(Result6,"% 10.9E\n",h_U[i]);
	}
	fclose(Result6);

	FILE *Result7;
	Result7=fopen("vt_matrix.txt","w");
	for(int i=0;i<(EF_NFCOIL+current_representation[25])*(EF_NFCOIL+current_representation[25]);i++)
	{
		fprintf(Result7,"% 10.9E\n",h_VT[i]);
	}
	fclose(Result7);
*///////////////////////////////////

		if (current_representation[0] == 1)
			cudaMemcpyAsync(d_U,h_U,(EF_DIAG+(int)current_representation[26])*(EF_DIAG+(int)current_representation[26])*sizeof(float),cudaMemcpyHostToDevice,0);
		if (current_representation[0] == 2)
			cudaMemcpyAsync(d_U,h_U,diag_num_all[0]*diag_num_all[0]*sizeof(float),cudaMemcpyHostToDevice,0);
		cudaMemcpyAsync(d_S,h_S,(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyHostToDevice,0);
		cudaMemcpyAsync(d_VT,h_VT,(EF_NFCOIL+(int)current_representation[25])*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyHostToDevice,0);
		cudaDeviceSynchronize();
		
		if (current_representation[0] == 1)
			changePos_2<<<(EF_DIAG+(int)current_representation[26]),(EF_DIAG+(int)current_representation[26])>>>(d_U,d_S,d_VT,cu_current_representation,cu_U_matrix,cu_S_vector,cu_VT_matrix);
		if (current_representation[0] == 2)
			changePos_3<<<diag_num_all[0],diag_num_all[0]>>>(d_U,d_S,d_VT,cu_current_representation,cu_U_matrix,cu_S_vector,cu_VT_matrix,cu_diag_num_all);
		
		if (current_representation[0] == 1)
			responsematrix_inverse_U<<<(EF_DIAG+(int)current_representation[26]),(EF_DIAG+(int)current_representation[26])>>>(cu_pro_diag,cu_U_matrix,cu_current_representation,cu_response_inverse_temp);
		if (current_representation[0] == 2)
			responsematrix_inverse_U_spline<<<diag_num_all[0],diag_num_all[0]>>>(cu_diag_right,cu_U_matrix,cu_current_representation,cu_response_inverse_temp,cu_diag_num_all);

       		responsematrix_inverse_S<<<1,(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp,cu_S_vector,cu_current_representation);	
		responsematrix_inverse_VT<<<(EF_NFCOIL+(int)current_representation[25]),(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp,cu_VT_matrix,test_right,cu_current_representation);

		currentRefresh<<<grid6,dim3(32,current_representation[25])>>>(test_right,effect_flux,cu_current_representation);
/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result8= (float *)malloc(ELI_NUM *sizeof(float));
			cudaMemcpy(result8, effect_flux, ELI_NUM *sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result8;
			Result8=fopen("effect_ip_current.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<ELI_NUM;i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result8,"%15.7e\n",result8[i]);
			fclose(Result8);
*////////////////////////////////////////////////////////////////////////////////////////////////////////
//	cublasDestroy(handle);
	}

	if (current_representation[0] == 3)
	{
		data_eli_init_fb<<<mesh,mesh>>>(cu_psi_gafile, effect_flux_fb, cu_num_se, cu_ffprime_gafile, cu_pprime_gafile, cu_r_bdry_gafile, cu_z_bdry_gafile, cu_numbdry_gafile, cu_psiaxis_gafile, cu_psibdry_gafile);

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
			double* result1= (double *)malloc(ELI_NUM*sizeof(double));
			cudaMemcpy(result1, effect_flux_fb, ELI_NUM*sizeof(double),cudaMemcpyDeviceToHost);
			FILE *Result1;
			double xxxyyy = 0;
			Result1=fopen("effect_flux.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<ELI_NUM;i++)
			{
				fprintf(Result1,"%15.7e\n",result1[i]);
				xxxyyy = xxxyyy + result1[i];
			}
			printf("PEFIT FIXED INIT PLASMA CURRENT1 = %f \n",xxxyyy);
			fclose(Result1);
*///////////////////////////////////////////////////////////////////////////////////////////////////////

		compute_converge_init_fb<<<mesh,mesh>>>(cu_result_fb+m*mesh1*mesh1,cu_result_compare_fb,cu_converge);

		boundflux_fb<<<1116,1024>>>(cu_boundresponse_fb, effect_flux_fb, cuLbuffer_fb);
		post_boundflux_fb<<<35,1024>>>(cuLbuffer_fb,cubuffer_fb);
		final_boundflux_fb<<<1,1024>>>(cubuffer_fb,cuBoundResult_fb);
		data_regroup_fb<<<mesh1,mesh1>>>(effect_flux_fb,cu_bfft_fb,cu_num_se);
		gs_solver_init_fb<<<mesh,mesh>>>(cuBoundResult_fb,cu_bfft_fb);
		fft_inver_fb<<<grid,threads>>>(cu_bfft_fb,cu_mat_fb,cu_ipflux_fb);
		tri_solver_fb<<<grid2,threads2>>>(cu_ipflux_fb,cu_act_trace_fb,cu_tri_main_fb,cu_bfft_fb);
		fft_inver_fb<<<grid,threads>>>(cu_mat_fb,cu_bfft_fb,cu_ipflux_fb);

/*////////////////////////////////////////////////////////////////////////////////////////////////////
			double* result12= (double *)malloc(mesh1*mesh1*sizeof(double));
			cudaMemcpy(result12,cu_ipflux_fb, mesh1*mesh1*sizeof(double),cudaMemcpyDeviceToHost);
	   		FILE *Result12;
			Result12=fopen("psi_ip.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh;i++)
				for(int j=0;j<mesh;j++)
					fprintf(Result12,"%15.7e\n",result12[i*mesh1+j]);
			fclose(Result12);
*///////////////////////////////////////////////////////////////////////////////////////////////////////

		ecFlux_fb<<<grid17,threads17>>>(cu_ecurrent_gafile,cu_ecgridresponse_fb,cu_ecflux_fb);

/*/////////////////////////////////////////////////////////////////////////////////////////////////////
			double* result14= (double *)malloc(mesh2*mesh2*sizeof(double));
			cudaMemcpy(result14, cu_ecflux_fb, mesh2*mesh2*sizeof(double),cudaMemcpyDeviceToHost);
	   		FILE *Result14;
			Result14=fopen("psi_ec.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh2;i++)
				for(int j=0;j<mesh2;j++)
					fprintf(Result14,"%15.7e\n",result14[i*mesh2+j]);
			fclose(Result14);
*///////////////////////////////////////////////////////////////////////////////////////////////////////

		data_process_fb_right<<<1,numbdry_gafile-1>>>(cu_ipflux_fb,cu_ecflux_fb,cu_psibdry_gafile,cu_r_bdry_gafile,cu_z_bdry_gafile,cu_pro_diag_fb,cu_fwtbdry);
		data_process_fb_left<<<EF_NFCOIL,numbdry_gafile-1>>>(cu_icresponse_fb,cu_r_bdry_gafile,cu_z_bdry_gafile,cu_dgfc_fb,cu_numbdry_gafile,cu_fwtbdry);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			double* result3 = (double *)malloc((numbdry_gafile-1)*sizeof(double));
			cudaMemcpy(result3, cu_pro_diag_fb, (numbdry_gafile-1)*sizeof(double),cudaMemcpyDeviceToHost);
			FILE *Result3;
			Result3=fopen("cu_pro_diag.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<(numbdry_gafile-1);i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result3,"% 10.9E\n",result3[i]);
			fclose(Result3);
*////////////////////////////////////////////////////////////////////////////////////////////////////////	

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			double* result4 = (double *)malloc(numbdry_gafile*EF_NFCOIL*sizeof(double));
			cudaMemcpy(result4, cu_dgfc_fb, (numbdry_gafile-1)*EF_NFCOIL*sizeof(double),cudaMemcpyDeviceToHost);
			FILE *Result4;
			Result4=fopen("responsematrix_2.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<(numbdry_gafile-1)*EF_NFCOIL;i++)
			{
				fwrite(&result4[i],8,1,Result4);
//				printf("Result4 = % 10.9E\n",result4[i]);
			}
			fclose(Result4);
*////////////////////////////////////////////////////////////////////////////////////////////////////////

		cudaMemcpy(h_A_fb, cu_dgfc_fb, (EF_DIAG+(int)current_representation[26])*(EF_NFCOIL+(int)current_representation[25])*sizeof(double),cudaMemcpyDeviceToHost);

		info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A', (EF_DIAG+(int)current_representation[26]), (EF_NFCOIL+(int)current_representation[25]), h_A_fb, (EF_DIAG+(int)current_representation[26]), h_S_fb, h_U_fb, (EF_DIAG+(int)current_representation[26]), h_VT_fb, (EF_NFCOIL+(int)current_representation[25]), superb_fb);

/*/////////////////////////////////
		FILE *Result5;
		Result5=fopen("s_vector.txt","w");
		for(int i=0;i<EF_NFCOIL+current_representation[25];i++)
		{
			fprintf(Result5,"% 10.9E\n",h_S_fb[i]);
			printf("singluar value %d  % 10.9E\n",i+1,h_S_fb[i]);
		}
		fclose(Result5);

		FILE *Result6;
		Result6=fopen("u_matrix.txt","w");
		for(int i=0;i<(EF_DIAG+(int)current_representation[26])*(EF_DIAG+(int)current_representation[26]);i++)
		{
			fprintf(Result6,"% 10.9E\n",h_U_fb[i]);
		}
		fclose(Result6);

		FILE *Result7;
		Result7=fopen("vt_matrix.txt","w");
		for(int i=0;i<(EF_NFCOIL+(int)current_representation[25])*(EF_NFCOIL+(int)current_representation[25]);i++)
		{
			fprintf(Result7,"% 10.9E\n",h_VT_fb[i]);
		}
		fclose(Result7);
*////////////////////////////////

		cudaMemcpy(cu_U_matrix_fb,h_U_fb,(EF_DIAG+(int)current_representation[26])*(EF_DIAG+(int)current_representation[26])*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(cu_S_vector_fb,h_S_fb,(EF_NFCOIL+(int)current_representation[25])*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(cu_VT_matrix_fb,h_VT_fb,(EF_NFCOIL+(int)current_representation[25])*(EF_NFCOIL+(int)current_representation[25])*sizeof(double),cudaMemcpyHostToDevice);
		cudaDeviceSynchronize();

		responsematrix_inverse_U_fb<<<(EF_DIAG+(int)current_representation[26]),(EF_DIAG+(int)current_representation[26])>>>(cu_pro_diag_fb,cu_U_matrix_fb,cu_current_representation,cu_response_inverse_temp_fb);
		responsematrix_inverse_S_fb<<<1,(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp_fb,cu_S_vector_fb,cu_current_representation);	
		responsematrix_inverse_VT_fb<<<(EF_NFCOIL+(int)current_representation[25]),(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp_fb,cu_VT_matrix_fb,test_right_fb,cu_current_representation);
		
//		cudaMemcpy(h_A_fb, test_right_fb, (EF_NFCOIL+(int)current_representation[25])*sizeof(double),cudaMemcpyDeviceToHost);

		for(int i=0;i<(EF_NFCOIL+(int)current_representation[25]);i++)
		{
			printf("test_right_fb %d = % 10.9E\n",h_A_fb[i]);
		}

		printf("P-EFIT fixed boundary init done! \n");	
	}

	if(OFFLINE_TEST == 0)
		shot_set();
    	return 0;
}

int pefit_run(int m, int s_key)
{
	small_loop=0;
	while(small_loop<LOOP_CONTROL)
	{
		if (current_representation[0] == 1 || current_representation[0] == 2)
		{	
			boundflux<<<1116,1024>>>(cu_boundresponse,effect_flux,cuLbuffer);
			post_boundflux<<<35,1024>>>(cuLbuffer,cubuffer);
			final_boundflux<<<1,1024>>>(cubuffer,cuBoundResult);
			data_regroup<<<mesh1,mesh1>>>(effect_flux,cu_bfft,cu_num_se);
			gs_solver_init<<<mesh,mesh>>>(cuBoundResult,cu_bfft);
			fft_inver<<<grid,threads>>>(cu_bfft,cu_mat,cu_result+m*mesh1*mesh1);
			tri_solver<<<grid2,threads2>>>(cu_result+m*mesh1*mesh1,cu_act_trace,cu_tri_main,cu_bfft);
			fft_inver<<<grid,threads>>>(cu_mat,cu_bfft,cu_result+m*mesh1*mesh1);
		
/*/////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result12= (float *)malloc(mesh1*mesh1*sizeof(float));
			cudaMemcpy(result12,cu_result+m*mesh1*mesh1, mesh1*mesh1*sizeof(float),cudaMemcpyDeviceToHost);
	   		FILE *Result12;
			Result12=fopen("psi_ip.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh;i++)
				for(int j=0;j<mesh;j++)
					fprintf(Result12,"%15.7e\n",result12[i*mesh1+j]);
			fclose(Result12);
*///////////////////////////////////////////////////////////////////////////////////////////////////////

			if(PF_RECONSTRUCT == 1)
			{
				pfFlux<<<grid7,threads7>>>(test_right,cu_icresponse,cu_icflux);
/*/////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result11= (float *)malloc(mesh2*mesh2*sizeof(float));
			cudaMemcpy(result11, cu_icflux, mesh2*mesh2*sizeof(float),cudaMemcpyDeviceToHost);
	   		FILE *Result11;
			Result11=fopen("psi_pf.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh2;i++)
				for(int j=0;j<mesh2;j++)
					fprintf(Result11,"%15.7e\n",result11[i*mesh2+j]);
			fclose(Result11);
*///////////////////////////////////////////////////////////////////////////////////////////////////////
				ecFlux<<<grid17,threads17>>>(test_diagnostic2,cu_ecgridresponse,cu_ecflux);
				Fluxadd<<<257,257>>>(cu_icflux,cu_ecflux,cu_result+m*mesh1*mesh1);

/*/////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result14= (float *)malloc(mesh2*mesh2*sizeof(float));
			cudaMemcpy(result14, cu_ecflux, mesh2*mesh2*sizeof(float),cudaMemcpyDeviceToHost);
	   		FILE *Result14;
			Result14=fopen("psi_ec.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh2;i++)
				for(int j=0;j<mesh2;j++)
					fprintf(Result14,"%15.7e\n",result14[i*mesh2+j]);
			fclose(Result14);
*///////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result13= (float *)malloc(mesh1*mesh1*sizeof(float));
			cudaMemcpy(result13,cu_result+m*mesh1*mesh1, mesh1*mesh1*sizeof(float),cudaMemcpyDeviceToHost);
	   		FILE *Result13;
			Result13=fopen("psi_all.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh;i++)
				for(int j=0;j<mesh;j++)
					fprintf(Result13,"%15.7e\n",result13[i*mesh1+j]);
			fclose(Result13);
///////////////////////////////////////////////////////////////////////////////////////////////////////

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result12= (float *)malloc(ELI_NUM*31*15*6*sizeof(float));
			cudaMemcpy(result12,b_green_mid, ELI_NUM*31*15*6*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result12;
			Result12=fopen("bgreenmid.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<ELI_NUM*31*15*6;i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result12,"%15.7e\n",result12[i]);
			fclose(Result12);
*////////////////////////////////////////////////////////////////////////////////////////////////////////
			
				expansion_finder<<<366,63>>>(cu_result+m*mesh1*mesh1,leas_flux);
				norm1<<<3,122>>>(leas_flux);;
				xpoint_brbz<<<dim3(2,2,3),dim3(2,2,1)>>>(leas_flux,cu_result+m*mesh1*mesh1);
				post_xpoint_brbz<<<3,6>>>(leas_flux);
				xpoint_pos<<<3,1>>>(leas_flux, cu_result+m*mesh1*mesh1, xpoint+m*12,cu_xlocation);
				limiterFlux<<<1,LIMITER_NUM>>>(cu_limiter,cu_result+m*mesh1*mesh1,xpoint+m*12,cuLimFlux);
				norm3<<<1,LIMITER_NUM>>>(xpoint+m*12,cuLimFlux);

//				printf(" P-EFIT 01 IS OK \n");
				if (current_representation[0] == 1)
					data_eli<<<mesh,mesh>>>(cu_result+m*mesh1*mesh1,effect_flux,xpoint+m*12,cuLimFlux,test_right,cu_bound_flux+m,cu_num_se,cu_current_representation);
				if (current_representation[0] == 2)
					data_eli_spline<<<mesh,mesh>>>(cu_result+m*mesh1*mesh1,effect_flux,xpoint+m*12,cuLimFlux,test_right,cu_bound_flux+m,cu_num_se,cu_current_representation);

				compute_converge<<<mesh,mesh>>>(cu_result+m*mesh1*mesh1,cu_result_compare);
				compute_converge_post<<<1,128>>>(cu_result_compare,xpoint+m*12,cu_bound_flux+m,cu_converge,cu_converge_error);

//			if(OFFLINE_TEST == 0)
//				c_error_comput<<<num_seg,1>>>(cu_inter_arg+s_key*4*17, cu_inter_pos+s_key*2*17,cu_result+m*4096,cu_bound_flux+m,cu_out_buf);
//			searchbdry<<<numbdry,32>>>(cu_result+m*4096,cu_xlocation,cu_bound_flux+m,cu_tani,cu_coti,xybdry);
//			forbetapli<<<63,63>>>(cu_result+m*4096,leas_flux,test_right,cu_bound_flux+m,xpoint,cu_xlocation);
//			postbetapli<<<1,numbdry>>>(cu_result+m*4096,leas_flux,xybdry);
/*///////////////////////////////////////////////////////////////////////////////////////////////////
                        float* result33= (float *)malloc(numbdry*2*sizeof(float));
                        cudaMemcpy(result33,xybdry, numbdry*2*sizeof(float),cudaMemcpyDeviceToHost);
                        FILE *Result33;
                        Result33=fopen("xybdry.txt","w");
                //      for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
                        for(int i=0;i<numbdry*2;i++)
                                 fprintf(Result33,"% 10.9E\n",result33[i]);
                        fclose(Result33);
*///////////////////////////////////////////////////////////////////////////////////////////////////////

				tridaig_mul<<<grid3,threads3>>>(effect_flux, cu_efDiagresponse, leas_flux, cu_current_representation);
				mul_postprocess<<<(int)current_representation[25],EF_DIAG>>>(leas_flux, responsematrix,cu_current_representation);

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result2= (float *)malloc(EF_DIAG*22*sizeof(float));
			cudaMemcpy(result2,cu_gfc, EF_DIAG*22*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result2;
			Result2=fopen("responsematrix.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<EF_DIAG*22;i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result2,"%15.7e\n",result2[i]);
			fclose(Result2);
*////////////////////////////////////////////////////////////////////////////////////////////////////////
//			data_process_ec<<<EF_DIAG,EF_ECCOIL>>>(cu_ec_g,test_diagnostic2,test_diagnostic1);
//			data_process_right<<<EF_DIAG,1>>>(fwt,test_diagnostic1,cu_diagbit1,test_ec,cu_weight,cu_pro_diag);

				data_process_left<<<EF_NFCOIL+current_representation[25],EF_DIAG>>>(fwt,cu_gfc,cu_dgfc,cu_current_representation);
//			rightHandSideInit<<<EF_NFCOIL+current_representation[25],EF_DIAG+EF_NFCOIL+current_representation[26]>>>(cu_dgfc,cu_pro_diag,right_side,cu_current_representation);
//			leftHandSideInit<<<dim3(EF_NFCOIL+current_representation[25],EF_NFCOIL+current_representation[25]),threads10>>>(cu_dgfc,temp_dmate,cu_current_representation);
//			leftInit<<<(EF_NFCOIL+current_representation[25])*(EF_NFCOIL+current_representation[25]),32>>>(temp_dmate,dmata);
			
//			truncate<<<UNKNOW_ALL,1>>>(dmata);

				if (current_representation[0] == 2)
				{
					if (JPHI_FIT == 1)
					{
						jphi_constrain_1<<<numbdry,128>>>(cu_result+m*mesh1*mesh1,cu_xlocation,cu_tani,cu_coti,xpoint+12*m,cu_num_jphi,cu_jphi_npsi,cu_jphi_surface);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result12= (float *)malloc(mjphi[0]*numbdry*2*sizeof(float));
			cudaMemcpy(result12, cu_jphi_surface, mjphi[0]*numbdry*2*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result12;
			Result12=fopen("jphi_surface.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mjphi[0]*numbdry*2;i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result12,"%15.7e\n",result12[i]);
			fclose(Result12);
////////////////////////////////////////////////////////////////////////////////////////////////////////
						jphi_constrain_2<<<mjphi[0],numbdry>>>(cu_result+m*mesh1*mesh1,cu_jphi_npsi,cu_jphi_surface,cu_current_representation,cu_spline_psi,cu_jphi_left);
////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result22= (float *)malloc(mjphi[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result22, cu_jphi_left, mjphi[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result22;
			Result22=fopen("jphi_left.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mjphi[0]*(EF_NFCOIL+(int)current_representation[25]);i++)
			//	for(int j=0;j<63;j++)
				fprintf(Result22,"%15.7e\n",result22[i]);
			fclose(Result22);
////////////////////////////////////////////////////////////////////////////////////////////////////////
					}
				}
				if (MSE_FIT == 1)
				{
			
				}
				if (current_representation[0] == 1)
				{
					changePos_1<<<(EF_NFCOIL+(int)current_representation[25]),(EF_DIAG+(int)current_representation[26])>>>(cu_dgfc,d_A,cu_current_representation);
				}
				if (current_representation[0] == 2)
				{
					change_Pos_2<<<(EF_NFCOIL+(int)current_representation[25]),EF_DIAG>>>(cu_dgfc,h_A_left,cu_current_representation,cu_diag_num_all,cu_pro_diag,cu_diag_right);
					if (MSE_FIT == 1)
					{
//						change_Pos_mse<<<,>>>();
					}
					if (JPHI_FIT == 1)
					{
						change_Pos_jphi<<<(EF_NFCOIL+(int)current_representation[25]),mjphi[0]>>>(cu_jphi_left,h_A_left,cu_j_phi,cu_j_phi_uncer,cu_diag_right,cu_current_representation,cu_diag_num_all,offset_dA_jphi);
					}
					if (PRES_FIT == 1)  
					{
						change_Pos_pres<<<(EF_NFCOIL+(int)current_representation[25]),mpress[0]>>>(cu_pres_left,h_A_left,cu_press,cu_press_uncer,cu_diag_right,cu_current_representation,cu_diag_num_all,offset_dA_pres,xpoint);
					}
				}
				cudaDeviceSynchronize();

////////////////////////////////////////////////////////////////////////////////////////////////////////
			float* result4 = (float *)malloc(diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float));
			cudaMemcpy(result4, h_A_left, diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
			FILE *Result4;
			Result4=fopen("responsematrix_new_run.txt","w");
			//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25]);i++)
			{
				fprintf(Result4,"%15.7e\n",result4[i]);
//				printf("Result4 = % 10.9E\n",result4[i]);
			}
			fclose(Result4);
////////////////////////////////////////////////////////////////////////////////////////////////////////


				if (current_representation[0] == 1)
					info = LAPACKE_sgesvd(LAPACK_COL_MAJOR, 'A', 'A', (EF_DIAG+(int)current_representation[26]), (EF_NFCOIL+(int)current_representation[25]), h_A, (EF_DIAG+(int)current_representation[26]), h_S, h_U, (EF_DIAG+(int)current_representation[26]), h_VT, (EF_NFCOIL+(int)current_representation[25]), superb);

				if (current_representation[0] == 2)
				{
					cudaMemcpy(h_A, h_A_left, diag_num_all[0]*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyDeviceToHost);
					info = LAPACKE_sgesvd(LAPACK_COL_MAJOR, 'A', 'A', diag_num_all[0], (EF_NFCOIL+(int)current_representation[25]), h_A, diag_num_all[0], h_S, h_U, diag_num_all[0], h_VT, (EF_NFCOIL+(int)current_representation[25]), superb);
				}
	
//			for(int i=0;i<EF_NFCOIL+current_representation[25];i++)
//				printf("singular value (%d) = % 10.9E \n",i+1,h_S[i]);
//			printf("\n\n");
				if (current_representation[0] == 1)
					cudaMemcpyAsync(d_U,h_U,(EF_DIAG+(int)current_representation[26])*(EF_DIAG+(int)current_representation[26])*sizeof(float),cudaMemcpyHostToDevice,0);
				if (current_representation[0] == 2)
					cudaMemcpyAsync(d_U,h_U,diag_num_all[0]*diag_num_all[0]*sizeof(float),cudaMemcpyHostToDevice,0);
				
				cudaMemcpyAsync(d_S,h_S,(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyHostToDevice,0);
				cudaMemcpyAsync(d_VT,h_VT,(EF_NFCOIL+(int)current_representation[25])*(EF_NFCOIL+(int)current_representation[25])*sizeof(float),cudaMemcpyHostToDevice,0);
				cudaDeviceSynchronize();
		
				if (current_representation[0] == 1)
				{
					changePos_2<<<(EF_DIAG+(int)current_representation[26]),(EF_DIAG+(int)current_representation[26])>>>(d_U,d_S,d_VT,cu_current_representation,cu_U_matrix,cu_S_vector,cu_VT_matrix);
					responsematrix_inverse_U<<<(EF_DIAG+(int)current_representation[26]),(EF_DIAG+(int)current_representation[26])>>>(cu_pro_diag,cu_U_matrix,cu_current_representation,cu_response_inverse_temp);
				}

				if (current_representation[0] == 2)
				{
					changePos_3<<<diag_num_all[0],diag_num_all[0]>>>(d_U,d_S,d_VT,cu_current_representation,cu_U_matrix,cu_S_vector,cu_VT_matrix,cu_diag_num_all);
					responsematrix_inverse_U_spline<<<diag_num_all[0],diag_num_all[0]>>>(cu_diag_right,cu_U_matrix,cu_current_representation,cu_response_inverse_temp,cu_diag_num_all);
				}				
				
        			responsematrix_inverse_S<<<1,(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp,cu_S_vector,cu_current_representation);	
				responsematrix_inverse_VT<<<(EF_NFCOIL+(int)current_representation[25]),(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp,cu_VT_matrix,test_right,cu_current_representation);		
				currentRefresh<<<grid6,dim3(32,current_representation[25])>>>(test_right,effect_flux,cu_current_representation);
			}
		}

		if (current_representation[0] == 3)
		{
			pfFlux_fb<<<grid7,threads7>>>(test_right_fb,cu_icresponse_fb,cu_icflux_fb);
			Fluxadd_fb<<<257,257>>>(cu_icflux_fb,cu_ecflux_fb,cu_ipflux_fb,cu_result_fb+m*mesh1*mesh1);

/*////////////////////////////////////////////////////////////////////////////////////////////////////
			double* result13= (double *)malloc(mesh1*mesh1*sizeof(double));
			cudaMemcpy(result13,cu_result_fb+m*mesh1*mesh1, mesh1*mesh1*sizeof(double),cudaMemcpyDeviceToHost);
	   		FILE *Result13;
			Result13=fopen("psi_all.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh;i++)
				for(int j=0;j<mesh;j++)
					fprintf(Result13,"%15.7e\n",result13[i*mesh1+j]);
			fclose(Result13);
*///////////////////////////////////////////////////////////////////////////////////////////////////////

			expansion_finder_fb<<<366,63>>>(cu_result_fb+m*mesh1*mesh1,leas_flux_fb);
			norm1_fb<<<3,122>>>(leas_flux_fb);
			xpoint_brbz_fb<<<dim3(2,2,3),dim3(2,2,1)>>>(leas_flux_fb,cu_result_fb+m*mesh1*mesh1);
			post_xpoint_brbz_fb<<<3,6>>>(leas_flux_fb);
			xpoint_pos_fb<<<3,1>>>(leas_flux_fb, cu_result_fb+m*mesh1*mesh1, xpoint_fb+m*12,cu_xlocation_fb);
			limiterFlux_fb<<<1,LIMITER_NUM>>>(cu_limiter,cu_result_fb+m*mesh1*mesh1,xpoint_fb+m*12,cuLimFlux_fb);
			norm3_fb<<<1,LIMITER_NUM>>>(xpoint_fb+m*12,cuLimFlux_fb);

			data_eli_fb<<<mesh,mesh>>>(cu_result_fb+m*mesh1*mesh1,effect_flux_fb,cu_num_se,cu_ffprime_gafile,cu_pprime_gafile,cu_r_bdry_gafile,cu_z_bdry_gafile,xpoint_fb+m*12,cuLimFlux_fb,cu_bound_flux_fb+m);
			compute_converge_fb<<<mesh,mesh>>>(cu_result_fb+m*mesh1*mesh1,cu_result_compare_fb);
			compute_converge_post_fb<<<1,128>>>(cu_result_compare_fb,xpoint_fb+m*12,cu_bound_flux_fb+m,cu_converge,cu_converge_error_fb);

			boundflux_fb<<<1116,1024>>>(cu_boundresponse_fb, effect_flux_fb, cuLbuffer_fb);
			post_boundflux_fb<<<35,1024>>>(cuLbuffer_fb,cubuffer_fb);
			final_boundflux_fb<<<1,1024>>>(cubuffer_fb,cuBoundResult_fb);
			data_regroup_fb<<<mesh1,mesh1>>>(effect_flux_fb,cu_bfft_fb,cu_num_se);
			gs_solver_init_fb<<<mesh,mesh>>>(cuBoundResult_fb,cu_bfft_fb);
			fft_inver_fb<<<grid,threads>>>(cu_bfft_fb,cu_mat_fb,cu_ipflux_fb);
			tri_solver_fb<<<grid2,threads2>>>(cu_ipflux_fb,cu_act_trace_fb,cu_tri_main_fb,cu_bfft_fb);
			fft_inver_fb<<<grid,threads>>>(cu_mat_fb,cu_bfft_fb,cu_ipflux_fb);

/*////////////////////////////////////////////////////////////////////////////////////////////////////
			double* result15= (double *)malloc(mesh1*mesh1*sizeof(double));
			cudaMemcpy(result15,cu_ipflux_fb, mesh1*mesh1*sizeof(double),cudaMemcpyDeviceToHost);
	   		FILE *Result15;
			Result15=fopen("psi_ip_1.txt","w");
		//	for(i=0;i<=((mesh+1)*(mesh+1)-1);i++)
			for(int i=0;i<mesh;i++)
				for(int j=0;j<mesh;j++)
					fprintf(Result15,"%15.7e\n",result15[i*mesh1+j]);
			fclose(Result15);
*///////////////////////////////////////////////////////////////////////////////////////////////////////


			data_process_fb<<<1,numbdry_gafile-1>>>(cu_ipflux_fb+m*mesh1*mesh1,cu_ecflux_fb,cu_icflux_fb,cu_r_bdry_gafile,cu_z_bdry_gafile,cu_pro_diag_fb,cu_fwtbdry);

			responsematrix_inverse_U_fb<<<(EF_DIAG+(int)current_representation[26]),(EF_DIAG+(int)current_representation[26])>>>(cu_pro_diag_fb,cu_U_matrix_fb,cu_current_representation,cu_response_inverse_temp_fb);
			responsematrix_inverse_S_fb<<<1,(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp_fb,cu_S_vector_fb,cu_current_representation);	
			responsematrix_inverse_VT_fb<<<(EF_NFCOIL+(int)current_representation[25]),(EF_NFCOIL+(int)current_representation[25])>>>(cu_response_inverse_temp_fb,cu_VT_matrix_fb,test_right_fb,cu_current_representation);
		}

		small_loop++;
	}
	return 0;
}
//

void *pcs_com(void* arg)
{
	if(OFFLINE_TEST == 1)
	{
//		printf("read rfm thread exit\n");
		pthread_exit(NULL);

	}
/*	pthread_mutex_lock(&mut);
	cpu_set_t mask ;
	CPU_ZERO(&mask);
	CPU_SET(0,&mask);
	int ret = 0;

	 ret = pthread_setaffinity_np(pthread_self(),sizeof(mask),(const cpu_set_t*)&mask );
	 if(ret < 0)
	{
		printf("pthread_setaffinity_np err \n");
		exit(1);
	}

//	printf("locking\n");
	float* send_buffer= (float*)malloc(100*sizeof(float));
	float* first_buffer = (float*)malloc(75*sizeof(float));
	int niter = 0;
	while(1)
	{
		if(east_rfm_read_float(OFFSET1+time_offset*sizeof(float),first_buffer,75))
    	{
    		printf("rfm read error\n");
    		exit(-1);
    	}
		cTime = first_buffer[0];
		if(cTime >= next_ctime)
		{
			s_key++;
			next_ctime = keyTimes[s_key];
		}
		printf("\nctime: %d %f\n", time_offset, cTime);
		if(cTime >= TSTART)
		{
			s_key = s_key-1;
			break;
		}
		niter++;
		if(niter>1000000)
		{
			printf("shot error");
			exit(1);
		}
	}
//	return 0;
	memcpy(host_buffer, first_buffer, 75*sizeof(float));
	cTime = first_buffer[0];
	pthread_mutex_unlock(&mut);
//		printf("unlocking\n");
	while(1)
	{
		if(east_rfm_read_float(OFFSET1+time_offset*sizeof(float),first_buffer,75))
		{
			printf("rfm read error\n");
			exit(-1);
		}
		if(cTime>TSTART+0.01)
		{
			if(east_rfm_write_float(OFFSET2,send_buffer,50))
			{
				printf("rfm write error\n");
				exit(-1);
			}
		}
		pthread_mutex_lock(&mut);
		memcpy(host_buffer, first_buffer, 75*sizeof(float));
		memcpy(send_buffer, out_buf, 50*sizeof(float));
		cTime = first_buffer[0];
		if(cTime >= next_ctime)
		{
			s_key++;
			next_ctime = keyTimes[s_key];
//			printf("ctime:%f, s_key:%d\n",cTime, s_key);
		}
		pthread_mutex_unlock(&mut);
		if(cTime>=TEND)
		{
			printf("read rfm thread exit\n");
			pthread_exit(NULL);
		}
	}
	free(first_buffer);
	free(send_buffer);
*/
	return (void*)0;
}

void *pefit_main(void* arg)
{
    	cpu_set_t mask ;
    	CPU_ZERO(&mask);
    	CPU_SET(2,&mask);
    	int ret = 0;

    	ret =pthread_setaffinity_np(pthread_self(),sizeof(mask),(const cpu_set_t*)&mask );
    	if(ret < 0)
    	{
        	printf("pthread_setaffinity_np err \n");
        	exit(1);
    	}
	int out_m = 0;
	int m = 0;
	usleep(1000);
	if(OFFLINE_TEST == 1)
	{
/////////////////////////////////////////////////////////////
		int iteration;
		slice_end = 0;
                
//		clock_t start, finish;
//              float costtime;
//              start = clock();
//		{
/////////////////////////////////////////////////////////////
		for(iteration=0;iteration<150;iteration++)
		{
			printf("iteration = %d  ",iteration);
			pefit_run(0,0);
			cudaDeviceSynchronize();
			if (converge[0] == 1)
			{
				printf("iteration = %d  ",iteration);
				break;
			}
		}
/////////////////////////////////////////////////////////////
//		}
//		iteration=i;
//		finish = clock();
//              costtime = (float)(finish - start) / CLOCKS_PER_SEC;
		
//		for(int i=0;i<UNKNOW_ALL;i++)
//				printf("%15.7e\n" , hright[i]);
//		printf("\n");

		result_output(0);

		slice_end = 1;

//		printf("grid size: %d\n",129);
//              printf("iteration: %d\n",iteration);
//              printf("time: %f\n",costtime);
/////////////////////////////////////////////////////////////
//		post_process(1);
		pthread_exit(NULL);
	}

	while(out_m<LOOP)
	{
		m = out_m/20;
		pthread_mutex_lock(&mut);
//		printf("pefit:locking\n");
		cudaMemcpy(de_buffer,host_buffer,75*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(out_buf,cu_out_buf,50*sizeof(float),cudaMemcpyDeviceToHost);
		
		time_values[m] = cTime;
		//gettimeofday(&end,NULL);
		//time_values[m] = (end.tv_sec-start.tv_sec)*1000000+end.tv_usec-start.tv_usec;
		//time_key = time_values[m]/TIME_STRIDE;
		if(cTime>=TEND)
		{
			printf("%d loops done in one second\n",out_m);
			pthread_mutex_unlock(&mut);
//		printf("pefit:unlocking\n");
			break;
		}
//		printf("pefit:abunlocking\n");
		pthread_mutex_unlock(&mut);
		pefit_run(m,s_key);
		out_m++;
	}
	for(int i=0;i<UNKNOW_ALL;i++)
		printf("%15.7e\n" , hright[i]);
	printf("\n");
//	post_process(m);
	pthread_exit(NULL);
	return (void*)0;
}

int main(int argc,char *argv[])
{
	cudaEvent_t cu_start, cu_stop;
	float cu_time;
	cudaEventCreate(&cu_start);
	cudaEventCreate(&cu_stop);
//	cudaEventRecord( cu_start, 0 );

	clock_t start_all, finish_all;
        float costtime_all;
        start_all = clock();
	{
	
	init(argc,argv);  //	
	
	printf("init is ok!\n");

	if(OFFLINE_TEST == 0)
		next_ctime = keyTimes[0];
	s_key = 0;
//	test_diagnostic1 = de_buffer + 1;

	clock_t start, finish;
        float costtime;
        start = clock();
	{
		cudaEventRecord( cu_start, 0 );
		for(no_slice=0;no_slice<islice;no_slice++)
		{
			test_diagnostic1 = cu_diagnostic1 + no_slice*EF_DIAG;
			test_diagnostic2 = cu_diagnostic2 + no_slice*EF_ECCOIL; 
			test_psiref = cu_psiref + no_slice;
			test_btor = cu_btor + no_slice;
			host_btor = btor + no_slice;
			test_diagbit = cu_diagbit1 + no_slice*EF_DIAG;
			test_weight = cu_weight + no_slice*EF_DIAG;
			diagnostic1 = diagnostic1 + EF_DIAG;

			
			pefit_init(argc,argv,no_slice);
			printf("pefit_init is ok!\n");

//			printf("\n thread create begin\n");

			if((pthread_create(&read_sig_t, NULL, pcs_com, NULL)) != 0)
				printf("thread 1 create failed\n");
			else
//				printf("thread 1 create succed\n");

			if((pthread_create(&pefit_run_t, NULL, pefit_main, NULL)) != 0)
				printf("thread 2 create failed\n");
			else
//				printf("thread 2 create succed\n");

			struct timeval start;
    			struct timeval end;

//			gettimeofday(&start,NULL);

			pthread_join(read_sig_t,NULL);
			pthread_join(pefit_run_t,NULL);
		
			while(1)
			{
				if (slice_end == 1)
				{
//					printf("%d.%d\n",int(ishot[no_slice]),int(itime[no_slice]));
					printf("%d\n",int(itime[no_slice]));
					break;
				}
			}
		}
		cudaEventRecord( cu_stop, 0 );
		cudaEventSynchronize( cu_stop );
		cudaEventElapsedTime( &cu_time, cu_start, cu_stop );
		printf("cu_time_iteration: %f\n",cu_time);
		cudaEventDestroy( cu_start );
		cudaEventDestroy( cu_stop );
	
	}
	finish = clock();
        costtime = (float)(finish - start) / CLOCKS_PER_SEC;	
	printf("time_iteration: %f\n",costtime);

	}
	finish_all = clock();
        costtime_all = (float)(finish_all - start_all) / CLOCKS_PER_SEC;	
	printf("time_all: %f\n",costtime_all);

	return 0;
}
