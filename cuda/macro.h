#ifndef _MACRO_H_
#define _MACRO_H_

#define mesh 255
#define mesh1 256
#define mesh2 257
#define BLOCK_SIZE 32
#define p_pm 8
#define pi 3.14159265358979324
#define ipi 0.00000003183098862
#define ELI_NUM 35712
#define EF_NUM_GRIDPTS_FULLSET_INVESSEL 35712 // if EF_GRID_257
#define EF_NSILOP 44
#define EF_MAGPRI 76
#define EF_NROGOW 1 // if USE_2_ROGOWSKIS define = 1
#define EF_DIAG 121
//#define NUM_CU_ROWS 4
#define EF_NFCOIL 18
#define POLY_NUM 5
#define UNKNOW_ALL 23
#define RTOZ 0.2822265625
#define PF_RECONSTRUCT 1 //if reconstruct PF current define = 1
#define LIMITER_NUM 115
#define LOOP 1000000
#define INVR 75.29411764705882    //  1/(2*delR)
#define INVZ 40.0                  //   1/(2*delZ)
#define INVR2 75.29411764705882    // has not been modefied for 257
#define INVAREA 1219.047619
#define ELI_CORE 1116
#define ELI_REST 0
#define MAX_ITERATION 80
#define TIME_STRIDE 1000
#define LOOP_CONTROL 1
#define OFFSET1 0
#define OFFSET2 0X5000
#define TSTART 1
#define TEND 8.5
#define OFFLINE_TEST 1
#define EF_ECCOIL 6
#define delR 0.006640625
#define delZ 0.0125
#define igridr 0.84
#define igridz -1.6
#define egridr 2.54
#define nwnh 66049
#define uncertain 0.03
#define truncation 0.000000000001
#define capa 0.02
#define cafa 0.0025
#define numbdry 96
#define tmu 10000000
#define numlim 87
#define darea 0.0000830078125
#define truncate_factor 0.000001
#define MSE_FIT 0
#define PRES_FIT 1
#define JPHI_FIT 1
#define tension 0.5
#endif
