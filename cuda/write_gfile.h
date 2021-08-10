#include<stdio.h>
#include<math.h>
//#include<fftw3.h>
#include<malloc.h>

int write_gfile_function(int num, float* data, FILE* gfile, int factor)
{
	int n,m,l;
	FILE* result;
	int i,f;

	n=num;
	result = gfile;
	f = factor;

	m = n/5;
	l = n - m*5;

	for(i=0;i<m;i++)
	{
		fprintf(result,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",f*data[5*i],f*data[5*i+1],f*data[5*i+2],f*data[5*i+3],f*data[5*i+4]);
	}
	
	if (l != 0)
	{
		for(i=0;i<l;i++)
		{
			fprintf(result,"% 10.9E",f*data[m*5+i]);
		}
		fprintf(result,"\n");
	}

	return 0;
}

int write_gfile_function_double(int num, double* data, FILE* gfile, int factor)
{
	int n,m,l;
	FILE* result;
	int i,f;

	n=num;
	result = gfile;
	f = factor;

	m = n/5;
	l = n - m*5;

	for(i=0;i<m;i++)
	{
		fprintf(result,"% 10.9E% 10.9E% 10.9E% 10.9E% 10.9E\n",f*data[5*i],f*data[5*i+1],f*data[5*i+2],f*data[5*i+3],f*data[5*i+4]);
	}
	
	if (l != 0)
	{
		for(i=0;i<l;i++)
		{
			fprintf(result,"% 10.9E",f*data[m*5+i]);
		}
		fprintf(result,"\n");
	}

	return 0;
}

/////////////////////////////////////////////////////////////////

int write_afile_function(int num, float* data, FILE* gfile)
{
	int n,m,l;
	FILE* result;
	int i;

	n=num;
	result = gfile;

	m = n/4;
	l = n - m*4;

	for(i=0;i<m;i++)
	{
		fprintf(result," % 9.9E% 9.9E% 9.9E% 9.9E\n",data[4*i],data[4*i+1],data[4*i+2],data[4*i+3]);
	}
	
	if (l != 0)
	{
		fprintf(result," ");
		for(i=0;i<l;i++)
		{
			fprintf(result,"% 9.9E",data[m*4+i]);
		}
		fprintf(result,"\n");
	}

	return 0;
}

int write_afile_function_double(int num, double* data, FILE* gfile)
{
	int n,m,l;
	FILE* result;
	int i;

	n=num;
	result = gfile;

	m = n/4;
	l = n - m*4;

	for(i=0;i<m;i++)
	{
		fprintf(result," % 9.9E% 9.9E% 9.9E% 9.9E\n",data[4*i],data[4*i+1],data[4*i+2],data[4*i+3]);
	}
	
	if (l != 0)
	{
		fprintf(result," ");
		for(i=0;i<l;i++)
		{
			fprintf(result,"% 9.9E",data[m*4+i]);
		}
		fprintf(result,"\n");
	}

	return 0;
}

