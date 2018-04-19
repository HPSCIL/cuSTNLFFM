

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "gdal.h"
#include "gdal_priv.h"
#include "gdalwarper.h"
#include <stdio.h>
#include<iostream>
#include "fusion.h"
#define num_thread 256
#define num_block 128
//__device__  void MultiplymBynandnByp(float *Mult1, float *Mult2, float *output, int m, int n, int p,int Idx,int Win_size1)
//{
//    int i,j,k;
//	for(i=0;i<m;i++)
//		for(j=0;j<p;j++)
//		{
//			output[i*p+j]=0;
//			for(k=0;k<n;k++)
//				output[i*p+j]+=Mult1[i*n+k+Idx*2*Win_size1*Win_size1]*Mult2[j+k*p+Idx*2*Win_size1*Win_size1];
//		}
//}
//__device__  void MatrixTransposemByn(float *input ,int m, int n, float *output,int Idx,int Win_size1)
//{
//    for(int i=0;i<n;i++)
//	for(int j=0;j<m;j++)
//	{
//		output[i*m+j+Idx*2*Win_size1*Win_size1] = input[j*n+i+Idx*2*Win_size1*Win_size1];
//	}
//}
//__device__ void CalcuRela(int num, float &a1,float &b1 ,float gamma,float *AA,float *BB,float *temp1,int Idx,int Win_size1,float B[2],float B1[2],float temp2[4],float temp3[2],float temp11[2],float temp22[4],float temp33[4],float tempB1[2],float tempB11[2], float output[4])
// {
//	 
//	 int i,j,k;
//	 a1=0,b1=0;
//	 for( i=0;i<2;i++)
//	 {
//		 for( j=0;j<num;j++)
//		 {
//			 temp1[i*num+j+Idx*2*Win_size1*Win_size1] = AA[j*2+i+Idx*2*Win_size1*Win_size1];
//		 }
//	 }
//	 for(i=0;i<2;i++)
//	 {
//		 for(j=0;j<2;j++)
//		 {
//			 temp2[i*2+j]=0;
//			 for(k=0;k<num;k++)
//				 temp2[i*2+j]+=temp1[i*num+k+Idx*2*Win_size1*Win_size1]*AA[j+k*2+Idx*2*Win_size1*Win_size1];
//		 }
//	 }
//	 // MatrixTransposemByn(AA,num,2,temp1,Idx,Win_size1);//(x,1)'
//	 // MultiplymBynandnByp(temp1,AA,temp2,2,num,2,Idx,Win_size1);//(x,1)'(x,1)
//	 for(i=0;i<2;i++)
//	 {
//		 for( j=0;j<1;j++)
//		 {
//			 temp11[i*1+j] = B1[j*2+i];
//		 }
//	 }
//	 //	 MultiplymBynandnByp(temp11,B,temp22,2,1,2,Idx,Win_size1);
//	 for(i=0;i<2;i++)
//	 {
//		 for(j=0;j<2;j++)
//		 {
//			 temp22[i*2+j]=0;
//			 for(k=0;k<1;k++)
//				 temp22[i*2+j]+=temp11[i*1+k]*B[j+k*2];
//		 }
//	 }
//	 for(i=0;i<2;i++)
//	 {
//		 for(j=0;j<2;j++)
//		 {
//			 temp33[i*2+j]=temp2[i*2+j]+temp22[i*2+j];
//		 }
//	 }
//
//	 for( i=0;i<4;i++)
//	 {
//		 output[i]=temp33[i];
//	 }
//	 output[2]=temp33[0]*temp33[3]-temp33[1]*temp33[2];
//	 output[0]=temp33[3]/output[2];
//	 output[3]=temp33[0]/output[2];
//	 output[1]=-temp33[2]/output[2];
//	 output[2]=-temp33[1]/output[2];
//	 for( i=0;i<4;i++)
//	 {
//		 temp33[i]=output[i];
//	 }
//	 // MultiplymBynandnByp(temp1,input2,temp3,2,num,1);
//	 for(i=0;i<2;i++)
//	 {
//		 for(j=0;j<1;j++)
//		 {
//			 temp3[i*1+j]=0;
//			 for(k=0;k<num;k++)
//				 temp3[i*1+j]+=temp1[i*num+k+Idx*2*Win_size1*Win_size1]*BB[j+k*1+Idx*Win_size1*Win_size1];
//		 }
//	 }
//	 //MatrixsumByp(temp3,temp11,tempB11,2,1);
//	 for(i=0;i<1;i++)
//	 {
//		 for(j=0;j<2;j++)
//		 {
//			 tempB11[i*2+j]=temp3[i*2+j]+temp11[i*2+j];
//		 }
//	 }
//	 //	 MultiplymBynandnByp(temp33,tempB11,Relationship,2,2,1);
//	 for(i=0;i<2;i++)
//	 {
//		 for(j=0;j<1;j++)
//		 {
//			 output[i+j]=0;
//			 for(k=0;k<2;k++)
//				 output[i+j]+=temp33[i*2+k]*tempB11[j+k];
//		 }
//	 }
//	 a1=output[0];
//	 b1=output[1];
// }
//__device__ void Location_pp(int *Location_P,float **BufferIn11,float **BufferIn22,float **BufferIn55,int i1,int j1,int rmin,int rmax,int smin,int smax,int nExWidth,int Height ,int Win_size1,int b,int BandNum,int &n)
//{
//	double threshold_d[10];
//	for (int i=0;i<BandNum;i++)
//	{
//		threshold_d[i]=0.01*pow(2,BufferIn11[i][ i1+nExWidth*j1]);
//	}
//
//	for(int r1=rmin;r1<=rmax;r1++)
//	{
//		for(int s1=smin;s1<=smax;s1++)
//		{  
//
//			int Result1=0;	
//			for(int i=0;i<BandNum;i++)
//			{  
//				if(fabs(BufferIn11[i][ r1+nExWidth*s1]-BufferIn11[i][ i1+nExWidth*j1])<=threshold_d[i])//ɸѡ
//				{
//					Result1++;
//				}
//				else
//					break;
//			}	
//
//			if(Result1==BandNum )
//			{	
//
//				double T1=fabs(BufferIn55[b][r1+nExWidth*s1]-BufferIn22[b][r1+nExWidth*s1]);
//				double S1=fabs(BufferIn11[b][r1+nExWidth*s1]-BufferIn22[b][r1+nExWidth*s1]);
//				if(S1<fabs(BufferIn11[b][i1+nExWidth*j1]-BufferIn22[b][i1+nExWidth*j1])+0.005||fabs(T1-fabs(BufferIn55[b][i1+nExWidth*j1]-BufferIn22[b][i1+nExWidth*j1]))<0.005)
//
//				{
//					Location_P[n+(b*nExWidth*Height+i1+nExWidth*j1)*100]=r1+nExWidth*s1;
//					n++;
//				}
//			}
//		}
//	}
//}
//__global__ void blending(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, float **BufferOut,int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int *location_p)
//{
//	const int tid=threadIdx.x;
//	const int bid=blockIdx.x;
//	const int Idx=num_thread*bid+tid;
//	int Result1=0,m=0;
//	int b=Idx/(Height*Width);
//	int j=(Idx-(Height*Width)*b)/Width;
//	int i=(Idx-(Height*Width)*b)%Width;
//	int rmin,rmax,smin,smax;
//	int n1=0;
//	float T1=0,S1=0;
//	if(b<BandNum)
//	{
//		if(i-Win_size1/2<=0)
//			rmin=0;
//		else
//			rmin = i-Win_size1/2;
//
//		if(i+Win_size1/2>=Width-1)
//			rmax = Width-1;
//		else
//			rmax = i+Win_size1/2;
//
//		if(j-Win_size1/2<=0)
//			smin=0;
//		else
//			smin = j-Win_size1/2;
//
//		if(j+Win_size1/2>=Height-1)
//			smax = Height-1;
//		else
//			smax = j+Win_size1/2;
//		//Location_pp(location_p,BufferIn11,BufferIn22,BufferIn55,i, j, rmin, rmax, smin, smax, Width,Height,Win_size1, b, BandNum,n1);
//		float threshold_d[10];
//		int r1=rmin,s1=smin;
//		for ( m=0;m<BandNum;m++)
//		{
//			threshold_d[m]=0.01*pow((float)2.0,BufferIn11[m][ i+Width*j]);
//		}
//
//		for(r1=rmin;r1<=rmax;r1++)
//		{
//			for(s1=smin;s1<=smax;s1++)
//			{  
//
//				Result1=0;	
//				for( m=0;m<BandNum;m++)
//				{  
//					if(fabs(BufferIn11[m][ r1+Width*s1]-BufferIn11[m][ i+Width*j])<=threshold_d[m])//ɸѡ
//					{
//						Result1++;
//					}
//					else
//						break;
//				}	
//
//				if(Result1==BandNum )
//				{	
//
//					T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
//					S1=fabs(BufferIn11[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
//					if(S1<fabs(BufferIn11[b][i+Width*j]-BufferIn22[b][i+Width*j])+0.005||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn22[b][i+Width*j]))<0.005)
//
//					{
//						location_p[n1+(b*Width*Height+i+Width*j)*100]=r1+Width*s1;
//						n1++;
//					}
//				}
//			}
//		}
//		float weight1=0;
//		float Lst1=0;
//		float Total_w1=0;   
//		float Average1=0;
//		int k_h=0;
//		float S=0;
//		for( k_h=0;k_h<n1;k_h++)
//		{
//
//			Lst1=BufferIn11[b][location_p[k_h+(b*Width*Height+i+Width*j)*100]]+BufferIn33[b][location_p[k_h+(b*Width*Height+i+Width*j)*100]]-BufferIn22[b][location_p[k_h+(b*Width*Height+i+Width*j)*100]];
//			S=fabs(BufferIn22[b][location_p[k_h+(b*Width*Height+i+Width*j)*100]]-BufferIn33[b][i+Width*j]);
//			weight1=exp(-(S)/h1);
//			if(Lst1>0)
//			{
//				Average1+=weight1*Lst1;
//				Total_w1+=weight1;
//			}
//			else
//			{
//				Average1+=0;
//				Total_w1+=0;
//			}
//
//		}
//		BufferOut[b][j*Width+i]=Average1/Total_w1;
//		//BufferOut[b][j*Width+i]=0;
//		if(BufferOut[b][j*Width+i]<=0)
//		{
//			BufferOut[b][j*Width+i]=BufferIn33[b][j*Width+i];
//		}
//	}
//}
//__global__ void blending1(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, float **BufferOut,int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int *location_p,float *Changed_BufferIn11,float *Changed_BufferIn33)
//{
//	const int tid=threadIdx.x;
//	const int bid=blockIdx.x;
//	const int Idx=num_thread*bid+tid;
//	int Result1=0,m=0;
//	int b=Idx/(Height*Width);
//	int j=(Idx-(Height*Width)*b)/Width;
//	int i=(Idx-(Height*Width)*b)%Width;
//	int rmin,rmax,smin,smax;
//	int n1=0;
//	float T1=0,S1=0;
//	float threshold_d[10];
//	float weight1=0;
//	float Lst1=0;
//	float Total_w1=0;
//	float Total_w2=0;
//	float Average2=0;
//	float Average1=0;
//	int k_h=0;
//	float S=0;
//	int r1,s1;
//	float Lst2;
//	float weight2;
//	double Aver11;
//	double Aver22;
//	float T1_weight;
//	float T2_weight;
//	float revise1,revise2;
//	float revise_w1,revise_w2;
//	for(int kkk=Idx;kkk<Width*Height*BandNum;kkk=kkk+num_thread*num_block)
//	{
//		Result1=0,m=0;
//		b=kkk/(Height*Width);
//		j=(kkk-(Height*Width)*b)/Width;
//		i=(kkk-(Height*Width)*b)%Width;
//		n1=0;
//		revise1=0,revise2=0,revise_w1=0,revise_w2;
//		T1=0,S1=0;
//		weight1=0;
//		weight2=0;
//		Lst1=0;
//		Lst2=0;
//		Total_w1=0;  
//		Total_w2=0;
//		Average1=0;
//		Average2=0;
//		k_h=0;
//		S=0;
//		Aver11=0;
//		Aver22=0;
//		if(b<BandNum)
//		{
//			if(i-Win_size1/2<=0)
//				rmin=0;
//			else
//				rmin = i-Win_size1/2;
//
//			if(i+Win_size1/2>=Width-1)
//				rmax = Width-1;
//			else
//				rmax = i+Win_size1/2;
//
//			if(j-Win_size1/2<=0)
//				smin=0;
//			else
//				smin = j-Win_size1/2;
//
//			if(j+Win_size1/2>=Height-1)
//				smax = Height-1;
//			else
//				smax = j+Win_size1/2;
//			//Location_pp(location_p,BufferIn11,BufferIn22,BufferIn55,i, j, rmin, rmax, smin, smax, Width,Height,Win_size1, b, BandNum,n1);
//			r1=rmin,s1=smin;
//			for ( m=0;m<BandNum;m++)
//			{
//				threshold_d[m]=0.01*pow(2,BufferIn11[m][ i+Width*j]);
//			}
//
//			for(r1=rmin;r1<=rmax;r1++)
//			{
//				for(s1=smin;s1<=smax;s1++)
//				{  
//
//					Result1=0;	
//					for( m=0;m<BandNum;m++)
//					{  
//						if(fabs(BufferIn11[m][ r1+Width*s1]-BufferIn11[m][ i+Width*j])<=threshold_d[m])//ɸѡ
//						{
//							Result1++;
//						}
//						else
//							break;
//					}	
//
//					if(Result1==BandNum )
//					{	
//
//						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
//						S1=fabs(BufferIn11[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
//						if(S1<fabs(BufferIn11[b][i+Width*j]-BufferIn22[b][i+Width*j])+0.005||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn22[b][i+Width*j]))<0.005)
//
//						{
//							location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
//							n1++;
//						}
//					}
//				}
//			}
//			for( k_h=0;k_h<n1;k_h++)
//			{
//
//				Lst1=BufferIn11[b][location_p[k_h+Idx*Win_size1*Win_size1]]+BufferIn55[b][location_p[k_h+Idx*Win_size1*Win_size1]]-BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]];
//				S=fabs(BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]]-BufferIn55[b][i+Width*j]);
//				weight1=exp(-(S)/h1);
//				revise1+= weight1*BufferIn11[b][location_p[k_h+Idx*Win_size1*Win_size1]];
//				revise_w1+= weight1;
//				if(Lst1>0)
//				{
//					Average1+=weight1*Lst1;
//					Total_w1+=weight1;
//				}
//				else
//				{
//					Average1+=0;
//					Total_w1+=0;
//				}
//
//			}
//
//
//
//
//
//
//			for ( m=0;m<BandNum;m++)
//			{
//				threshold_d[m]=0.01*pow(2,BufferIn33[m][ i+Width*j]);
//			}
//			n1=0;
//			for(r1=rmin;r1<=rmax;r1++)
//			{
//				for(s1=smin;s1<=smax;s1++)
//				{  
//
//					Result1=0;	
//					for( m=0;m<BandNum;m++)
//					{  
//						if(fabs(BufferIn33[m][ r1+Width*s1]-BufferIn33[m][ i+Width*j])<=threshold_d[m])//ɸѡ
//						{
//							Result1++;
//						}
//						else
//							break;
//					}	
//
//					if(Result1==BandNum )
//					{	
//
//						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn44[b][r1+Width*s1]);
//						S1=fabs(BufferIn33[b][r1+Width*s1]-BufferIn44[b][r1+Width*s1]);
//						if(S1<fabs(BufferIn33[b][i+Width*j]-BufferIn44[b][i+Width*j])+0.005||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn44[b][i+Width*j]))<0.005)
//
//						{
//							location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
//							n1++;
//						}
//					}
//				}
//			}
//			for( k_h=0;k_h<n1;k_h++)
//			{
//
//				Lst2=BufferIn33[b][location_p[k_h+Idx*Win_size1*Win_size1]]+BufferIn55[b][location_p[k_h+Idx*Win_size1*Win_size1]]-BufferIn44[b][location_p[k_h+Idx*Win_size1*Win_size1]];
//				S=fabs(BufferIn44[b][location_p[k_h+Idx*Win_size1*Win_size1]]-BufferIn55[b][i+Width*j]);
//				weight2=exp(-(S)/h1);
//				revise2+= weight2*BufferIn33[b][location_p[k_h+Idx*Win_size1*Win_size1]];
//				revise_w2+= weight2;
//				if(Lst2>0)
//				{
//					Average2+=weight2*Lst2;
//					Total_w2+=weight2;
//				}
//				else
//				{
//					Average2+=0;
//					Total_w2+=0;
//				}
//
//			}
//			for(int r1=rmin;r1<=rmax;r1++)
//			{
//				for(int s1=smin;s1<=smax;s1++)
//				{  
//					Aver11+=fabs(BufferIn22[b][r1+Width*s1]-BufferIn55[b][r1+Width*s1]);
//					Aver22+=fabs(BufferIn44[b][r1+Width*s1]-BufferIn55[b][r1+Width*s1]);	
//				}
//			}
//			T1_weight=1/Aver11/(1/Aver11+1/Aver22);
//			T2_weight=1/Aver22/(1/Aver11+1/Aver22);	
//			BufferOut[b][j*Width+i]=T1_weight*Average1/Total_w1+T2_weight*Average2/Total_w2;
//			//BufferOut[b][j*Width+i]=0;
//			if(BufferOut[b][j*Width+i]<=0||BufferOut[b][j*Width+i]==NULL)
//			{
//				BufferOut[b][j*Width+i]=T1_weight*revise1/revise_w1+T2_weight*revise2/revise_w2;
//			}
//		}
//	}
//}
__global__ void blending2(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, float **BufferOut,int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int *location_p,float *Changed_BufferIn11,float *Changed_BufferIn33,float *GausKernel,int p_Para,float d_Para,int pacthSize )
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int Result1=0,m=0;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int rmin,rmax,smin,smax;
	int n1=0;
	float T1=0,S1=0;
	float threshold_d[10];
	float weight1=0;
	float Lst1=0;
	float Total_w1=0;
	float Total_w2=0;
	float Average2=0;
	float Average1=0;
	int k_h=0;
	float S=0;
	int r1,s1;
	float Lst2;
	float weight2;
	double Aver11;
	double Aver22;
	float T1_weight;
	float T2_weight;
	float revise1,revise2;
	float revise_w1,revise_w2;
	float value;
	float sumGS;
	for(int kkk=Idx;kkk<Width*Height*BandNum;kkk=kkk+num_thread*num_block)
	{
		Result1=0,m=0;
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		n1=0;
		revise1=0;
		revise2=0;
		revise_w1=0;
		revise_w2=0;
		T1=0;
		S1=0;
		weight1=0;
		weight2=0;
		Lst1=0;
		Lst2=0;
		Total_w1=0;  
		Total_w2=0;
		Average1=0;
		Average2=0;
		k_h=0;
		S=0;
		Aver11=0;
		Aver22=0;
		if(b<BandNum)
		{
			if(fabs(BufferIn55[b][j*Width+i])<1e-6)
			{
				BufferOut[b][j*Width+i]=0;
			}
			else
			{
			if(i-Win_size1/2<=0)
				rmin=0;
			else
				rmin = i-Win_size1/2;

			if(i+Win_size1/2>=Width-1)
				rmax = Width-1;
			else
				rmax = i+Win_size1/2;

			if(j-Win_size1/2<=0)
				smin=0;
			else
				smin = j-Win_size1/2;

			if(j+Win_size1/2>=Height-1)
				smax = Height-1;
			else
				smax = j+Win_size1/2;
			//Location_pp(location_p,BufferIn11,BufferIn22,BufferIn55,i, j, rmin, rmax, smin, smax, Width,Height,Win_size1, b, BandNum,n1);
			r1=rmin,s1=smin;
			for ( m=0;m<BandNum;m++)
			{
				if (p_Para == 1||p_Para == 2)       //prodType == 1表示为反射率产品,prodType == 2表示为指数类产品
		       {
				threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]);
				}
				else if(p_Para == 3)
				{
					threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]/1000);
				}
			}
			float thresSpecHomo = sqrt(L_err*L_err+M_err*M_err);
	               float thresTempDiff = sqrt(2.0)*M_err;
			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum;m++)
					{  
						if(fabs(BufferIn11[m][ r1+Width*s1]-BufferIn11[m][ i+Width*j])<=threshold_d[m])//ɸѡ
						{
							Result1++;
						}
						else
							break;
					}	

					if(Result1==BandNum )
					{	

						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						S1=fabs(BufferIn11[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						if(S1<fabs(BufferIn11[b][i+Width*j]-BufferIn22[b][i+Width*j])+thresSpecHomo||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn22[b][i+Width*j]))<thresTempDiff)

						{
							location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
							n1++;
						}
					}
				}
			}
			sumGS=0;
			if(pacthSize!=1)
			{
			for(m=0;m<=pacthSize/2;m++)
			{
				value=1/((2.0*(double)m+1)*(2*(double)m+1));
				for (int r=-m; r<=m; r++)
				{
					for (int c=-m; c<=m; c++)
					{
						GausKernel[Idx*pacthSize*pacthSize+(pacthSize/2+r)*pacthSize+(pacthSize/2+c)] += value;			
					}
				}
			}
			for (m=0; m<pacthSize*pacthSize; m++)
			{
				GausKernel[Idx*pacthSize*pacthSize+m] = GausKernel[Idx*pacthSize*pacthSize+m]/(pacthSize/2);
				sumGS+= GausKernel[Idx*pacthSize*pacthSize+m];
			}

			for (m=0; m<pacthSize*pacthSize; m++)
			{
				GausKernel[Idx*pacthSize*pacthSize+m] = GausKernel[Idx*pacthSize*pacthSize+m]/sumGS;
			}
			}
			for( k_h=0;k_h<n1;k_h++)
			{

				Lst1=Changed_BufferIn11[b*Height*Width+location_p[k_h+Idx*Win_size1*Win_size1]];
				if((p_Para==1&&Lst1>0&&Lst1<=1)||(p_Para==2&&Lst1>=-1&&Lst1<=1)||(p_Para==3&&Lst1!=0))
				{
				if(pacthSize==1)
				{
				S=fabs(BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]]-BufferIn55[b][i+Width*j]);
				}
				else
				{
					int tarPRow, tarPCol;
					int simPRow, simPCol;
					int tarPLoc, simPLoc, kerPLoc;
					float sumCa= 0;
					int row=location_p[k_h+Idx*Win_size1*Win_size1]/Width;
					int colume=location_p[k_h+Idx*Win_size1*Win_size1]%Width;
					for (int rowOffset=-pacthSize/2; rowOffset<=pacthSize/2; rowOffset++)
					{
						for (int colOffset=-pacthSize/2; colOffset<=pacthSize/2; colOffset++)
						{
							if (j+rowOffset < 0)
								tarPRow = -(j+rowOffset)-1;
							else if (j+rowOffset > Height-1)
								tarPRow = 2*Height-1-(j+rowOffset);
							else
								tarPRow = j+rowOffset;

							if (i+colOffset < 0)
								tarPCol = -(i+colOffset)-1;
							else if (i+colOffset > Width-1)
								tarPCol = 2*Width-1-(i+colOffset);
							else
								tarPCol = i+colOffset;


							/* 对相似像元图块中的像元的行列范围进行判断，对超出图幅行列范围的像元进行填补 */
							if (row+rowOffset < 0)
								simPRow = -(row+rowOffset)-1;
							else if (row+rowOffset > Height-1)
								simPRow = 2*Height-1-(row+rowOffset);
							else
								simPRow = row+rowOffset;

							if (colume+colOffset < 0)
								simPCol = -(colume+colOffset)-1;
							else if (colume+colOffset > Width-1)
								simPCol = 2*Width-1-(colume+colOffset);
							else
								simPCol = colume+colOffset;

							/* 计算块的相似度 */
							tarPLoc = tarPRow * Width + tarPCol;
							simPLoc = simPRow * Width + simPCol;
							kerPLoc = (pacthSize/2+rowOffset) *pacthSize +(pacthSize/2+colOffset);
							sumCa+= GausKernel[kerPLoc+Idx*pacthSize*pacthSize]*fabs(BufferIn22[b][simPLoc]-BufferIn55[b][tarPLoc]);
						}
					}
					S=sumCa;

				}
				weight1=exp(-(S)/(h1*h1));
				revise1+= weight1*BufferIn11[b][location_p[k_h+Idx*Win_size1*Win_size1]];
				revise_w1+= weight1;
					Average1+=weight1*Lst1;
					Total_w1+=weight1;

				}

			}






			for ( m=0;m<BandNum;m++)
			{
				if (p_Para == 1||p_Para == 2)       //prodType == 1表示为反射率产品,prodType == 2表示为指数类产品
		       {
				threshold_d[m]=d_Para*pow((float)2,BufferIn33[m][ i+Width*j]);
				}
				else if(p_Para == 3)
				{
					threshold_d[m]=d_Para*pow((float)2,BufferIn33[m][ i+Width*j]/1000);
				}
			}
			n1=0;
			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum;m++)
					{  
						if(fabs(BufferIn33[m][ r1+Width*s1]-BufferIn33[m][ i+Width*j])<=threshold_d[m])//ɸѡ
						{
							Result1++;
						}
						else
							break;
					}	
				
					if(Result1==BandNum )
					{	

						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn44[b][r1+Width*s1]);
						S1=fabs(BufferIn33[b][r1+Width*s1]-BufferIn44[b][r1+Width*s1]);
						if(S1<fabs(BufferIn33[b][i+Width*j]-BufferIn44[b][i+Width*j])+thresSpecHomo||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn44[b][i+Width*j]))<thresTempDiff)

						{
							location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
							n1++;
						}
					}
				}
			}
			for( k_h=0;k_h<n1;k_h++)
			{

				Lst2=Changed_BufferIn33[b*Height*Width+location_p[k_h+Idx*Win_size1*Win_size1]];
				if((p_Para==1&&Lst2>0&&Lst2<=1)||(p_Para==2&&Lst2>=-1&&Lst2<=1)||(p_Para==3&&Lst2!=0))
				{
				if(pacthSize==1)
				{
				S=fabs(BufferIn44[b][location_p[k_h+Idx*Win_size1*Win_size1]]-BufferIn55[b][i+Width*j]);
				}
				else
				{
					int tarPRow, tarPCol;
					int simPRow, simPCol;
					int tarPLoc, simPLoc, kerPLoc;
					float sumCa= 0;
					int row=location_p[k_h+Idx*Win_size1*Win_size1]/Width;
					int colume=location_p[k_h+Idx*Win_size1*Win_size1]%Width;
					for (int rowOffset=-pacthSize/2; rowOffset<=pacthSize/2; rowOffset++)
					{
						for (int colOffset=-pacthSize/2; colOffset<=pacthSize/2; colOffset++)
						{
							if (j+rowOffset < 0)
								tarPRow = -(j+rowOffset)-1;
							else if (j+rowOffset > Height-1)
								tarPRow = 2*Height-1-(j+rowOffset);
							else
								tarPRow = j+rowOffset;

							if (i+colOffset < 0)
								tarPCol = -(i+colOffset)-1;
							else if (i+colOffset > Width-1)
								tarPCol = 2*Width-1-(i+colOffset);
							else
								tarPCol = i+colOffset;


							/* 对相似像元图块中的像元的行列范围进行判断，对超出图幅行列范围的像元进行填补 */
							if (row+rowOffset < 0)
								simPRow = -(row+rowOffset)-1;
							else if (row+rowOffset > Height-1)
								simPRow = 2*Height-1-(row+rowOffset);
							else
								simPRow = row+rowOffset;

							if (colume+colOffset < 0)
								simPCol = -(colume+colOffset)-1;
							else if (colume+colOffset > Width-1)
								simPCol = 2*Width-1-(colume+colOffset);
							else
								simPCol = colume+colOffset;

							/* 计算块的相似度 */
							tarPLoc = tarPRow * Width + tarPCol;
							simPLoc = simPRow * Width + simPCol;
							kerPLoc = (pacthSize/2+rowOffset) *pacthSize +(pacthSize/2+colOffset);
							sumCa+= GausKernel[kerPLoc+Idx*pacthSize*pacthSize]*fabs(BufferIn44[b][simPLoc]-BufferIn55[b][tarPLoc]);
						}
					}
					S=sumCa;

				}
				weight2=exp(-(S)/(h1*h1));
				revise2+= weight2*BufferIn33[b][location_p[k_h+Idx*Win_size1*Win_size1]];
				revise_w2+= weight2;
					Average2+=weight2*Lst2;
					Total_w2+=weight2;
				}

			}
			for(int r1=rmin;r1<=rmax;r1++)
			{
				for(int s1=smin;s1<=smax;s1++)
				{  
					Aver11+=fabs(BufferIn22[b][r1+Width*s1]-BufferIn55[b][r1+Width*s1]);
					Aver22+=fabs(BufferIn44[b][r1+Width*s1]-BufferIn55[b][r1+Width*s1]);	
				}
			}
			T1_weight=1/Aver11/(1/Aver11+1/Aver22);
			T2_weight=1/Aver22/(1/Aver11+1/Aver22);	
			BufferOut[b][j*Width+i]=T1_weight*Average1/Total_w1+T2_weight*Average2/Total_w2;
			//BufferOut[b][j*Width+i]=0;
			if(BufferOut[b][j*Width+i]<=0||BufferOut[b][j*Width+i]==NULL)
			{
				BufferOut[b][j*Width+i]=T1_weight*revise1/revise_w1+T2_weight*revise2/revise_w2;
			}
			}
		}
	}
}
__global__ void blending3(float **BufferIn11,float **BufferIn22,float **BufferIn55, float **BufferOut,int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int *location_p,float *Changed_BufferIn11,float *GausKernel,int p_Para,float d_Para,int pacthSize )
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int Result1=0,m=0;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int rmin,rmax,smin,smax;
	int n1=0;
	float T1=0,S1=0;
	float threshold_d[10];
	float weight1=0;
	float Lst1=0;
	float Total_w1=0;
	float Average1=0;
	int k_h=0;
	float S=0;
	int r1,s1;
	float revise1;
	float revise_w1;
	 float sumGS;
	 float value;
	for(int kkk=Idx;kkk<Width*Height*BandNum;kkk=kkk+num_thread*num_block)
	{
		value=0;
		 sumGS=0;
		Result1=0,m=0;
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		n1=0;
		revise1=0;
		revise_w1=0;
		T1=0;
		S1=0;
		weight1=0;
		Lst1=0;
		Total_w1=0;  
		Average1=0;
		k_h=0;
		S=0;
		if(b<BandNum)
		{
			if(fabs(BufferIn55[b][j*Width+i])<1e-6)
			{
				BufferOut[b][j*Width+i]=0;
			}
			else
			{
			if(i-Win_size1/2<=0)
				rmin=0;
			else
				rmin = i-Win_size1/2;

			if(i+Win_size1/2>=Width-1)
				rmax = Width-1;
			else
				rmax = i+Win_size1/2;

			if(j-Win_size1/2<=0)
				smin=0;
			else
				smin = j-Win_size1/2;

			if(j+Win_size1/2>=Height-1)
				smax = Height-1;
			else
				smax = j+Win_size1/2;
			//Location_pp(location_p,BufferIn11,BufferIn22,BufferIn55,i, j, rmin, rmax, smin, smax, Width,Height,Win_size1, b, BandNum,n1);
			r1=rmin,s1=smin;
			for ( m=0;m<BandNum;m++)
			{
				if(p_Para==1||p_Para==2)
				{
				threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]);
				}
				else if(p_Para==3)
				{
					threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]/1000);
				}
			}

			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum;m++)
					{  
						if(fabs(BufferIn11[m][ r1+Width*s1]-BufferIn11[m][ i+Width*j])<=threshold_d[m])//ɸѡ
						{
							Result1++;
						}
						else
							break;
					}	

					if(Result1==BandNum )
					{	

						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						S1=fabs(BufferIn11[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						if(S1<fabs(BufferIn11[b][i+Width*j]-BufferIn22[b][i+Width*j])+sqrt(M_err*M_err+L_err*L_err)||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn22[b][i+Width*j]))<M_err*sqrt(2.0))

						{
							location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
							n1++;
						}
					}
				}
			}
			sumGS=0;
			if(pacthSize!=1)
			{
			for(m=0;m<=pacthSize/2;m++)
			{
				value=1.0/((2.0*(double)m+1)*(2*(double)m+1));
				for (int r=-m; r<=m; r++)
				{
					for (int c=-m; c<=m; c++)
					{
						GausKernel[Idx*pacthSize*pacthSize+(pacthSize/2+r)*pacthSize+(pacthSize/2+c)] += value;			
					}
				}
			}
			for (m=0; m<pacthSize*pacthSize; m++)
			{
				GausKernel[Idx*pacthSize*pacthSize+m] = GausKernel[Idx*pacthSize*pacthSize+m]/(pacthSize/2);
				sumGS+= GausKernel[Idx*pacthSize*pacthSize+m];
			}

			for (m=0; m<pacthSize*pacthSize; m++)
			{
				GausKernel[Idx*pacthSize*pacthSize+m] = GausKernel[Idx*pacthSize*pacthSize+m]/sumGS;
			}
			}
			for( k_h=0;k_h<n1;k_h++)
			{

				Lst1=Changed_BufferIn11[b*Height*Width+location_p[k_h+Idx*Win_size1*Win_size1]];
				if((p_Para==1&&Lst1>0&&Lst1<=1)||(p_Para==2&&Lst1>=-1&&Lst1<=1)||(p_Para==3&&Lst1!=0))
				{
				if(pacthSize==1)
				{
				S=fabs(BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]]-BufferIn55[b][i+Width*j]);
				}
				else
				{
					int tarPRow, tarPCol;
					int simPRow, simPCol;
					int tarPLoc, simPLoc, kerPLoc;
					float sumCa= 0;
					int row=location_p[k_h+Idx*Win_size1*Win_size1]/Width;
					int colume=location_p[k_h+Idx*Win_size1*Win_size1]% Width;
					for (int rowOffset=-pacthSize/2; rowOffset<=pacthSize/2; rowOffset++)
					{
						for (int colOffset=-pacthSize/2; colOffset<=pacthSize/2; colOffset++)
						{
							if (j+rowOffset < 0)
								tarPRow = -(j+rowOffset)-1;
							else if (j+rowOffset > Height-1)
								tarPRow = 2*Height-1-(j+rowOffset);
							else
								tarPRow = j+rowOffset;

							if (i+colOffset < 0)
								tarPCol = -(i+colOffset)-1;
							else if (i+colOffset > Width-1)
								tarPCol = 2*Width-1-(i+colOffset);
							else
								tarPCol = i+colOffset;


							/* 对相似像元图块中的像元的行列范围进行判断，对超出图幅行列范围的像元进行填补 */
							if (row+rowOffset < 0)
								simPRow = -(row+rowOffset)-1;
							else if (row+rowOffset > Height-1)
								simPRow = 2*Height-1-(row+rowOffset);
							else
								simPRow = row+rowOffset;

							if (colume+colOffset < 0)
								simPCol = -(colume+colOffset)-1;
							else if (colume+colOffset > Width-1)
								simPCol = 2*Width-1-(colume+colOffset);
							else
								simPCol = colume+colOffset;

							/* 计算块的相似度 */
							tarPLoc = tarPRow * Width + tarPCol;
							simPLoc = simPRow * Width + simPCol;
							kerPLoc = (pacthSize/2+rowOffset) *pacthSize +(pacthSize/2+colOffset);
							sumCa+= GausKernel[kerPLoc+Idx*pacthSize*pacthSize]*fabs(BufferIn22[b][simPLoc]-BufferIn55[b][tarPLoc]);
						}
					}
					S=sumCa;

				}
				weight1=exp(-(S)/(h1*h1));
				revise1+= weight1*BufferIn11[b][location_p[k_h+Idx*Win_size1*Win_size1]];
				revise_w1+= weight1;
					Average1+=weight1*Lst1;
					Total_w1+=weight1;
			}

			}



				
			BufferOut[b][j*Width+i]=Average1/Total_w1;
			//BufferOut[b][j*Width+i]=0;
			if(BufferOut[b][j*Width+i]<0||BufferOut[b][j*Width+i]==NULL)
			{
				BufferOut[b][j*Width+i]=revise1/revise_w1;
			}
		}
		}
	}
}
__global__ void limit_a_CalcuRela(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int *location_p,float *Changed_BufferIn11,float *Changed_BufferIn33,int p_Para,float d_Para)
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int Result1=0,m=0;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int rmin,rmax,smin,smax;
	int n1=0,n2=0;
	float T1=0,S1=0;
   float threshold_d[10];
	//float B1[2]={1,0};
	//float B2[2]={gamma,0};
	//float temp2[4],temp22[4],temp33[4],tempB1[2],tempB11[2];
	//float temp3[2],temp11[2];
	int k_h=0;
	//float output[4];
	int r1,s1;
		float sumx,sumy,sumxy,sumx2;
	float aa,bb,aa2,bb2;
	for(int kkk=Idx;kkk<Width*Height*BandNum;kkk=kkk+num_thread*num_block)
	{
	//	B1[0]=1,B1[1]=0;
	//	B2[0]=1,B2[1]=0;
		aa=0,bb=0;
		aa2=0,bb2=0;
		Result1=0,m=0;
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		n1=0;
		n2=0;
		T1=0,S1=0;
		k_h=0;
			sumx=0;
			sumy=0;
		     sumxy=0;
			sumx2=0;
		if(b<BandNum)
		{
			if(fabs(BufferIn11[b][j*Width+i])<1e-6&&fabs(BufferIn33[b][j*Width+i])<1e-6)
			{
				Changed_BufferIn11[kkk]=0;
			     Changed_BufferIn33[kkk]=0;
			}
			else
			{
			if(i-Win_size1/2<=0)
				rmin=0;
			else
				rmin = i-Win_size1/2;

			if(i+Win_size1/2>=Width-1)
				rmax = Width-1;
			else
				rmax = i+Win_size1/2;

			if(j-Win_size1/2<=0)
				smin=0;
			else
				smin = j-Win_size1/2;

			if(j+Win_size1/2>=Height-1)
				smax = Height-1;
			else
				smax = j+Win_size1/2;
			//Location_pp(location_p,BufferIn11,BufferIn22,BufferIn55,i, j, rmin, rmax, smin, smax, Width,Height,Win_size1, b, BandNum,n1);
			r1=rmin,s1=smin;
			for ( m=0;m<BandNum;m++)
			{
				if (p_Para == 1||p_Para == 2)       //prodType == 1表示为反射率产品,prodType == 2表示为指数类产品
		       {
				threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]);
				}
				else if(p_Para==3)
				{
					threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]/1000);
				}

			}

			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum;m++)
					{  
						if(fabs(BufferIn11[m][ r1+Width*s1]-BufferIn11[m][ i+Width*j])<=threshold_d[m])//ɸѡ
						{
							Result1++;
						}
						else
							break;
					}	

					if(Result1==BandNum )
					{	

						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						S1=fabs(BufferIn11[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						if(S1<fabs(BufferIn11[b][i+Width*j]-BufferIn22[b][i+Width*j])+sqrt(M_err*M_err+L_err*L_err)||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn22[b][i+Width*j]))<M_err*sqrt(2.0))

						{
							location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
							n1++;
						}
					}
				}
			}
			for ( m=0;m<BandNum;m++)
			{
				if (p_Para == 1||p_Para == 2)       //prodType == 1表示为反射率产品,prodType == 2表示为指数类产品
		       {
				threshold_d[m]=d_Para*pow((float)2,BufferIn33[m][ i+Width*j]);
				}
				else if(p_Para==3)
				{
					threshold_d[m]=d_Para*pow((float)2,BufferIn33[m][ i+Width*j]/1000);
				}
			}
			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum;m++)
					{  
						if(fabs(BufferIn33[m][ r1+Width*s1]-BufferIn33[m][ i+Width*j])<=threshold_d[m])//ɸѡ
						{
							Result1++;
						}
						else
							break;
					}	

					if(Result1==BandNum )
					{	

						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn44[b][r1+Width*s1]);
						S1=fabs(BufferIn33[b][r1+Width*s1]-BufferIn44[b][r1+Width*s1]);
						if(S1<fabs(BufferIn33[b][i+Width*j]-BufferIn44[b][i+Width*j])+sqrt(M_err*M_err+L_err*L_err)||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn44[b][i+Width*j]))<M_err*sqrt(2.0))
						{
							location_p[n1+n2+Idx*Win_size1*Win_size1]=r1+Width*s1;
							n2++;
						}
					}
				}
			}
			if(n1>5&&n2>5)
			{
				sumx=0;
				sumy=0;
				sumxy=0;
				sumx2=0;
				for(k_h=0;k_h<n1;k_h++)
				{
					sumxy+=BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn55[b][location_p[k_h+Idx*Win_size1*Win_size1]];
					sumx+=BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]];
					sumy+=BufferIn55[b][location_p[k_h+Idx*Win_size1*Win_size1]];
					sumx2+=BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]];
				}
				sumxy+=gamma;
				sumx2+=gamma;
				aa=(n1*sumxy-sumx*sumy)/(n1*sumx2-sumx*sumx);
				bb=sumy/n1-aa*sumx/n1; 
				sumx=0;
				sumy=0;
				sumxy=0;
				sumx2=0;
				for(k_h=0;k_h<n2;k_h++)
				{
					sumxy+=BufferIn44[b][location_p[n1+k_h+Idx*Win_size1*Win_size1]]*BufferIn55[b][location_p[n1+k_h+Idx*Win_size1*Win_size1]];
					sumx+=BufferIn44[b][location_p[n1+k_h+Idx*Win_size1*Win_size1]];
					sumy+=BufferIn55[b][location_p[n1+k_h+Idx*Win_size1*Win_size1]];
					sumx2+=BufferIn44[b][location_p[n1+k_h+Idx*Win_size1*Win_size1]]*BufferIn44[b][location_p[n1+k_h+Idx*Win_size1*Win_size1]];
				}
				sumxy+=gamma;
				sumx2+=gamma;
				aa2=(n2*sumxy-sumx*sumy)/(n2*sumx2-sumx*sumx);
				bb2=sumy/n2-aa2*sumx/n2; 
			}
			else
			{
				sumx=0;
				sumy=0;
				sumxy=0;
				sumx2=0;
				for(r1=rmin;r1<=rmax;r1++)
				{
					for( s1=smin;s1<=smax;s1++)
					{
						sumxy+=BufferIn22[b][r1+Width*s1]*BufferIn55[b][r1+Width*s1];
						sumx+=BufferIn22[b][r1+Width*s1];
						sumy+=BufferIn55[b][r1+Width*s1];
						sumx2+=BufferIn22[b][r1+Width*s1]*BufferIn22[b][r1+Width*s1];
					}
				}
				sumxy+=gamma;
				sumx2+=gamma;
				aa=((rmax-rmin+1)*(smax-smin+1)*sumxy-sumx*sumy)/((rmax-rmin+1)*(smax-smin+1)*sumx2-sumx*sumx);
				bb=sumy/((rmax-rmin+1)*(smax-smin+1))-aa*sumx/((rmax-rmin+1)*(smax-smin+1)); 
				sumx=0;
				sumy=0;
				sumxy=0;
				sumx2=0;
				for(r1=rmin;r1<=rmax;r1++)
				{
					for( s1=smin;s1<=smax;s1++)
					{
						sumxy+=BufferIn44[b][r1+Width*s1]*BufferIn55[b][r1+Width*s1];
						sumx+=BufferIn44[b][r1+Width*s1];
						sumy+=BufferIn55[b][r1+Width*s1];
						sumx2+=BufferIn44[b][r1+Width*s1]*BufferIn44[b][r1+Width*s1];
					}
				}
				sumxy+=gamma;
				sumx2+=gamma;
				aa2=((rmax-rmin+1)*(smax-smin+1)*sumxy-sumx*sumy)/((rmax-rmin+1)*(smax-smin+1)*sumx2-sumx*sumx);
				bb2=sumy/((rmax-rmin+1)*(smax-smin+1))-aa2*sumx/((rmax-rmin+1)*(smax-smin+1)) ; 
			}
			Changed_BufferIn11[kkk]=BufferIn11[b][j*Width+i]*aa+bb;
			Changed_BufferIn33[kkk]=BufferIn33[b][j*Width+i]*aa2+bb2;
			}
		}
	}
}
__global__ void No_limit_a_CalcuRela(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, int Height,int Width,int BandNum,float *Changed_BufferIn11,float *Changed_BufferIn33)
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	for(int kkk=Idx;kkk<Width*Height*BandNum;kkk=kkk+num_thread*num_block)
	{
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		Changed_BufferIn11[kkk]=BufferIn11[b][j*Width+i]+BufferIn55[b][j*Width+i]-BufferIn22[b][j*Width+i];
		Changed_BufferIn33[kkk]=BufferIn33[b][j*Width+i]+BufferIn55[b][j*Width+i]-BufferIn44[b][j*Width+i];
	}
}
__global__ void No_limit_a_CalcuRela2(float **BufferIn11,float **BufferIn22,float **BufferIn55, int Height,int Width,int BandNum,float *Changed_BufferIn11)
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	for(int kkk=Idx;kkk<Width*Height*BandNum;kkk=kkk+num_thread*num_block)
	{
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		Changed_BufferIn11[kkk]=BufferIn11[b][j*Width+i]+BufferIn55[b][j*Width+i]-BufferIn22[b][j*Width+i];
	}
}
__global__ void limit_a_CalcuRela2(float **BufferIn11,float **BufferIn22,float **BufferIn55, int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int *location_p,float *Changed_BufferIn11,int p_Para,float d_Para)
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int Result1=0,m=0;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int rmin,rmax,smin,smax;
	int n1=0;
	float T1=0,S1=0;
   float threshold_d[10];
	//float B1[2]={1,0};
	//float B2[2]={gamma,0};
	//float temp2[4],temp22[4],temp33[4],tempB1[2],tempB11[2];
	//float temp3[2],temp11[2];
	int k_h=0;
	//float output[4];
	int r1,s1;
		float sumx,sumy,sumxy,sumx2;
	float aa,bb;
	for(int kkk=Idx;kkk<Width*Height*BandNum;kkk=kkk+num_thread*num_block)
	{
	//	B1[0]=1,B1[1]=0;
	//	B2[0]=1,B2[1]=0;
		aa=0,bb=0;
		Result1=0,m=0;
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		n1=0;
		T1=0,S1=0;
		k_h=0;
			sumx=0;
			sumy=0;
		     sumxy=0;
			sumx2=0;
		if(b<BandNum)
		{
			if(fabs(BufferIn11[b][j*Width+i])<1e-6)
			{
				Changed_BufferIn11[kkk]=0;
			}
			else
			{
			if(i-Win_size1/2<=0)
				rmin=0;
			else
				rmin = i-Win_size1/2;

			if(i+Win_size1/2>=Width-1)
				rmax = Width-1;
			else
				rmax = i+Win_size1/2;

			if(j-Win_size1/2<=0)
				smin=0;
			else
				smin = j-Win_size1/2;

			if(j+Win_size1/2>=Height-1)
				smax = Height-1;
			else
				smax = j+Win_size1/2;
			//Location_pp(location_p,BufferIn11,BufferIn22,BufferIn55,i, j, rmin, rmax, smin, smax, Width,Height,Win_size1, b, BandNum,n1);
			r1=rmin,s1=smin;
			for ( m=0;m<BandNum;m++)
			{
				if(p_Para==1||p_Para==2)
				{
				threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]);
				}
				else if(p_Para==3)
				{
					threshold_d[m]=d_Para*pow((float)2,BufferIn11[m][ i+Width*j]/1000);
				}
			}

			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum;m++)
					{  
						if(fabs(BufferIn11[m][ r1+Width*s1]-BufferIn11[m][ i+Width*j])<=threshold_d[m])//ɸѡ
						{
							Result1++;
						}
						else
							break;
					}	

					if(Result1==BandNum )
					{	

						T1=fabs(BufferIn55[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						S1=fabs(BufferIn11[b][r1+Width*s1]-BufferIn22[b][r1+Width*s1]);
						if(S1<fabs(BufferIn11[b][i+Width*j]-BufferIn22[b][i+Width*j])+sqrt(M_err*M_err+L_err*L_err)||fabs(T1-fabs(BufferIn55[b][i+Width*j]-BufferIn22[b][i+Width*j]))<M_err*sqrt(2.0))

						{
							location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
							n1++;
						}
					}
				}
			}
			if(n1>5)
			{
				sumx=0;
				sumy=0;
				sumxy=0;
				sumx2=0;
				for(k_h=0;k_h<n1;k_h++)
				{
					sumxy+=BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn55[b][location_p[k_h+Idx*Win_size1*Win_size1]];
					sumx+=BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]];
					sumy+=BufferIn55[b][location_p[k_h+Idx*Win_size1*Win_size1]];
					sumx2+=BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn22[b][location_p[k_h+Idx*Win_size1*Win_size1]];
				}
				sumxy+=gamma;
				sumx2+=gamma;
				aa=(n1*sumxy-sumx*sumy)/(n1*sumx2-sumx*sumx);
				bb=sumy/n1-aa*sumx/n1;
			}
			else
			{
				sumx=0;
				sumy=0;
				sumxy=0;
				sumx2=0;
				for(r1=rmin;r1<=rmax;r1++)
				{
					for( s1=smin;s1<=smax;s1++)
					{
						sumxy+=BufferIn22[b][r1+Width*s1]*BufferIn55[b][r1+Width*s1];
						sumx+=BufferIn22[b][r1+Width*s1];
						sumy+=BufferIn55[b][r1+Width*s1];
						sumx2+=BufferIn22[b][r1+Width*s1]*BufferIn22[b][r1+Width*s1];
					}
				}
				sumxy+=gamma;
				sumx2+=gamma;
				aa=((rmax-rmin+1)*(smax-smin+1)*sumxy-sumx*sumy)/((rmax-rmin+1)*(smax-smin+1)*sumx2-sumx*sumx);
				bb=sumy/((rmax-rmin+1)*(smax-smin+1))-aa*sumx/((rmax-rmin+1)*(smax-smin+1)); 
			}
			Changed_BufferIn11[kkk]=BufferIn11[b][j*Width+i]*aa+bb;
			}
		}
	}
}
 int runtest1(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, float **BufferOut,int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int p_Para,float d_Para,int patchSize)
{
	float **dev_BufferIn11,**dev_BufferIn22,**dev_BufferIn33,**dev_BufferIn44,**dev_BufferIn55,**dev_BufferOut,*Changed_BufferIn11,*Changed_BufferIn33;
	float **a,**f,**c,**d,**e,**out;
	//float *AA,*BB;
	a = (float**)malloc(BandNum*sizeof(float*));
	f = (float**)malloc(BandNum*sizeof(float*));
	c = (float**)malloc(BandNum*sizeof(float*));
	d = (float**)malloc(BandNum*sizeof(float*));
	e = (float**)malloc(BandNum*sizeof(float*));
	out=(float**)malloc(BandNum*sizeof(float*));
	for(int b=0;b<BandNum;b++)
	{
		cudaMalloc((void**)&a[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&f[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&c[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&d[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&e[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&out[b],Height*Width*sizeof(float));

	}
	//int num_block= Height* Width*BandNum/num_thread+1;
	int *Location_P;
	float *GausKernel;
	cudaMalloc((void***)&dev_BufferIn11,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn22,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn33,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn44,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn55,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferOut,sizeof(float*)*BandNum);
	cudaMalloc((void**)&Location_P,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	cudaMalloc((void**)&GausKernel,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	//cudaMalloc((void**)&AA,sizeof(float)*2*Win_size1*Win_size1*num_block*num_thread);
	//cudaMalloc((void**)&temp1,sizeof(float)*2*Win_size1*Win_size1*num_block*num_thread);
//	cudaMalloc((void**)&BB,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	cudaMalloc((void**)&Changed_BufferIn11,sizeof(float)*BandNum*Height*Width);
	cudaMalloc((void**)&Changed_BufferIn33,sizeof(float)*BandNum*Height*Width);
	cudaMemset(GausKernel,0,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	//cudaMalloc((void**)&Location_P,sizeof(float)*100*Width*Height*BandNum);
	for(int g=0;g<BandNum;g++)
	{
		cudaMemcpy(a[g], BufferIn11[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(f[g], BufferIn22[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(c[g], BufferIn33[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(d[g], BufferIn44[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(e[g], BufferIn55[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
	}
	cudaMemcpy(dev_BufferIn11,a,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn22,f,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn33, c,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn44, d,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn55, e,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferOut,out,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	if(flag==0)
	{
		No_limit_a_CalcuRela<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn33,dev_BufferIn44,dev_BufferIn55, Height, Width,BandNum,Changed_BufferIn11,Changed_BufferIn33);
	}
	else
	{
		limit_a_CalcuRela<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn33,dev_BufferIn44,dev_BufferIn55, Height, Width,  Win_size1, flag,L_err,M_err,h1,BandNum,gamma,Location_P,Changed_BufferIn11,Changed_BufferIn33,p_Para,d_Para);
	}
	blending2<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn33,dev_BufferIn44,dev_BufferIn55,dev_BufferOut, Height, Width,  Win_size1, flag,L_err,M_err,h1,BandNum,gamma,Location_P,Changed_BufferIn11,Changed_BufferIn33,GausKernel,p_Para,d_Para,patchSize);
	for(int g=0;g<BandNum;g++)
	{
		cudaMemcpy(BufferOut[g],out[g],Height*Width*sizeof(float),cudaMemcpyDeviceToHost);
	}
	for(int g=0;g<BandNum;g++)
	{
		cudaFree(a[g]);
		cudaFree(f[g]);
		cudaFree(c[g]);
		cudaFree(d[g]);
		cudaFree(e[g]);
		cudaFree(out[g]);
	}
	cudaFree(Changed_BufferIn11);
	cudaFree(Changed_BufferIn33);
	cudaFree(Location_P);
//	cudaFree(AA);
	//cudaFree(BB);
	//cudaFree(temp1);
	return 0;
}
  int runtest2(float **BufferIn11,float **BufferIn22,float **BufferIn55, float **BufferOut,int Height,int Width, int Win_size1,int flag,double L_err,double M_err,double h1,int BandNum,double gamma,int p_Para,float d_Para,int patchSize)
{
	float **dev_BufferIn11,**dev_BufferIn22,**dev_BufferIn55,**dev_BufferOut,*Changed_BufferIn11;
	float **a,**d,**e,**out;
	//float *AA,*BB;
	a = (float**)malloc(BandNum*sizeof(float*));
	d = (float**)malloc(BandNum*sizeof(float*));
	e = (float**)malloc(BandNum*sizeof(float*));
	out=(float**)malloc(BandNum*sizeof(float*));
	for(int b=0;b<BandNum;b++)
	{
		cudaMalloc((void**)&a[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&d[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&e[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&out[b],Height*Width*sizeof(float));

	}
	//int num_block= Height* Width*BandNum/num_thread+1;
	int *Location_P;
	float *GausKernel;

	cudaMalloc((void**)&GausKernel,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	cudaMalloc((void***)&dev_BufferIn11,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn22,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn55,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferOut,sizeof(float*)*BandNum);
	cudaMalloc((void**)&Location_P,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	//cudaMalloc((void**)&AA,sizeof(float)*2*Win_size1*Win_size1*num_block*num_thread);
	//cudaMalloc((void**)&temp1,sizeof(float)*2*Win_size1*Win_size1*num_block*num_thread);
//	cudaMalloc((void**)&BB,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	cudaMalloc((void**)&Changed_BufferIn11,sizeof(float)*BandNum*Height*Width);
	//cudaMalloc((void**)&Changed_BufferIn33,sizeof(float)*BandNum*Height*Width);
	//cudaMalloc((void**)&Location_P,sizeof(float)*100*Width*Height*BandNum);
	cudaMemset(GausKernel,0,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	for(int g=0;g<BandNum;g++)
	{
		cudaMemcpy(a[g], BufferIn11[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(d[g], BufferIn22[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(e[g], BufferIn55[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
	}
	cudaMemcpy(dev_BufferIn11,a,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn22,d,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn55, e,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferOut,out,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	if(flag==0)
	{
		No_limit_a_CalcuRela2<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn55, Height, Width,BandNum,Changed_BufferIn11);
	}
	else
	{
	limit_a_CalcuRela2<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn55, Height, Width,  Win_size1, flag,L_err,M_err,h1,BandNum,gamma,Location_P,Changed_BufferIn11,p_Para,d_Para);
	}
	blending3<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn55,dev_BufferOut, Height, Width,  Win_size1, flag,L_err,M_err,h1,BandNum,gamma,Location_P,Changed_BufferIn11,GausKernel,p_Para,d_Para,patchSize);
	for(int g=0;g<BandNum;g++)
	{
		cudaMemcpy(BufferOut[g],out[g],Height*Width*sizeof(float),cudaMemcpyDeviceToHost);
	}
	for(int g=0;g<BandNum;g++)
	{
		cudaFree(a[g]);
		cudaFree(d[g]);
		cudaFree(e[g]);
		cudaFree(out[g]);
	}
	cudaFree(Changed_BufferIn11);
	//cudaFree(Changed_BufferIn33);
	cudaFree(Location_P);
//	cudaFree(AA);
	//cudaFree(BB);
	//cudaFree(temp1);
	return 0;
}
 void IDWInterpolation(float **bufferFusResult, int rowTar, int colTar, int band, int height, int width)
{
	int numValPix = 0, halfWin = 1;
	double distIndex, sumDistIndex, sumWeigPred;

	while (numValPix <3)
	{
		numValPix = 0;
		sumDistIndex = 0;
		sumWeigPred = 0;

		/* 根据待插值像元的位置信息和搜索半径，确定搜索窗口的行列范围 */
		int rowMin = rowTar-halfWin<0 ? 0 : rowTar-halfWin;
		int colMin = colTar-halfWin<0 ? 0 : colTar-halfWin;
		int rowMax = rowTar+halfWin>height-1 ? height-1: rowTar+halfWin;
		int colMax = colTar+halfWin>width-1 ? width-1: colTar+halfWin;
		
		/* 遍历窗口内所有像元，对有效像元进行计数。使用有效像元信息，通过IDW对目标像元进行插值。*/
		for (int rowValPix = rowMin; rowValPix <= rowMax; rowValPix++)
		{
			for (int colValPix = colMin; colValPix <= colMax; colValPix++)
			{
				if(!_isnan(bufferFusResult[band][rowValPix*width+colValPix]))
				{
					distIndex = 1.0/(double)((rowTar-rowValPix)*(rowTar-rowValPix)+(colTar-colValPix)*(colTar-colValPix));
					sumWeigPred += distIndex * bufferFusResult[band][rowValPix*width+colValPix];
					sumDistIndex += distIndex;
					numValPix ++;
				}
			}
		}
		
		halfWin++;    // 更改搜索窗口半径
	}
	bufferFusResult[band][rowTar*width+colTar] = sumWeigPred / sumDistIndex;

}
 void runtest_one(float **BufferIn11,float **BufferIn22,float **BufferIn55,float **BufferOut,int Height,int Width,int Win_size1,int flag,float L_err,float M_err,float Para_h,int BandNum,float gamma,int p_Para,float d_Para ,int patchSize)
{
	int maxnum;
	size_t ff,tt;
	//cudaSetDevice(0);
	cudaMemGetInfo(&ff, &tt);
	maxnum=(ff-sizeof(float)*Win_size1*Win_size1*num_block*num_thread*2)/(BandNum*sizeof(float)*5);
	int sub_height=maxnum/Width-Win_size1;
//	sub_height=100;
	int kk=0;
	int i,j;
	float **sub_BufferIn11,**sub_BufferIn22,**sub_BufferIn55,**sub_out;
	for(int heiht_all=0;heiht_all<Height;heiht_all+=sub_height)
	{
		int task_start=kk*sub_height;
		int task_end;
		if((kk+1)*sub_height-Height<=0)
			task_end=(kk+1)*sub_height-1;
		else
			task_end=Height-1; 
		int data_start,data_end;
		if(task_start-Win_size1+1<=0)
			data_start= 0;
		else
			data_start=task_start-Win_size1+1;
		if(task_end+Win_size1-1>=Height-1)
			data_end=Height-1;
		else
			data_end=task_end+Win_size1-1;
		int data_height=data_end-data_start+1;
		sub_BufferIn11=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn22=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn55=(float**)malloc(BandNum*sizeof(float*));
		sub_out=(float**)malloc(BandNum*sizeof(float*));
		for(int b=0;b<BandNum;b++)
		{
			sub_BufferIn11[b]=new float[data_height*Width];
			sub_BufferIn22[b]=new float[data_height*Width];
			sub_BufferIn55[b]=new float[data_height*Width];
			sub_out[b]=new float[data_height*Width];
		}
		int copy;
		for(int k=0;k<BandNum;k++)
		{
			copy=0;
			for( i=data_start;i<=data_end;i++)
			{
				for( j=0;j<Width;j++)
				{
					sub_BufferIn11[k][copy*Width+j]=BufferIn11[k][i*Width+j];
					sub_BufferIn22[k][copy*Width+j]=BufferIn22[k][i*Width+j];
					sub_BufferIn55[k][copy*Width+j]=BufferIn55[k][i*Width+j];
				}
				copy++;
			}
		}
		int current=task_start-data_start;
		runtest2(sub_BufferIn11,sub_BufferIn22,sub_BufferIn55,sub_out,data_height,Width,Win_size1,flag, L_err, M_err, Para_h,BandNum,1.0,p_Para,d_Para,patchSize);
		
		for(int k=0;k<BandNum;k++)
		{
			current=task_start-data_start;
			for(int i=task_start;i<=task_end;i++)
			{
				for(int j=0;j<Width;j++)
				{
					BufferOut[k][i*Width+j]=sub_out[k][current*Width+j];
				}
				current++;
			}
		}
		for(int g=0;g<BandNum;g++)
	{
		delete sub_BufferIn11[g];
		delete sub_BufferIn22[g];
		delete sub_BufferIn55[g];
		delete sub_out[g];
		/*cudaFree(dev_BufferIn11[g]);
		cudaFree(dev_BufferIn22[g]);
		cudaFree(dev_BufferIn33[g]);
		cudaFree(dev_BufferIn44[g]);
		cudaFree(dev_BufferIn55[g]);
		cudaFree(dev_BufferOut[g]);*/
	}
		kk++;
	}
	for(int b=0; b<BandNum; b++)
	{ 
		for(int j=0; j<Height; j++)
		{
			for(int i=0; i<Width;i++)
			{
				if (_isnan(BufferOut[b][j*Width+i]))
				{
					/* 对无效值位置处进行IDW插值 */
					IDWInterpolation(BufferOut, j, i, b, Height, Width);
				}
			}
		}
	}
}
 void runtest(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55,float **BufferOut,int Height,int Width,int Win_size1,int flag,float L_err,float M_err,float Para_h,int BandNum,float gamma,int p_Para,float d_Para ,int patchSize)
{
	int maxnum;
	size_t ff,tt;
	//cudaSetDevice(0);
	cudaMemGetInfo(&ff, &tt);
	maxnum=(ff-sizeof(float)*Win_size1*Win_size1*num_block*num_thread*2)/(BandNum*sizeof(float)*8);
	int sub_height=maxnum/Width-Win_size1;
//	sub_height=100;
	int kk=0;
	int i,j;
	float **sub_BufferIn11,**sub_BufferIn22,**sub_BufferIn33,**sub_BufferIn44,**sub_BufferIn55,**sub_out;
	for(int heiht_all=0;heiht_all<Height;heiht_all+=sub_height)
	{
		int task_start=kk*sub_height;
		int task_end;
		if((kk+1)*sub_height-Height<=0)
			task_end=(kk+1)*sub_height-1;
		else
			task_end=Height-1; 
		int data_start,data_end;
		if(task_start-Win_size1+1<=0)
			data_start= 0;
		else
			data_start=task_start-Win_size1+1;
		if(task_end+Win_size1-1>=Height-1)
			data_end=Height-1;
		else
			data_end=task_end+Win_size1-1;
		int data_height=data_end-data_start+1;
		sub_BufferIn11=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn22=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn33=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn44=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn55=(float**)malloc(BandNum*sizeof(float*));
		sub_out=(float**)malloc(BandNum*sizeof(float*));
		for(int b=0;b<BandNum;b++)
		{
			sub_BufferIn11[b]=new float[data_height*Width];
			sub_BufferIn22[b]=new float[data_height*Width];
			sub_BufferIn33[b]=new float[data_height*Width];
			sub_BufferIn44[b]=new float[data_height*Width];
			sub_BufferIn55[b]=new float[data_height*Width];
			sub_out[b]=new float[data_height*Width];
		}
		int copy;
		for(int k=0;k<BandNum;k++)
		{
			copy=0;
			for( i=data_start;i<=data_end;i++)
			{
				for( j=0;j<Width;j++)
				{
					sub_BufferIn11[k][copy*Width+j]=BufferIn11[k][i*Width+j];
					sub_BufferIn22[k][copy*Width+j]=BufferIn22[k][i*Width+j];
					sub_BufferIn33[k][copy*Width+j]=BufferIn33[k][i*Width+j];
					sub_BufferIn44[k][copy*Width+j]=BufferIn44[k][i*Width+j];
					sub_BufferIn55[k][copy*Width+j]=BufferIn55[k][i*Width+j];
				}
				copy++;
			}
		}
		int current=task_start-data_start;
		runtest1(sub_BufferIn11,sub_BufferIn22,sub_BufferIn33,sub_BufferIn44,sub_BufferIn55,sub_out,data_height,Width,Win_size1,flag, L_err, M_err, Para_h,BandNum,1.0,p_Para,d_Para,patchSize);
		
		for(int k=0;k<BandNum;k++)
		{
			current=task_start-data_start;
			for(int i=task_start;i<=task_end;i++)
			{
				for(int j=0;j<Width;j++)
				{
					BufferOut[k][i*Width+j]=sub_out[k][current*Width+j];
				}
				current++;
			}
		}
		for(int g=0;g<BandNum;g++)
	{
		delete sub_BufferIn11[g];
		delete sub_BufferIn22[g];
		delete sub_BufferIn33[g];
		delete sub_BufferIn44[g];
		delete sub_BufferIn55[g];
		delete sub_out[g];
		/*cudaFree(dev_BufferIn11[g]);
		cudaFree(dev_BufferIn22[g]);
		cudaFree(dev_BufferIn33[g]);
		cudaFree(dev_BufferIn44[g]);
		cudaFree(dev_BufferIn55[g]);
		cudaFree(dev_BufferOut[g]);*/
	}
		kk++;
	}
	for(int b=0; b<BandNum; b++)
	{ 
		for(int j=0; j<Height; j++)
		{
			for(int i=0; i<Width;i++)
			{
				if (_isnan(BufferOut[b][j*Width+i]))
				{
					/* 对无效值位置处进行IDW插值 */
					IDWInterpolation(BufferOut, j, i, b, Height, Width);
				}
			}
		}
	}
}
 void Re_fusion(CuLayer *psensor,PARAMETER *par)
 {
	 int c;
	long now1 = clock();
	for(c=0;c<par->NUM_PREDICTIONS;c++)
	{
		psensor[2*(par->NUM_PAIRS+c)].Read(psensor[2*(par->NUM_PAIRS+c)].outpath);
		psensor[2*(par->NUM_PAIRS+c)+1].resize(psensor[0].getWidth(),psensor[0].getHeight(),psensor[0].getbandCount());
		if(par->NUM_PAIRS==2)
		{
			runtest(psensor[0].getData(),psensor[2].getData(),psensor[1].getData(),psensor[3].getData(),psensor[2*(par->NUM_PAIRS+c)].getData(),psensor[2*(par->NUM_PAIRS+c)+1].getData(),psensor[0].getHeight(),psensor[0].getWidth(),par->WIN_SIZE,par->r,par->L_ERR,par->M_ERR,par->h,psensor[0].getbandCount(),par->gamma,par->p,par->d,par->pathSize);
			//char* driverName = "GTiff";
			psensor[2*(par->NUM_PAIRS+c)+1].setGeoTransform(psensor[0].getGeoTransform());
			psensor[2*(par->NUM_PAIRS+c)+1].setProjection(psensor[0].getProjection());
			psensor[2*(par->NUM_PAIRS+c)+1].Write(psensor[2*(par->NUM_PAIRS+c)+1].outpath,par->G_Type);
		}
		else if(par->NUM_PAIRS==1)
		{
			runtest_one(psensor[0].getData(),psensor[1].getData(),psensor[2*(par->NUM_PAIRS+c)].getData(),psensor[2*(par->NUM_PAIRS+c)+1].getData(),psensor[0].getHeight(),psensor[0].getWidth(),par->WIN_SIZE,par->r,par->L_ERR,par->M_ERR,par->h,psensor[0].getbandCount(),par->gamma,par->p,par->d,par->pathSize);
			psensor[2*(par->NUM_PAIRS+c)+1].setGeoTransform(psensor[0].getGeoTransform());
			psensor[2*(par->NUM_PAIRS+c)+1].setProjection(psensor[0].getProjection());
			psensor[2*(par->NUM_PAIRS+c)+1].Write(psensor[2*(par->NUM_PAIRS+c)+1].outpath,par->G_Type);
		}
		psensor[2*(par->NUM_PAIRS+c)].ddata();
		psensor[2*(par->NUM_PAIRS+c)+1].ddata();
		//delete &psensor[2*(par->NUM_PAIRS+c)];
		//delete &psensor[2*(par->NUM_PAIRS+c)+1];
	}
	 printf("GPU运行时间为：%dms\n", int(((double)(clock() - now1)) / CLOCKS_PER_SEC * 1000));
 }
//void Re_fusion2(const char * BufferIn0,const char * BufferIn1,const char * BufferIn2,const char * BufferIn3,const char * BufferIn4,const char * BufferOut,int win_size,int flag,float L_err,float M_err,float Para_h)
//{
//	GDALAllRegister();
//	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO"); 
//	GDALDataset *Landsat0 = (GDALDataset*) GDALOpen(BufferIn0,GA_ReadOnly);
//	int width,height,BandNum;
//	width = Landsat0->GetRasterXSize();
//	height = Landsat0->GetRasterYSize();
//	BandNum = Landsat0->GetRasterCount();
//	float** BufferLandsat_0 = new float*[BandNum];
//	int b,k;
//	for( b=0;b<BandNum;b++)
//	{
//		BufferLandsat_0[b] = new float[width*height];	
//	}
//	
//	for( k=0;k<BandNum;k++)
//	{
//		GDALRasterBand* hInBand1 = Landsat0->GetRasterBand(k+1);
//		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferLandsat_0[k],width,height,GDT_Float32,0,0);
//	}	
//
//	GDALDataset *MODIS0 = (GDALDataset*) GDALOpen(BufferIn1,GA_ReadOnly);
//	float** BufferModis_0 = new float*[BandNum];
//	for( b=0;b<BandNum;b++)
//	{
//		BufferModis_0[b] = new float[width*height];	
//	}
//	
//	for( k=0;k<BandNum;k++)
//	{
//		GDALRasterBand* hInBand1 = MODIS0->GetRasterBand(k+1);
//		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferModis_0[k],width,height,GDT_Float32,0,0);		
//	}	
//
//	GDALDataset *Landsat1 = (GDALDataset*) GDALOpen(BufferIn2,GA_ReadOnly);
//	float** BufferLandsat_1 = new float*[BandNum];
//	for( b=0;b<BandNum;b++)
//	{
//		BufferLandsat_1[b] = new float[width*height];	
//	}
//	
//	for(k=0;k<BandNum;k++)
//	{
//		GDALRasterBand* hInBand1 = Landsat1->GetRasterBand(k+1);
//		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferLandsat_1[k],width,height,GDT_Float32,0,0);
//	}	
//
//	GDALDataset *MODIS1 = (GDALDataset*) GDALOpen(BufferIn3,GA_ReadOnly);
//	float** BufferModis_1 = new float*[BandNum];
//	for( b=0;b<BandNum;b++)
//	{
//		BufferModis_1[b] = new float[width*height];	
//	}
//	
//	for( k=0;k<BandNum;k++)
//	{
//		GDALRasterBand* hInBand1 = MODIS1->GetRasterBand(k+1);
//		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferModis_1[k],width,height,GDT_Float32,0,0);		
//	}	
//	
//	GDALDataset *MODIS2 = (GDALDataset*) GDALOpen(BufferIn4,GA_ReadOnly);
//	
//	float** BufferModis_2 = new float*[BandNum];
//	for( b=0;b<BandNum;b++)
//	{
//		BufferModis_2[b] = new float[width*height];	
//	}
//	
//	for( k=0;k<BandNum;k++)
//	{
//		GDALRasterBand* hInBand1 = MODIS2->GetRasterBand(k+1);
//		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferModis_2[k],width,height,GDT_Float32,0,0);		
//		
//	}
//	
//	GDALDataset *LandsatDs;
//	char* driverName = "GTiff";
//	GDALDriver *pDriver = (GDALDriver*)GDALGetDriverByName(driverName);
//	LandsatDs = pDriver->Create(BufferOut,width,height,BandNum,GDT_Float64,NULL);
//	double* geos=new double[6];
//	Landsat0->GetGeoTransform(geos);
//	LandsatDs->SetGeoTransform(geos);
//	LandsatDs->SetProjection(Landsat0->GetProjectionRef());
//	
//	float** BufferOutColor = new float*[BandNum];
//	for( b=0;b<BandNum;b++)
//	{
//		BufferOutColor[b] = new float[width*height];
//	}
//	//e.Blending2(BufferLandsat_0,BufferModis_0,BufferLandsat_1,BufferModis_1,BufferModis_2,BufferOutColor,height,width,win_size,flag, L_err, M_err, Para_h,BandNum,1.0);
//	long now1 = clock();
//	 runtest(BufferLandsat_0,BufferModis_0,BufferLandsat_1,BufferModis_1,BufferModis_2,BufferOutColor,height,width,win_size,flag, L_err, M_err, Para_h,BandNum,1.0);
//	 printf("GPU运行时间为：%dms\n", int(((double)(clock() - now1)) / CLOCKS_PER_SEC * 1000));
//	for (b=0;b<BandNum;b++)
//	{
//		GDALRasterBand* HOut = LandsatDs->GetRasterBand(b+1);
//		HOut->RasterIO(GF_Write,0,0,width,height,BufferOutColor[b],width,height,GDT_Float32,0,0);
//	}
//	GDALClose(Landsat0);
//	GDALClose(MODIS0);
//	GDALClose(Landsat1);
//	GDALClose(MODIS1);
//	GDALClose(MODIS2);
//	GDALClose(LandsatDs);
//
//	for (b=0;b<BandNum;b++)
//	{
//		delete []BufferLandsat_0[b];
//		delete []BufferModis_0[b];
//		delete []BufferLandsat_1[b];
//		delete []BufferModis_1[b];
//		delete []BufferModis_2[b];
//		delete []BufferOutColor[b];
//	}
//	delete []BufferLandsat_0;
//	delete [] BufferModis_0;
//	delete []BufferLandsat_1;
//	delete [] BufferModis_1;
//	delete [] BufferModis_2;
//	delete [] BufferOutColor;
//}
//int main()
//{
//	const char* modFile1="D:\\cuda\\shikong\\软件\\测试数据\\M_2002_01_04.tif";
//    const char* modFile2="D:\\cuda\\shikong\\软件\\测试数据\\M_2002_02_21.tif";
//	const char* tifFile1="D:\\cuda\\shikong\\软件\\测试数据\\L_2002_01_04.tif";
//	const char* tifFile2="D:\\cuda\\shikong\\软件\\测试数据\\L_2002_02_21.tif";
//
//	const char* modFile0="D:\\cuda\\shikong\\软件\\测试数据\\M_2002_02_12.tif";
//	const char* out="D:\\cuda\\shikong\\软件\\测试数据\\kk_Llick_2001_11_02.tif";
//	Re_fusion2(tifFile1, modFile1,tifFile2,modFile2,modFile0,out,51,1,0.0028,0.0028,0.03);
//	cudaDeviceReset();
//	return 0;
//}