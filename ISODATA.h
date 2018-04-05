// ISODATA.h: interface for the CISODATA class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ISODATA_H__F4CFB88C_3155_4ADD_93C6_F271E587249B__INCLUDED_)
#define AFX_ISODATA_H__F4CFB88C_3155_4ADD_93C6_F271E587249B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include<vector>
#include<iostream>
#include "math.h"
using namespace std;
struct distance{
	double d;
	int i;
	int j;
};//Record the distance and label between clusters in merging process
class CISODATA  
{
protected:
    vector<BYTE *> m_data;//Input data
	int *m_ClassLabel;//The label of cluster after classification
	int m_n;//The dimension of data(number of bands)
	long m_N;//The number of members(pixels) in a cluster
	int m_cN;//The minimum number of members(pixels) to form a cluster
	double m_cS;//The threshold of standard deviation of a cluster
	double m_cC;//The threshold of distance between cluster centers
	int m_L;//The maximum pairs of clusters that are allowed to be merged in each interation
	long m_I;//The maximum number of iterations
    double m_M;//The coefficient for splitting
	int m_K;//The required number of clusters
    vector<vector<double> > m_z;//Cluster Center vector
    vector<double> m_de;//The vector of the average distance between all members(pixels) in a cluster and the cluster center
	double m_dm;//The average distance between all members(pixels) in a cluster and the cluster center
	vector<vector<double> > m_s; //The standard deviation vector of each cluster center
	vector<long> m_Ne;//The number of members(pixels) in each cluster
	int m_initC;//The initial number of clusters
public:
	CISODATA();
	CISODATA(vector<BYTE *> data,int n,int N,int C,int K,int cN,double cS,double cC,int L,int I,double m=0.5);
	void InitData(vector<BYTE *> data,int n,int N,int C,int K,int cN,double cS,double cC,int L,long I,double m=0.5);
	
    void InitCenter();
	int  *ISODATA();
	double Distance(int i,int j);
	int MinD(int i);
	double ClassDistance(int i,int j);
	long* NumOfEachClass();
    vector<vector<double> > ClusterCenter();
	vector<vector<double> > RMSofEachClass();
	vector<vector<double> > AllClassDistances();
	int NumOfClasses();
	double *MeanDistance();
	double TotalMeanDistance();
	//void ReCalculateCenter();
	//void CalDistance();
	virtual ~CISODATA();
    
};

#endif // !defined(AFX_ISODATA_H__F4CFB88C_3155_4ADD_93C6_F271E587249B__INCLUDED_)
