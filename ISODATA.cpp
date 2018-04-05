// ISODATA.cpp: implementation of the CISODATA class.
//
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "ISODATA.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CISODATA::CISODATA()
{
	m_ClassLabel=NULL;
}

CISODATA::CISODATA(vector<BYTE *> data,int n,int N,int C,int K,int cN,double cS,double cC,int L,int I,double m)
{
    m_data=data;m_cC=cC;m_cN=cN;m_K=K;
	m_cS=cS;m_n=n;m_L=L;m_M=m;m_N=N;m_initC=C;
	m_ClassLabel=new int[N];
	m_z.resize(K);
	for(int i=0;i<K;i++)
		m_z[i].resize(n);
}
CISODATA::~CISODATA()
{
     delete []m_ClassLabel;
	 m_data.clear();
}
/////////////////////////////////////////////////////////////////////////////
/* Name£ºInitData
   Function£ºInitialize the parameters in ISODATA algorithm
   Return£ºvoid
*/

void CISODATA::InitData(vector<BYTE *> data,int n,int N,int C,int K,int cN,double cS,double cC,int L,long I,double m)
{
	m_data=data;m_cC=cC;m_cN=cN;m_K=K;
	m_cS=cS;m_n=n;m_L=L;m_M=m;m_N=N;m_initC=C;
	m_I=I;
	m_ClassLabel=new int[N];
	m_s.resize(m_n);
}
/////////////////////////////////////////////////////////////////////////////
/*Name£ºInitCenter
  Function£ºInitialized cluster centers
  Return£ºvoid
*/
void CISODATA::InitCenter()
{
	m_z.resize(m_n);
	for(int i=0;i<m_n;i++)
	{
		m_z[i].resize(m_initC);
		for(int j=0;j<m_initC;j++)
			m_z[i][j]=m_data[i][j];
	}
	return;
}
/////////////////////////////////////////////////////////////////////////////
/*Name£ºDistance
  Function£ºCalculate the distance between a member(pixel) in all bands and the cluster center
  Parameters£ºi is the member(pixel) label£¬j is the cluster center label
  Return£ºDistance
*/
double CISODATA::Distance(int i,int j)
{
	double d=0;
	for(int k=0;k<m_n;k++)
		d+=(m_data[k][i]-m_z[k][j])*(m_data[k][i]-m_z[k][j]);
	d=sqrt(d);
    return d;
}
////////////////////////////////////////////////////////////////////////////
/*Name£ºClassDistance
  Function£ºCalculate the distance between two cluster centers
  Parameters£ºi,j are the lables of the cluster centers
  Return£ºDistance
*/
double CISODATA::ClassDistance(int i,int j)
{
	double d=0;
	for(int k=0;k<m_n;k++)
		d+=(m_z[k][i]-m_z[k][j])*(m_z[k][i]-m_z[k][j]);
	d=sqrt(d);
	
    return d;
}
/////////////////////////////////////////////////////////////////////////////
/*Name£ºMinD
  Function£ºCalculate the label of a cluster center with the minimum distance to a member(pixel)
  Parameters£ºi is the label of member(pixel)
  Return£ºThe label of cluster center
*/
int CISODATA::MinD(int i)
{
	int min;
	double d,dmin;
	dmin=Distance(i,0);
	min=0;
	for(int j=1;j<m_z[0].size();j++)
	{
		d=Distance(i,j);
		if(d<dmin)
		{
			dmin=d;
			min=j;
			
		}
	    
	}
	
	return min;
}

int* CISODATA::ISODATA()
{
	int C=m_initC;
	int i,j,k;
	long num=0;//The number of iterations
	int flag=0;//Mark whether to stop the iteration or not
	vector<int> S;//The largest component
	vector<double> smax;
	vector<int> IsCombine;
	vector<double> dist;//The distance between two cluster centers
	vector<int> dist_i;//The label for each cluster center
	vector<int> dist_j;
	
	if(m_cN>=2*m_N/m_K) 
	{
		cout<<"Please adjust parameters cN and K!"<<endl;
		return NULL;
	}
    InitCenter();
	
    m_Ne.resize(C*2);
	m_de.resize(C*2);
	S.resize(C*2);
	smax.resize(C*2);
	IsCombine.resize(C*2);
	dist.resize(C*(C-1)/2);
	dist_i.resize(C*(C-1)/2);
	dist_j.resize(C*(C-1)/2);

    while(1)
	{
		if(C>S.size())
		{
			S.resize(S.size()*2);
			m_Ne.resize(S.size());
			m_de.resize(S.size());
			smax.resize(S.size());
			IsCombine.resize(S.size());
		}
		num++;
		////////////
		for(i=0;i<C;i++)
			m_Ne[i]=0;
		for(i=0;i<m_N;i++)
		{
			m_ClassLabel[i]=MinD(i);
			m_Ne[m_ClassLabel[i]]++;
		}
		
		int r=m_z[0].size();
		for(i=0;i<C;i++)
			if(m_Ne[i]<m_cN)
			{
				
				for(j=0;j<m_n;j++)
				    m_z[j].erase(m_z[j].begin()+i);
				C--;
				flag=1;
				break;
			}
	    if(flag==1)
		    {num--;flag=0;continue;}
		for(i=0;i<m_n;i++)
			for(j=0;j<C;j++)
				m_z[i][j]=0;
		for(i=0;i<m_n;i++)
		    for(j=0;j<m_N;j++)
			    m_z[i][m_ClassLabel[j]]+=m_data[i][j];
		for(i=0;i<m_n;i++)
		    for(j=0;j<C;j++)
				m_z[i][j]/=m_Ne[j];
		//Stop iteration
		if(num>m_I) 
		{
			//Calculate the standard deviation vector when stopping iteration
			for(i=0;i<m_n;i++)
			{
				m_s[i].resize(C);
				for(j=0;j<C;j++)
					m_s[i][j]=0;
			}
			for(k=0;k<m_N;k++)
			    for(i=0;i<m_n;i++)
					{
						m_s[i][m_ClassLabel[k]]+=(m_data[i][k]-m_z[i][m_ClassLabel[k]])*(m_data[i][k]-m_z[i][m_ClassLabel[k]]);
					}
			for(i=0;i<m_n;i++)
				for(j=0;j<C;j++)
					m_s[i][j]=sqrt(m_s[i][j]/m_Ne[j]);
			//Calculate the average distance between members(pixels) and the cluster center in each cluster when stopping interation
            for(i=0;i<C;i++)
			    m_de[i]=0;
		    for(i=0;i<m_N;i++)
			    m_de[m_ClassLabel[i]]+=Distance(i,m_ClassLabel[i]);
		    m_dm=0;
		    for(i=0;i<C;i++)
			{
			    m_dm+=m_de[i];
			    m_de[i]/=m_Ne[i];
			}
		    m_dm/=m_N;
			break;                                    
		}
        
        if(C<m_K/2||(num%2==1&&C<=m_K))//The splitting process
		{   
			//Calculate the average distance between members(pixels) and the cluster center in each cluster when stopping interation
            for(i=0;i<C;i++)
			    m_de[i]=0;
		    for(i=0;i<m_N;i++)
			    m_de[m_ClassLabel[i]]+=Distance(i,m_ClassLabel[i]);
		    m_dm=0;
		    for(i=0;i<C;i++)
			{
			    m_dm+=m_de[i];
			    m_de[i]/=m_Ne[i];
			}
		    m_dm/=m_N;
			//Calculate the stadard deviation vector of each cluster
			for(i=0;i<m_n;i++)
			{
				m_s[i].resize(C);
				for(j=0;j<C;j++)
					m_s[i][j]=0;
			}
			for(k=0;k<m_N;k++)
			    for(i=0;i<m_n;i++)
					{
						m_s[i][m_ClassLabel[k]]+=(m_data[i][k]-m_z[i][m_ClassLabel[k]])*(m_data[i][k]-m_z[i][m_ClassLabel[k]]);
					}
			for(i=0;i<m_n;i++)
				for(j=0;j<C;j++)
					m_s[i][j]=sqrt(m_s[i][j]/m_Ne[j]);
			for(j=0;j<C;j++)
			{
				smax[j]=m_s[0][j];
				S[j]=0;
				for(i=1;i<m_n;i++)
					if(m_s[i][j]>smax[j])
					{
						smax[j]=m_s[i][j];
						S[j]=i;
					}
			}
			int t;
			for(t=0,j=0;t<C;j++,t++)
			{
				if(smax[j]>m_cS)
				{
					double zhu=smax[j];
					if((m_de[j]>m_dm&&m_Ne[j]>2*(m_cN+1))||C<=m_K/2)
					{
						m_z[S[j]][t]+=m_M*smax[j];
						for(i=0;i<m_n;i++)
							m_z[i].insert(m_z[i].begin()+t,m_z[i][t]);
						m_z[S[j]][t+1]-=2*m_M*smax[j];
						C++;
						t++;
					}
				}
			}
            
		}
		else//The merging process
		{
			
			
			dist.clear();
		    dist_i.clear();
			dist_j.clear();
			
		    double D;//Record distance
			int t,k,znew;//znew is the new merging center£¬IsCombine marks whether to merge or not
			for(i=0;i<C;i++)
				IsCombine[i]=0;
			k=0;
		  
			for(i=0;i<C;i++)
				for(j=i+1;j<C;j++)
				{
					D=ClassDistance(i,j);
					if(D<m_cC)
					{
						if(k==0) 
						{
							dist.push_back(D);
							dist_i.push_back(i);
							dist_j.push_back(j);
							k++;
							
						}
					    else
						{
							for(t=0;t<k;t++)
								if(D<dist[t])
									break;
								{
									dist.insert(dist.begin()+t,D);
									dist_i.insert(dist_i.begin()+t,i);
									dist_j.insert(dist_j.begin()+t,j);
									k++;
								}
						}
					}
				}
		
			t=0;//The number of merging pairs
			k=C;//The number of cluster before mering
            for(i=0;i<dist_i.size();i++)
			{	
				if(IsCombine[dist_i[i]]==0&&IsCombine[dist_j[i]]==0)
				{
					
					for(j=0;j<m_n;j++)
					{
						znew=(m_Ne[dist_i[i]]*m_z[j][dist_j[i]]+m_Ne[dist_j[i]]*m_z[j][dist_j[i]])/(m_Ne[dist_i[i]]+m_Ne[dist_j[i]]);
						m_z[j].push_back(znew);
					}
					
					C++;
				    IsCombine[dist_i[i]]=IsCombine[dist_j[i]]=1;
					t++;
					if(t>m_L) break;
				}
			}
		   
			for(j=k-1;j>=0;j--)
				if(IsCombine[j]==1)
				{
					for(i=0;i<m_n;i++)
						m_z[i].erase(m_z[i].begin()+j);
					C--;
				}
		}
		
	}
	return m_ClassLabel;
}
//////////////////////////////////////////////////////////////////////////////
/* Name£ºNumOfClass
   Function£ºReturn the number of clusters after classification
   Return£ºThe number of clusters
*/
long* CISODATA::NumOfEachClass()
{
	if(m_z.size()==0) return 0;
	long *num;
	num=new long[m_z[0].size()];
	for(int i=0;i<m_z[0].size();i++)
		num[i]=m_Ne[i];
    return num;
}
///////////////////////////////////////////////////////////////////////////////
/* Name£ºClusterCenter
   Function: Return the cluster center of all clusters
   Return£ºThe vector of cluster centers
*/

vector<vector<double> >CISODATA::ClusterCenter()
{
	return m_z;
}
///////////////////////////////////////////////////////////////////////////////
/*Name£ºRMSoEachClass
  Function£ºReturn the stadard deviation vector of each cluster
  Return£ºThe vector of standard deviation
*/
vector<vector<double> >CISODATA::RMSofEachClass()
{
	return m_s;
}
//////////////////////////////////////////////////////////////////////////////
/*Name£ºAllClassDistance
  Function: Return distances between all clusters
  Return£ºThe two dimensional vector of all distances
*/
vector<vector<double> >CISODATA::AllClassDistances()
{
	vector<vector<double> > distance;
	int C=m_z[0].size(),i,j;
	distance.resize(C);
	
	for(i=0;i<C;i++)
		distance[i].resize(C);
    for(i=0;i<C;i++)
	{
		for(j=i+1;j<C;j++)
		{
			distance[i][j]=ClassDistance(i,j);
			distance[j][i]=distance[i][j];
		}
		distance[i][i]=0;
	}
	return distance;
}
///////////////////////////////////////////////////////////////////////////	
/*Name:NumofClasses
  Function£ºReturnt the number of clusters
  Return£ºThe number of clusters
*/
int CISODATA::NumOfClasses()
{
	return m_z[0].size();
}

///////////////////////////////////////////////////////////////////////////
/*Name£ºMeanDistance
  Function£ºCalculate the average distance between all members(pixels) and the cluster center for each cluster after classification
  Return£ºThe pointer to distance vector
*/
double *CISODATA::MeanDistance()
{
	double *md=new double[NumOfClasses()];
	for(int i=0;i<NumOfClasses();i++)
		md[i]=m_de[i];
	return md;
}
/*Name£ºTotalMeanDistance
  Function£ºThe total average distance between all samples and the cluster center in one cluster
  Return£ºDistance
*/
double CISODATA::TotalMeanDistance()
{
	return m_dm;
}