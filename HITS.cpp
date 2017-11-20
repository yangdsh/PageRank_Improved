//���ݲֿ��������ھ� ��� 1400012898
//HITS

//���ݼ��������߼��������õ�����
//�㣺����
//����ߣ���������ù�ϵ 

#include <iostream>
#include <fstream>  
#include <string>
#include <cstring>
#include <map> 
#include <algorithm>
#include <vector>
#include <windows.h>
using namespace std;

//����������
#define maxiter 1000
//������ֵ 
#define tiny 1e-10
//����ֵ��ָ��
#define tinyIndex 10

//������ 
#define maxauthor 17000
//������
#define maxedge 600000 
//��������������� 
#define maxindegree 7500
//��������������� 
#define maxoutdegree 7500
//������ߵĸ���  
#define C 0.15
//������󳤶�
#define maxlength 100 

//��Ͷ�Ӧ��Ȩֵ����pair 
typedef pair<int, double> PAIR;

//��ŵ����ֵ�ӳ�� 
map<int, string> isIndex;
//���ֵ���ŵ�ӳ�� 
map<string, int> siIndex;
//���ڶ�Ȩֵ��������ݽṹ 
vector<PAIR> prList ;
//��¼�����ù�ϵ��ϡ����� 
unsigned short citedBy[maxauthor][maxindegree] ; 
unsigned short citedBy2[maxauthor][maxoutdegree] ; 
//��¼�����޳��ߵĵ� 
int zero[maxauthor] ;
int zero2[maxauthor] ;

//�������������޳��ߵ��� 
int ntotal, etotal, cntz, cntz2;
//��¼ÿ����ĳ����� 
int cnt_out[maxauthor] ;
//��¼ÿ���������� 
int cnt_in[maxauthor] ;
//��¼��ľ�Ȩֵ����Ȩֵ 
double val[maxauthor], newval[maxauthor] ;
//��¼���Hubֵ 
double hub[maxauthor], newhub[maxauthor] ;

//���ε���֮��Ĳ�ֵ���ֱ��Ǿɲ�ֵ���²�ֵ 
double olddif[maxauthor], dif[maxauthor] ; 
//ÿһ�ε���������Ȩֵ�仯��
double change[maxiter] ;
//ͳ��ÿ������ֵ�ĵ�ı��
int number[maxiter];
//ͳ��ÿ������ֵ�ĵ���¾�ֵ
double v[maxiter], vn[maxiter] ;
//ͳ�Ƶ�����
int iterate[tinyIndex+1] ;
//ͳ��ʱ��
double time[tinyIndex+1];
//CPUƵ��
LARGE_INTEGER m_nFreq; 
//ʱ�������ʼʱ�䣬����ʱ�� 
LARGE_INTEGER m_nBeginTime, nEndTime;

//�ȽϺ��� 
struct CmpByValue {  
  	bool operator()(const PAIR& lhs, const PAIR& rhs) {  
    	return lhs.second > rhs.second;  
  	}  
};
 
void BuildIndex() 
{ 
  	 cout << "start to build index" << endl; 
  	 //�򿪴洢�������ߵ����ֵ��ļ� 
 	 FILE *f1 = fopen("input/author_ids.txt", "r"); 
 	 //�򿪴洢���߼�����ù�ϵ���ļ� 
 	 FILE *f2 = fopen("input/author-citation-nonself-network.txt", "r");
 	 if(f1 == NULL || f2 == NULL) 
 	 { 
  	   	cout << "open input file error" << endl ; 
  	   	return ; 
     }
     
     //��ȡ�������֣��������ֻ�����������ź�����֮���˫��ӳ��
	 //����ʵ�ʼ����п���ʹ�ñ������authority 
	 int i ;
  	 for( i = 1; i <= maxauthor; i ++ ) 
  	 { 
	   	int index ;
	   	char name[maxlength] ;
	   	if( fscanf( f1, "%d ", &index ) == 0 ) break ;
	   	if( fgets( name, 100, f1 ) == 0 ) break ;
	   	name[strlen(name)-1] = 0 ; //ȥ�����з� 
	   	isIndex[index] = name ;
	   	siIndex[name] = index ;
     } 
     //�õ��������� 
     ntotal = i - 1 ;
     
	 //��ȡ�������ù�ϵ���������ֻ�
  	 for( i = 1; i <= maxedge; i ++ ) 
  	 { 
   	   	char temp[2*maxlength+5], name1[maxlength], name2[maxlength] ;
   	   	int j = 0 ;
   	   	if( fgets( temp, 2*maxlength+5, f2 ) == 0 ) break ;
   	   	
   	   	//����temp�ַ�������"="Ϊ��־��ǰ��ֱ����������ߵ����� 
   	   	for( j = 0 ; temp[j+1] != '=' ; j ++ )
   	   	{
	   		name1[j] = temp[j] ;
	   	}
	   	name1[j] = 0 ;
	   	int j0 = j+5 ; //��λ���ǵڶ����������ֵ���ʼλ�� 
	   	for( j = j0 ; temp[j-1] != 0 ; j ++ )
   	   	{
	   		name2[j-j0] = temp[j] ;
	   	}
	   	name2[j-j0-2] = 0 ;
	   	
	   	//������߲������ݿ��У��򲻼�¼�������ù�ϵ 
	   	if( siIndex.find(name1) == siIndex.end() ) continue ;
	   	if( siIndex.find(name2) == siIndex.end() ) continue ;
	   	int index1 = siIndex[name1] ;
	   	int index2 = siIndex[name2] ;
	   
	   	//��¼�������ù�ϵ
	   	//index2��index1���ã�index1ָ��index2 
	   	//index2�������cnt_in[index2]��1 
	  	//index1�ĳ�����cnt_out[index1]��1 
       	citedBy[index2][cnt_in[index2]] = index1 ;
       	citedBy2[index1][cnt_out[index1]] = index2 ;
       	cnt_in[index2] ++ ;
       	cnt_out[index1] ++ ;
     }
     //�õ�������
     etotal = i - 1 ;
     fclose(f1) ;
     fclose(f2) ;
 	 cout<<"finish building the index"<<endl; 
}  

void MakeStatistics()
{
	cout << "number of nodes: " << ntotal <<endl ;
  	cout << "number of edges: " << etotal << endl ;
  	
	int maxnum = 0, maxnum2 = 0 ;
	//��ʼ״̬��ÿ�����Ȩֵ��ͬ
	//���е�Ȩֵ�ܺ���1 
  	for( int i = 1 ; i <= ntotal ; i ++ )
  	{
	  	if( cnt_in[i] > maxnum ) maxnum = cnt_in[i] ;
	  	if( cnt_out[i] > maxnum2 ) maxnum2 = cnt_out[i] ;
	}
	cout << "max in_edge of a node: " << maxnum << endl ;
	cout << "max out_edge of a node: " << maxnum2 << endl ;
     
     //ͳ��û�г��ߵĵ������
	 //��¼û�г��ߵĵ� 
     for( int i = 1; i <= ntotal; i ++ ) 
  	 {
        if( cnt_out[i] == 0 ) 
        {
            zero[cntz] = i ;
            cntz ++ ;
        }
        if( cnt_in[i] == 0 ) 
        {
            zero2[cntz2] = i ;
            cntz2 ++ ;
        }
     }
}

void InitValue()
{
	//���е�Ȩֵ�ܺ���1	
	//ÿ����ĳ�ʼֵ��ͬ
	for( int i = 1; i <= ntotal; i ++ )
	{
  		val[i] = 1/double(ntotal) ;
  		hub[i] = val[i] ;
	}
}

int CheckChange(int iter)
{
	//�����ľ��� 
	static float accuracy = 0.001 ;
	//���ȵ�ָ�� 
	static int power = 3 ;
	//authorityֵ�����仯�� 
	change[iter] = 0 ;
		
		for( int i = 1; i <= ntotal; i ++ )
        {
			//ͳ�����仯��
   			if( val[i] - newval[i] > change[iter] ) 
   			{
			   number[iter] = i ;
			   v[iter] = val[i] ;
			   vn[iter] = newval[i] ;
			   change[iter] = val[i] - newval[i] ;
			}
   			else if( newval[i] - val[i] > change[iter] ) 
			{
			   number[iter] = i ;
			   v[iter] = val[i] ;
			   vn[iter] = newval[i] ;
			   change[iter] = newval[i] - val[i] ;
		    }
        }  		
   		//�ڴﵽһ���ľ���ʱ����ʱ
   		if( change[iter] < accuracy )
   		{
   			accuracy /= 10 ;
			iterate[power] = iter ;
   			QueryPerformanceCounter(&nEndTime);  
            time[power++] = (double)(nEndTime.QuadPart-m_nBeginTime.QuadPart) / m_nFreq.QuadPart ;
		}
        //�ж��Ƿ�������������
        if( change[iter] < tiny ) return 1 ;
        else return 0 ;
}

void HITS()
{
	for( int iter = 1; iter < maxiter; iter ++ ) 
  	{ 
   	  	cout << "iterate: " << iter << endl ;
   	  	//���е��޳��ߵĵ�hubֵ�ܺ�
		double sum_val0 = 0, sum_val02 = 0 ;
		for( int j = 0 ; j < cntz ; j ++ )
    		sum_val0 += hub[zero[j]] ;
    	
    	//����ÿ�������ֵ 
		for( int i = 1; i <= ntotal; i ++ ) 
   		{
			//��һ���ֵ�authorityֵ����������ߵ�����
			//����������ߣ���һ������C���׵ģ��ڶ��������޳��ߵĵ㹱�׵�
			//���е������е����ߣ��޳��ߵĵ������е�����
    	  	newval[i] = C/ntotal + (1-C)*sum_val0/ntotal  ;
		   	
			//�ڶ����ֵ�authorityֵ��ÿ���������������ߡ���ȡ��   			 
			for( int j = 0; j < cnt_in[i]; j ++ )
    	  	{
				//��i�����authorityֵ += ��������߸��� * ��i��ĵ�j����ߵ�hubֵ / ���׵�ĳ����� 
				newval[i] += (1-C) * hub[citedBy[i][j]] / cnt_out[citedBy[i][j]] ;
            }        
        }
        if( CheckChange(iter) == 1 ) break ;
        //����Ȩֵȡ����Ȩֵ 
        for( int i = 1 ; i <= ntotal ; i ++ )
			val[i] = newval[i] ;
			
		for( int j = 0 ; j < cntz2 ; j ++ )
    		sum_val02 += val[zero2[j]] ;
        for( int i = 1; i <= ntotal; i ++ ) 
   		{
			//��һ���ֵ�hubֵ����������ߵ�����
			//����������ߣ���һ������C���׵ģ��ڶ��������޳��ߵĵ㹱�׵�
			//���е������е����ߣ��޳��ߵĵ������е�����
    	  	newhub[i] = C/ntotal + (1-C)*sum_val02/ntotal  ;
		   	
			//�ڶ����ֵ�hubֵ��ÿ���������������ߡ���ȡ��   			 
			for( int j = 0; j < cnt_out[i]; j ++ )
    	  	{
				//��i�����authorityֵ += ��������߸��� * ��i��ĵ�j����ߵ�authorityֵ / ���׵�ĳ����� 
				newhub[i] += (1-C) * val[citedBy2[i][j]] / cnt_in[citedBy2[i][j]] ;
            }        
        }
        for( int i = 1 ; i <= ntotal ; i ++ )
        	hub[i] = newhub[i] ;
    }
} 

void MakeResult()
{
	//ͳ�����е�Ȩֵ�ܺͣ����������ȷ�����е�Ȩֵ��Ӧ����1 
    double allsum = 0 ;
    for( int i = 1; i <= ntotal; i ++ ) 
   	{ 
   		//������ֵ���������� 
		prList.push_back(make_pair(i, val[i])) ;
		allsum += val[i] ;
   	}
   	cout << "the sum of authority: " << allsum << endl ;
    
    //����Ȩֵ�Ե�������� 
    sort( prList.begin(), prList.end(), CmpByValue() ) ;	 
	
	//�򿪽���ļ� 
  	cout << "writing result" << endl; 
  	FILE * f3, * f4 ;
  	
    f3 = fopen("output/HITS.txt", "w"); 
   	if( f3 == NULL ) 
	   cout << "open output file error" << endl ;
   	
   	//�ѱ��ת�����������֣�������ֺ�authorityֵ
  	for( int i = 0 ; i < ntotal ; i ++ )
  	{
		fprintf( f3, "%s\t%f\n", isIndex[prList[i].first].c_str(), prList[i].second ) ;
	}
	fclose(f3) ;
	
	//���ÿ�ε���ʱ��������� 	
	f4 = fopen("output/statistics_HITS.txt", "w"); 
   	if( f4 == NULL ) 
	    cout << "open output file error" << endl ;
	    
	fprintf( f4, "����ָ��\tʱ��\n" ) ;
	for( int power = 3 ; power <= tinyIndex ; power ++ )
	{
		fprintf( f4, "%d\t%llf\n", power, time[power] ) ;
	}
	fprintf( f4, "����ָ��\t��������\n" ) ;
	for( int power = 3 ; power <= tinyIndex ; power ++ )
	{
		fprintf( f4, "%d\t%d\n", power, iterate[power] ) ;
	}
	fprintf( f4, "��������\tֵ�仯���ĵ�\t��ֵ\t��ֵ\t�仯��\n" ) ;
	for( int i = 1 ; i <= maxiter ; i ++ )
  	{
		fprintf( f4, "%d\t%d\t%.30llf\t%.30llf\t%.30llf\n", i, number[i], v[i], vn[i], change[i] ) ;
		if(change[i] < tiny) break ;
	}
	fclose(f4) ;
	prList.clear() ;
}

int main(int argc, char *argv[]) 
{  
	QueryPerformanceFrequency(&m_nFreq); // ��ȡʱ������ 
	
	//���ݶ�ȡ��Ԥ���� 
	BuildIndex() ; 
	//ͳ������ 	
	MakeStatistics() ;
	
	QueryPerformanceCounter(&m_nBeginTime); // ��ȡʱ�Ӽ���
	//��ʼ��authorityֵ 
	InitValue() ;
	HITS() ; 
	QueryPerformanceCounter(&nEndTime);  // ��ȡʱ�Ӽ���
	//�������ʱ 
    cout << (double)(nEndTime.QuadPart-m_nBeginTime.QuadPart) / m_nFreq.QuadPart << endl;
    
	//������֤���������
	MakeResult() ;
	
 	return 0;
}
