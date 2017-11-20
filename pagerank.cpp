//���ݲֿ��������ھ� ��� 1400012898
//pagerank

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
//��¼�ߺͱ�Ȩ�ĳ��ܾ��� 
float matrix[maxauthor][maxauthor] ;
//��¼�����޳��ߵĵ� 
int zero[maxauthor] ;

//�������������޳��ߵ��� 
int ntotal, etotal, cntz;
//��¼ÿ����ĳ����� 
int cnt_out[maxauthor] ;
//��¼ÿ���������� 
int cnt_in[maxauthor] ;
//��¼��ĳ��߹�������������
int cnt_out_in[maxauthor] ;
//��¼���ֵ�Ƿ��Ѿ�����
bool fix[maxauthor] ; 
//��¼�ߵ�Ȩֵ
float weight[maxauthor][maxindegree] ;
//��¼��ľ�Ȩֵ����Ȩֵ 
double val[maxauthor], newval[maxauthor] ;
//��¼���Hubֵ 
double hub[maxauthor] ;

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
	 //����ʵ�ʼ����п���ʹ�ñ������pagerank 
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
       	cnt_in[index2] ++ ;
       	cnt_out[index1] ++ ;
       	
       	//��ʾindex1������index2 
       	matrix[index2][index1] += 1 ;
     }
     //�õ�������
     etotal = i - 1 ;
     fclose(f1) ;
     fclose(f2) ;
 	 cout<<"finish building the index"<<endl; 
}  

void MakeStatistics(int weighted)
{
	cout << "number of nodes: " << ntotal <<endl ;
  	cout << "number of edges: " << etotal << endl ;
  	
	int maxnum = 0 ;
	//��ʼ״̬��ÿ�����Ȩֵ��ͬ
	//���е�Ȩֵ�ܺ���1 
  	for( int i = 1 ; i <= ntotal ; i ++ )
  	{
	  	if( cnt_in[i] > maxnum ) maxnum = cnt_in[i] ;
	}
	cout << "maxedge of a node: " << maxnum << endl ;
	
	for( int i = 1; i <= ntotal; i ++ )
		for( int j = 0; j < cnt_in[i]; j ++ )
			cnt_out_in[citedBy[i][j]] += cnt_in[i] ;
	
	if( weighted == 0 )
		for( int i = 1; i <= ntotal; i ++ ) 
			for( int j = 0; j < cnt_in[i]; j ++ ) 
				weight[i][j] = (double)1 / cnt_out[citedBy[i][j]] ;
	else
		for( int i = 1; i <= ntotal; i ++ ) 
			for( int j = 0; j < cnt_in[i]; j ++ )
				weight[i][j] = (double)cnt_in[i] / cnt_out_in[citedBy[i][j]] ;
	
	//���ݳ�������������߸��ʣ������ù�ϵ�����ת�ƾ��� 
     for( int i = 1; i <= ntotal; i ++ )
     {
     	for( int j = 1; j <= ntotal; j ++ )
     	{
		 	if(matrix[i][j] != 0) 
			{
				if( weighted == 0 ) matrix[i][j] = matrix[i][j] * (1-C) / cnt_out[j] ;
				else matrix[i][j] = matrix[i][j] * (1-C) * cnt_in[i] / cnt_out_in[j] ;
			}
		 	matrix[i][j] += (double)C / ntotal ;
		}
	 }
     
     //ͳ��û�г��ߵĵ������
	 //��¼û�г��ߵĵ� 
     for( int i = 1; i <= ntotal; i ++ ) 
  	 {
        if( cnt_out[i] == 0 ) 
        {
            zero[cntz] = i ;
            cntz ++ ;
            
			//��û�г��ߵĵ������е����ߣ�����ת�ƾ���  
            for( int j = 1; j <= ntotal; j ++ )
			  	matrix[j][i] += (double)(1-C)/ntotal ;
        }
     }
}

void InitValue( int Init )
{
	//���е�Ȩֵ�ܺ���1
	
	//ÿ����ĳ�ʼֵ��ͬ
	if( Init == 0 )
	{
		for( int i = 1; i <= ntotal; i ++ )
  		{
  			val[i] = 1/double(ntotal) ;
		}
	}
	//��ʼֵ������������� 
	if( Init == 1 )
	{
		int sum = 0 ;
		for( int i = 1; i <= ntotal; i ++ )
		{
			sum += cnt_in[i] + 1 ;
		}
		for( int i = 1; i <= ntotal; i ++ )
  		{
  			val[i] = (cnt_in[i]+1)/double(sum) ;
		}
	}
}

int CheckChange(int iter)
{
	//�����ľ��� 
	static float accuracy = 0.001 ;
	//���ȵ�ָ�� 
	static int power = 3 ;
	//PRֵ�����仯�� 
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
		    if( newval[i] - val[i] < tiny / 10
				&& -newval[i] + val[i] < tiny / 10 ) fix[i] = 1 ;
        	//����Ȩֵȡ����Ȩֵ 
            val[i] = newval[i] ;
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

//���ö�λ���Ʒ���������
void ExtraPolation()
{
	double lamda ;
	double sum1 = 0, sum2 = 0 ; 
	//������ʷֵ�Ĳ�ֵ������ֵ 
	for( int i = 1; i <= ntotal; i ++ )
	{
		sum1 += olddif[i] * olddif[i] ;
		sum2 += olddif[i] * dif[i] ;
	}
	lamda = sum2/sum1 ;
	for( int i = 1; i <= ntotal; i ++ )
	{
		newval[i] = (lamda * val[i] - newval[i]) / (lamda-1) ;
	}
}

void PageRank_naive()
{
	cout << "pagerank naive" << endl ;
	
	for( int iter = 1; iter <= maxiter/4; iter ++ ) 
  	{ 
   	  	cout << "iterate: " << iter << endl ;
   	  	//���ļ����Ǿ���������ĳ˷� 
		for( int i = 1; i <= ntotal; i ++ )
		{
			newval[i] = 0 ;
			for( int j = 1; j <= ntotal; j ++ )
        		newval[i] += matrix[i][j] * val[j]; 
		}
		//����Ƿ����� 
		if( CheckChange(iter) == 1 ) break ;
    }
}

void PageRank_sparse(int Extra, int Fix)
{
    for( int iter = 1; iter < maxiter; iter ++ ) 
  	{ 
   	  	cout << "iterate: " << iter << endl ;
   	  	//���е��޳��ߵĵ�PRֵ�ܺ�
		double sum_val0 = 0 ;
		for( int j = 0 ; j < cntz ; j ++ )
    		sum_val0 += val[zero[j]] ;
    	
    	//����ÿ�������PRֵ 
		for( int i = 1; i <= ntotal; i ++ ) 
   		{ 
    	  	if( Fix && fix[i] ) continue ;
			//��һ���ֵ�PRֵ����������ߵ�����
			//����������ߣ���һ������C���׵ģ��ڶ��������޳��ߵĵ㹱�׵�
			//���е������е����ߣ��޳��ߵĵ������е����� 
			
			//����һ���涨ÿ���޳��ߵ�û�е��Լ������� 
			if( cnt_out[i] != 0 ) newval[i] = C/double(ntotal) + (1-C)*sum_val0/double(ntotal-1)  ;
    	  	else newval[i] = C/double(ntotal) + (1-C)*(sum_val0-val[i])/double(ntotal-1)  ;
    	  	
    	  	//���������涨ÿ���޳��ߵ��е��Լ�������
    	  	//newval[i] = C/ntotal + (1-C)*sum_val0/ntotal  ;
		   	
			//�ڶ����ֵ�PRֵ��ÿ���������������ߡ���ȡ�� 			 
			for( int j = 0; j < cnt_in[i]; j ++ ) 
    	  	{
				//��i�����PRֵ += ��������߸��� * ��i��ĵ�j����ߵ�Ȩֵ * ����PRֵ 
				newval[i] += (1-C) * val[citedBy[i][j]] * weight[i][j] ;
            }        
        }
        //���ö������Ʒ��������� 
        if( Extra == 1 )
        {
        	if( iter == 1 )	//��һ�β���ʹ�� 
        	{
        		for( int i = 1; i <= ntotal; i++ )
					dif[i] = newval[i] - val[i] ;
			}
			else
			{
				for( int i = 1; i <= ntotal; i++ )
				{
					olddif[i] = dif[i] ;
					dif[i] = newval[i] - val[i] ;
				}
				ExtraPolation() ;
			}
		}
        //����Ƿ����� 
        if( CheckChange(iter) == 1 ) break ;
    }
}

void HITS()
{
	for( int iter = 1; iter < maxiter; iter ++ ) 
  	{ 
   	  	cout << "iterate: " << iter << endl ;
   	  	//���е��޳��ߵĵ�hubֵ�ܺ�
		double sum_val0 = 0 ;
		for( int j = 0 ; j < cntz ; j ++ )
    		sum_val0 += val[zero[j]] ;
    	
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
				newval[i] += (1-C) * val[citedBy[i][j]] / cnt_out[citedBy[i][j]] ;
            }        
        }
        //����Ƿ����� 
        if( CheckChange(iter) == 1 ) break ;
    }
} 

void MakeResult(string name)
{
	//ͳ�����е�Ȩֵ�ܺͣ����������ȷ�����е�Ȩֵ��Ӧ����1 
    double allsum = 0 ;
    for( int i = 1; i <= ntotal; i ++ ) 
   	{ 
   		//������ֵ���������� 
		prList.push_back(make_pair(i, val[i])) ;
		allsum += val[i] ;
   	}
   	cout << "the sum of the pageranks: " << allsum << endl ;
    
    //����Ȩֵ�Ե�������� 
    sort( prList.begin(), prList.end(), CmpByValue() ) ;	 
	
	//�򿪽���ļ� 
  	cout << "writing result" << endl; 
  	FILE * f3, * f4 ;
  	
    f3 = fopen((string("output/PageRank_") + name).c_str(), "w"); 
   	if( f3 == NULL ) 
	   cout << "open output file error" << endl ;
   	
   	//�ѱ��ת�����������֣�������ֺ�pagerankֵ
  	for( int i = 0 ; i < ntotal ; i ++ )
  	{
		fprintf( f3, "%s\t%f\n", isIndex[prList[i].first].c_str(), prList[i].second ) ;
	}
	fclose(f3) ;
	
	//���ÿ�ε���ʱ��������� 	
	f4 = fopen((string("output/statistics_") + name).c_str(), "w"); 
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
    bool N = 0 , S = 0 , H = 0 , I = 0 , E = 0 , F = 0 , W = 0 ;
    string name ;
    if( argc < 2 )
    {
    	printf("-n for naive pagerank\n") ;
    	printf("-w for weighted edge\n") ;
    	printf("-s for sparse pagerank\n") ;
    	printf("-i for initial value\n") ;
    	printf("-e for extra polation\n") ;
    	printf("-f for fixable value\n") ;
    	return 0 ;
	}
	for (argc--, argv++; argc > 0; argc -= 1, argv += 1) {
        if (!strcmp(*argv, "-n"))
        {
            printf ("naive pagerank\n");
            N = 1 ;
            name += "naive" ;
        }
        else if (!strcmp(*argv, "-s"))
        {
            printf ("sparse matrix multiply pagerank\n");
            S = 1 ;
            name += "sparse" ;
        }
        else if (!strcmp(*argv, "-h"))
        {
            printf ("HITS algorithm\n");
            H = 1 ;
            name += "HITS" ;
        }
        if (!strcmp(*argv, "-i"))
        {
            printf ("with carefully designed initial value\n");
            I = 1 ;
            name += "_init" ;
        }
        if (!strcmp(*argv, "-e"))
        {
            printf ("with extra polation\n");
            E = 1 ;
            name += "_extra" ;
        }
        if (!strcmp(*argv, "-f"))
        {
            printf ("with fixable value\n");
            F = 1 ;
            name += "_fix" ;
        }
        if (!strcmp(*argv, "-w"))
        {
            printf ("with weighted edge\n");
            W = 1 ;
            name += "_weighted" ;
        }
    }
	
	QueryPerformanceFrequency(&m_nFreq); // ��ȡʱ������ 
	
	//���ݶ�ȡ��Ԥ���� 
	BuildIndex() ; 
	//ͳ������ 	
	MakeStatistics(W) ;
	
	QueryPerformanceCounter(&m_nBeginTime); // ��ȡʱ�Ӽ���
	//��ʼ��PRֵ 
	InitValue(I) ;
	if( N )
	{ 	
  		//����pagerank_naive 
    	PageRank_naive() ;
	}
	if( S )
	{
		PageRank_sparse(E, F) ;
	}
	QueryPerformanceCounter(&nEndTime);  // ��ȡʱ�Ӽ���
	//�������ʱ 
    cout << (double)(nEndTime.QuadPart-m_nBeginTime.QuadPart) / m_nFreq.QuadPart << endl;
    
	//������֤���������
	MakeResult(name + ".txt") ;
	
 	return 0;
}
