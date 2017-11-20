//数据仓库与数据挖掘 杨东升 1400012898
//pagerank

//数据集采用作者间论文引用的网络
//点：作者
//有向边：单向的引用关系 

#include <iostream>
#include <fstream>  
#include <string>
#include <cstring>
#include <map> 
#include <algorithm>
#include <vector>
#include <windows.h>
using namespace std;

//最大迭代次数
#define maxiter 1000
//收敛阈值 
#define tiny 1e-10
//收敛值的指数
#define tinyIndex 10

//最大点数 
#define maxauthor 17000
//最大边数
#define maxedge 600000 
//单个点的最大入边数 
#define maxindegree 7500
//随机游走的概率  
#define C 0.15
//名字最大长度
#define maxlength 100 

//点和对应的权值构成pair 
typedef pair<int, double> PAIR;

//编号到名字的映射 
map<int, string> isIndex;
//名字到编号的映射 
map<string, int> siIndex;
//用于对权值排序的数据结构 
vector<PAIR> prList ;
//记录被引用关系的稀疏矩阵 
unsigned short citedBy[maxauthor][maxindegree] ; 
//记录边和边权的稠密矩阵 
float matrix[maxauthor][maxauthor] ;
//记录所有无出边的点 
int zero[maxauthor] ;

//点数，边数，无出边点数 
int ntotal, etotal, cntz;
//记录每个点的出边数 
int cnt_out[maxauthor] ;
//记录每个点的入边数 
int cnt_in[maxauthor] ;
//记录点的出边关联点的入边总数
int cnt_out_in[maxauthor] ;
//记录点的值是否已经收敛
bool fix[maxauthor] ; 
//记录边的权值
float weight[maxauthor][maxindegree] ;
//记录点的旧权值和新权值 
double val[maxauthor], newval[maxauthor] ;
//记录点的Hub值 
double hub[maxauthor] ;

//两次迭代之间的差值，分别是旧差值与新差值 
double olddif[maxauthor], dif[maxauthor] ; 
//每一次迭代中最大的权值变化量
double change[maxiter] ;
//统计每代最大差值的点的编号
int number[maxiter];
//统计每代最大差值的点的新旧值
double v[maxiter], vn[maxiter] ;
//统计迭代数
int iterate[tinyIndex+1] ;
//统计时间
double time[tinyIndex+1];
//CPU频率
LARGE_INTEGER m_nFreq; 
//时间戳：开始时间，结束时间 
LARGE_INTEGER m_nBeginTime, nEndTime;

//比较函数 
struct CmpByValue {  
  	bool operator()(const PAIR& lhs, const PAIR& rhs) {  
    	return lhs.second > rhs.second;  
  	}  
};
 
void BuildIndex() 
{ 
  	 cout << "start to build index" << endl; 
  	 //打开存储所有作者的名字的文件 
 	 FILE *f1 = fopen("input/author_ids.txt", "r"); 
 	 //打开存储作者间的引用关系的文件 
 	 FILE *f2 = fopen("input/author-citation-nonself-network.txt", "r");
 	 if(f1 == NULL || f2 == NULL) 
 	 { 
  	   	cout << "open input file error" << endl ; 
  	   	return ; 
     }
     
     //读取作者名字，进行数字化，即建立编号和名字之间的双向映射
	 //这样实际计算中可以使用编号来算pagerank 
	 int i ;
  	 for( i = 1; i <= maxauthor; i ++ ) 
  	 { 
	   	int index ;
	   	char name[maxlength] ;
	   	if( fscanf( f1, "%d ", &index ) == 0 ) break ;
	   	if( fgets( name, 100, f1 ) == 0 ) break ;
	   	name[strlen(name)-1] = 0 ; //去掉换行符 
	   	isIndex[index] = name ;
	   	siIndex[name] = index ;
     } 
     //得到作者总数 
     ntotal = i - 1 ;
     
	 //读取作者引用关系，进行数字化
  	 for( i = 1; i <= maxedge; i ++ ) 
  	 { 
   	   	char temp[2*maxlength+5], name1[maxlength], name2[maxlength] ;
   	   	int j = 0 ;
   	   	if( fgets( temp, 2*maxlength+5, f2 ) == 0 ) break ;
   	   	
   	   	//解析temp字符串，以"="为标志，前后分别是两个作者的名字 
   	   	for( j = 0 ; temp[j+1] != '=' ; j ++ )
   	   	{
	   		name1[j] = temp[j] ;
	   	}
	   	name1[j] = 0 ;
	   	int j0 = j+5 ; //此位置是第二个作者名字的起始位置 
	   	for( j = j0 ; temp[j-1] != 0 ; j ++ )
   	   	{
	   		name2[j-j0] = temp[j] ;
	   	}
	   	name2[j-j0-2] = 0 ;
	   	
	   	//如果作者不在数据库中，则不记录此条引用关系 
	   	if( siIndex.find(name1) == siIndex.end() ) continue ;
	   	if( siIndex.find(name2) == siIndex.end() ) continue ;
	   	int index1 = siIndex[name1] ;
	   	int index2 = siIndex[name2] ;
	   
	   	//记录此条引用关系
	   	//index2被index1引用，index1指向index2 
	   	//index2的入边数cnt_in[index2]加1 
	  	//index1的出边数cnt_out[index1]加1 
       	citedBy[index2][cnt_in[index2]] = index1 ;
       	cnt_in[index2] ++ ;
       	cnt_out[index1] ++ ;
       	
       	//表示index1引用了index2 
       	matrix[index2][index1] += 1 ;
     }
     //得到边总数
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
	//初始状态，每个点的权值相同
	//所有点权值总和是1 
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
	
	//根据出边数和随机游走概率，从引用关系计算出转移矩阵 
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
     
     //统计没有出边的点的数量
	 //记录没有出边的点 
     for( int i = 1; i <= ntotal; i ++ ) 
  	 {
        if( cnt_out[i] == 0 ) 
        {
            zero[cntz] = i ;
            cntz ++ ;
            
			//令没有出边的点向所有点连边，更新转移矩阵  
            for( int j = 1; j <= ntotal; j ++ )
			  	matrix[j][i] += (double)(1-C)/ntotal ;
        }
     }
}

void InitValue( int Init )
{
	//所有点权值总和是1
	
	//每个点的初始值相同
	if( Init == 0 )
	{
		for( int i = 1; i <= ntotal; i ++ )
  		{
  			val[i] = 1/double(ntotal) ;
		}
	}
	//初始值与入边数成正比 
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
	//收敛的精度 
	static float accuracy = 0.001 ;
	//精度的指数 
	static int power = 3 ;
	//PR值的最大变化量 
	change[iter] = 0 ;
		
		for( int i = 1; i <= ntotal; i ++ )
        {
			//统计最大变化量
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
        	//用新权值取代旧权值 
            val[i] = newval[i] ;
        }  		
   		//在达到一定的精度时，计时
   		if( change[iter] < accuracy )
   		{
   			accuracy /= 10 ;
			iterate[power] = iter ;
   			QueryPerformanceCounter(&nEndTime);  
            time[power++] = (double)(nEndTime.QuadPart-m_nBeginTime.QuadPart) / m_nFreq.QuadPart ;
		}
        //判断是否满足收敛条件
        if( change[iter] < tiny ) return 1 ;
        else return 0 ;
}

//利用二位外推法加速收敛
void ExtraPolation()
{
	double lamda ;
	double sum1 = 0, sum2 = 0 ; 
	//根据历史值的差值推演新值 
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
   	  	//核心计算是矩阵和向量的乘法 
		for( int i = 1; i <= ntotal; i ++ )
		{
			newval[i] = 0 ;
			for( int j = 1; j <= ntotal; j ++ )
        		newval[i] += matrix[i][j] * val[j]; 
		}
		//检查是否收敛 
		if( CheckChange(iter) == 1 ) break ;
    }
}

void PageRank_sparse(int Extra, int Fix)
{
    for( int iter = 1; iter < maxiter; iter ++ ) 
  	{ 
   	  	cout << "iterate: " << iter << endl ;
   	  	//所有的无出边的点PR值总和
		double sum_val0 = 0 ;
		for( int j = 0 ; j < cntz ; j ++ )
    		sum_val0 += val[zero[j]] ;
    	
    	//计算每个点的新PR值 
		for( int i = 1; i <= ntotal; i ++ ) 
   		{ 
    	  	if( Fix && fix[i] ) continue ;
			//第一部分的PR值是由随机游走得来的
			//两种随机游走：第一种是由C贡献的，第二种是由无出边的点贡献的
			//所有点向所有点连边，无出边的点向所有点连边 
			
			//方案一：规定每个无出边点没有到自己的连边 
			if( cnt_out[i] != 0 ) newval[i] = C/double(ntotal) + (1-C)*sum_val0/double(ntotal-1)  ;
    	  	else newval[i] = C/double(ntotal) + (1-C)*(sum_val0-val[i])/double(ntotal-1)  ;
    	  	
    	  	//方案二：规定每个无出边点有到自己的连边
    	  	//newval[i] = C/ntotal + (1-C)*sum_val0/ntotal  ;
		   	
			//第二部分的PR值由每个点从它的所有入边“吸取” 			 
			for( int j = 0; j < cnt_in[i]; j ++ ) 
    	  	{
				//第i点的新PR值 += 不随机游走概率 * 第i点的第j条入边的权值 * 入点的PR值 
				newval[i] += (1-C) * val[citedBy[i][j]] * weight[i][j] ;
            }        
        }
        //采用二次外推法加速收敛 
        if( Extra == 1 )
        {
        	if( iter == 1 )	//第一次不能使用 
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
        //检查是否收敛 
        if( CheckChange(iter) == 1 ) break ;
    }
}

void HITS()
{
	for( int iter = 1; iter < maxiter; iter ++ ) 
  	{ 
   	  	cout << "iterate: " << iter << endl ;
   	  	//所有的无出边的点hub值总和
		double sum_val0 = 0 ;
		for( int j = 0 ; j < cntz ; j ++ )
    		sum_val0 += val[zero[j]] ;
    	
    	//计算每个点的新值 
		for( int i = 1; i <= ntotal; i ++ ) 
   		{
			//第一部分的authority值是由随机游走得来的
			//两种随机游走：第一种是由C贡献的，第二种是由无出边的点贡献的
			//所有点向所有点连边，无出边的点向所有点连边
    	  	newval[i] = C/ntotal + (1-C)*sum_val0/ntotal  ;
		   	
			//第二部分的authority值由每个点从它的所有入边“吸取”   			 
			for( int j = 0; j < cnt_in[i]; j ++ )
    	  	{
				//第i点的新authority值 += 不随机游走概率 * 第i点的第j条入边的hub值 / 贡献点的出边数 
				newval[i] += (1-C) * val[citedBy[i][j]] / cnt_out[citedBy[i][j]] ;
            }        
        }
        //检查是否收敛 
        if( CheckChange(iter) == 1 ) break ;
    }
} 

void MakeResult(string name)
{
	//统计所有点权值总和，如果计算正确，所有点权值和应该是1 
    double allsum = 0 ;
    for( int i = 1; i <= ntotal; i ++ ) 
   	{ 
   		//插入点和值，用于排序 
		prList.push_back(make_pair(i, val[i])) ;
		allsum += val[i] ;
   	}
   	cout << "the sum of the pageranks: " << allsum << endl ;
    
    //按照权值对点进行排序 
    sort( prList.begin(), prList.end(), CmpByValue() ) ;	 
	
	//打开结果文件 
  	cout << "writing result" << endl; 
  	FILE * f3, * f4 ;
  	
    f3 = fopen((string("output/PageRank_") + name).c_str(), "w"); 
   	if( f3 == NULL ) 
	   cout << "open output file error" << endl ;
   	
   	//把编号转换回作者名字，输出名字和pagerank值
  	for( int i = 0 ; i < ntotal ; i ++ )
  	{
		fprintf( f3, "%s\t%f\n", isIndex[prList[i].first].c_str(), prList[i].second ) ;
	}
	fclose(f3) ;
	
	//输出每次迭代时收敛的情况 	
	f4 = fopen((string("output/statistics_") + name).c_str(), "w"); 
   	if( f4 == NULL ) 
	    cout << "open output file error" << endl ;
	    
	fprintf( f4, "精度指数\t时间\n" ) ;
	for( int power = 3 ; power <= tinyIndex ; power ++ )
	{
		fprintf( f4, "%d\t%llf\n", power, time[power] ) ;
	}
	fprintf( f4, "精度指数\t迭代次数\n" ) ;
	for( int power = 3 ; power <= tinyIndex ; power ++ )
	{
		fprintf( f4, "%d\t%d\n", power, iterate[power] ) ;
	}
	fprintf( f4, "迭代次数\t值变化最大的点\t旧值\t新值\t变化量\n" ) ;
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
	
	QueryPerformanceFrequency(&m_nFreq); // 获取时钟周期 
	
	//数据读取和预处理 
	BuildIndex() ; 
	//统计数据 	
	MakeStatistics(W) ;
	
	QueryPerformanceCounter(&m_nBeginTime); // 获取时钟计数
	//初始化PR值 
	InitValue(I) ;
	if( N )
	{ 	
  		//计算pagerank_naive 
    	PageRank_naive() ;
	}
	if( S )
	{
		PageRank_sparse(E, F) ;
	}
	QueryPerformanceCounter(&nEndTime);  // 获取时钟计数
	//输出总用时 
    cout << (double)(nEndTime.QuadPart-m_nBeginTime.QuadPart) / m_nFreq.QuadPart << endl;
    
	//数据验证和排序输出
	MakeResult(name + ".txt") ;
	
 	return 0;
}
