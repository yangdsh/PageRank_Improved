//数据仓库与数据挖掘 杨东升 1400012898
//HITS

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
//单个点的最大出边数 
#define maxoutdegree 7500
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
unsigned short citedBy2[maxauthor][maxoutdegree] ; 
//记录所有无出边的点 
int zero[maxauthor] ;
int zero2[maxauthor] ;

//点数，边数，无出边点数 
int ntotal, etotal, cntz, cntz2;
//记录每个点的出边数 
int cnt_out[maxauthor] ;
//记录每个点的入边数 
int cnt_in[maxauthor] ;
//记录点的旧权值和新权值 
double val[maxauthor], newval[maxauthor] ;
//记录点的Hub值 
double hub[maxauthor], newhub[maxauthor] ;

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
	 //这样实际计算中可以使用编号来算authority 
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
       	citedBy2[index1][cnt_out[index1]] = index2 ;
       	cnt_in[index2] ++ ;
       	cnt_out[index1] ++ ;
     }
     //得到边总数
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
	//初始状态，每个点的权值相同
	//所有点权值总和是1 
  	for( int i = 1 ; i <= ntotal ; i ++ )
  	{
	  	if( cnt_in[i] > maxnum ) maxnum = cnt_in[i] ;
	  	if( cnt_out[i] > maxnum2 ) maxnum2 = cnt_out[i] ;
	}
	cout << "max in_edge of a node: " << maxnum << endl ;
	cout << "max out_edge of a node: " << maxnum2 << endl ;
     
     //统计没有出边的点的数量
	 //记录没有出边的点 
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
	//所有点权值总和是1	
	//每个点的初始值相同
	for( int i = 1; i <= ntotal; i ++ )
	{
  		val[i] = 1/double(ntotal) ;
  		hub[i] = val[i] ;
	}
}

int CheckChange(int iter)
{
	//收敛的精度 
	static float accuracy = 0.001 ;
	//精度的指数 
	static int power = 3 ;
	//authority值的最大变化量 
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

void HITS()
{
	for( int iter = 1; iter < maxiter; iter ++ ) 
  	{ 
   	  	cout << "iterate: " << iter << endl ;
   	  	//所有的无出边的点hub值总和
		double sum_val0 = 0, sum_val02 = 0 ;
		for( int j = 0 ; j < cntz ; j ++ )
    		sum_val0 += hub[zero[j]] ;
    	
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
				newval[i] += (1-C) * hub[citedBy[i][j]] / cnt_out[citedBy[i][j]] ;
            }        
        }
        if( CheckChange(iter) == 1 ) break ;
        //用新权值取代旧权值 
        for( int i = 1 ; i <= ntotal ; i ++ )
			val[i] = newval[i] ;
			
		for( int j = 0 ; j < cntz2 ; j ++ )
    		sum_val02 += val[zero2[j]] ;
        for( int i = 1; i <= ntotal; i ++ ) 
   		{
			//第一部分的hub值是由随机游走得来的
			//两种随机游走：第一种是由C贡献的，第二种是由无出边的点贡献的
			//所有点向所有点连边，无出边的点向所有点连边
    	  	newhub[i] = C/ntotal + (1-C)*sum_val02/ntotal  ;
		   	
			//第二部分的hub值由每个点从它的所有入边“吸取”   			 
			for( int j = 0; j < cnt_out[i]; j ++ )
    	  	{
				//第i点的新authority值 += 不随机游走概率 * 第i点的第j条入边的authority值 / 贡献点的出边数 
				newhub[i] += (1-C) * val[citedBy2[i][j]] / cnt_in[citedBy2[i][j]] ;
            }        
        }
        for( int i = 1 ; i <= ntotal ; i ++ )
        	hub[i] = newhub[i] ;
    }
} 

void MakeResult()
{
	//统计所有点权值总和，如果计算正确，所有点权值和应该是1 
    double allsum = 0 ;
    for( int i = 1; i <= ntotal; i ++ ) 
   	{ 
   		//插入点和值，用于排序 
		prList.push_back(make_pair(i, val[i])) ;
		allsum += val[i] ;
   	}
   	cout << "the sum of authority: " << allsum << endl ;
    
    //按照权值对点进行排序 
    sort( prList.begin(), prList.end(), CmpByValue() ) ;	 
	
	//打开结果文件 
  	cout << "writing result" << endl; 
  	FILE * f3, * f4 ;
  	
    f3 = fopen("output/HITS.txt", "w"); 
   	if( f3 == NULL ) 
	   cout << "open output file error" << endl ;
   	
   	//把编号转换回作者名字，输出名字和authority值
  	for( int i = 0 ; i < ntotal ; i ++ )
  	{
		fprintf( f3, "%s\t%f\n", isIndex[prList[i].first].c_str(), prList[i].second ) ;
	}
	fclose(f3) ;
	
	//输出每次迭代时收敛的情况 	
	f4 = fopen("output/statistics_HITS.txt", "w"); 
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
	QueryPerformanceFrequency(&m_nFreq); // 获取时钟周期 
	
	//数据读取和预处理 
	BuildIndex() ; 
	//统计数据 	
	MakeStatistics() ;
	
	QueryPerformanceCounter(&m_nBeginTime); // 获取时钟计数
	//初始化authority值 
	InitValue() ;
	HITS() ; 
	QueryPerformanceCounter(&nEndTime);  // 获取时钟计数
	//输出总用时 
    cout << (double)(nEndTime.QuadPart-m_nBeginTime.QuadPart) / m_nFreq.QuadPart << endl;
    
	//数据验证和排序输出
	MakeResult() ;
	
 	return 0;
}
