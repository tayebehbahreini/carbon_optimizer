#include <iostream>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include<iostream>
#include<list>
#include <random>
#include <omp.h>
using namespace std;
const double hours =4; //simulation  time: 4 hours 
const double rate =1; // scaling factor: increasing the number of submitted workloads by this factor
const int slot_duration =1;//1minutes 
const int T=hours*60/slot_duration;//240 time slots
int TTT=T-15;//to have a feasible solution, we do not schedule jobs arriveed after TTT
const int num_sc=10; // number of scenarios
const int N=rate*15000;//maximum number of jobs in a scenario within T
 
const int M=4;// number of data centers
const int C_range=rate*80;//
const  double avalablae_rate= .0005;
const int sigma_range=20;

double alpha=10;
double ff[M];
int Ns[num_sc];
 int Ind[T][num_sc];//ind[t][sc]:index of the first job arrive after t
struct job
{
    int id;
    int d;
    double l;
    double r;
    int a;
   double p[M];
}jobs[num_sc][N];//jobs[][t][sc]: set of jobs arrive at time t for scenario sc

bool customer_sorter (const job & s1, const job & s2)
    {
        return(s1.d <s2.d);
           
    }

struct center
{
    int id;
    double C;
    double F[T];
    double sigma;
    double av_F[T];
}centers[M];



int bar_x[N][M][T];
int alloc[N][M][num_sc];

int remain[num_sc][N];
int bar_c[num_sc][M][T];
//------------------------------------------------
//------------------------------------------------
struct center_av
{
    int id;
    double av_F;
};
bool customer_sorter1 (const center_av & s1, const center_av & s2)
    {
        return(s1.av_F <s2.av_F);
           
    }


//------------------------------------------------
//------------------------------------------------
double min (double a, double b)
{

    if(a<b)
    return a;
return b;
}

//------------------------------------------------
//------------------------------------------------
int max (int a, int b)
{

    if(a>b)
    return a;
return b;
}

//------------------------------------------------
//------------------------------------------------
void generate_sample() /*generate work load for various scenarios */
{
     ifstream myfile;
    ofstream myfile2;
    double fake=0;
   
    int N_sc[T][num_sc];
    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/workload/ar-Ns.txt");
   
  
    for( int ii=0;ii<TTT;ii++)  
        for(int s=0;s<num_sc;s++)
        {
             myfile>>N_sc[ii][s]; 
            N_sc[ii][s]*= rate;
        }
         
      for( int ii=TTT;ii<T;ii++)  
        for(int s=0;s<num_sc;s++)
         N_sc[ii][s]=0;
        
    
    myfile.close();
  
    double ls[N];
    int rs[N];
    int tmp=0;
    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/workload/l-r.txt");
    /* l-r gives a random sample for the values of l and r*/
    
    for( int ii=0;ii<N;ii++)  
        myfile>>rs[ii]>>tmp>>ls[ii];  
        
    
       
    
    myfile.close();
    //for each scenario in each time slot define arriving workloads
  
   
    int NN[T];
    
    for(int sc=0;sc<num_sc;sc++)   
    {
        int i=0;
        Ns[sc]=0;
  
        for(int t=0;t<T;t++)
        {
         
           for( int ii=0;ii<N_sc[t][sc];ii++)
            {
                int in=rand()%N;
                 jobs[sc][i].l=ceil(ls[in]/(60*slot_duration))+rand()%3;
                jobs[sc][i].l=min( jobs[sc][i].l,10);
              
                 jobs[sc][i].l=min(T-t, jobs[sc][i].l);
                jobs[sc][i].a=t;
                jobs[sc][i].r=min(5,rs[in]);//1+rand()%7;
                jobs[sc][i].d=min(T,jobs[sc][i].a +jobs[sc][i].l+alpha);
               i++; 
            }
           
            Ind[t][sc]=i;
            Ns[sc]+=N_sc[t][sc]; 
           
        }
       
     
       
        }


    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/forecast/F.txt");
    int t=0;
   ofstream myfile1;
    myfile1.open ("/Users/tayebehbahreini/Desktop/cpp/forecast/Fs.txt");
    int cnt=0;
   double sum[M];
   for(int j=0;j<M;j++)
   {
      
    sum[j]=0;
   }
   //using real forecast from quaserClient
   // /*
    while(!myfile.eof() &&t<hours ) //we have info for every   hour (for 8 hours)
    {
        myfile>>ff[0]>>ff[1]>>ff[2]>>ff[3];
        
        for(int tt=0;tt<60/slot_duration;tt++)//repeat each of them for 15 slots 
        {
           
            for(int j=0;j<M;j++)
            {
                centers[j].F[t*60+tt]=ff[j];
                centers[j].av_F[t*60+tt]=centers[j].F[t*60+tt]+sum[j];
                sum[j]=centers[j].av_F[t*60+tt];
            }
           
             for(int j=0;j<M;j++)
                myfile1<<ff[j]<<" ";
             
         
            myfile1<<endl;
           
        }
        t++;

    }
   // */
   //using random C.I
   /*
   for(int t=0;t<60*hours;t++)//repeat each of them for 15 slots 
        {
           
            for(int j=0;j<M;j++)
            {
                if(t%30==0)
                centers[j].F[t]=100*(2+rand()%3);
                else
                centers[j].F[t]=centers[j].F[t-1];

                centers[j].av_F[t]=centers[j].F[t]+sum[j];
                sum[j]=centers[j].av_F[t];
            }
           
             for(int j=0;j<M;j++)
                myfile1<<centers[j].F[t]<<" ";
             
         
            myfile1<<endl;
           
        }
     */
    myfile.close();
     myfile1.close();  
   
    for(int j=0;j<M;j++)
    {
        centers[j].id=j;
        centers[j].C=C_range;//cap[j]*avalablae_rate;//
        
        centers[j].sigma=sigma_range;
        for(int t=0;t<T;t++)
            for(int sc=0;sc<num_sc;sc++)
             bar_c[sc][j][t]=centers[j].C;
       
    }
  
}


//------------------------------------------------
//------------------------------------------------
double MILP_GAP(int tt)
{

int nn=Ind[tt][0];

    double lambda[N][M];
   
    for(int i=0;i<nn;i++)
        for(int j=0;j<M;j++)
        {
            lambda[i][j]=0;
           for(int sc=0;sc<num_sc;sc++)
                lambda[i][j]+=alloc[i][j][sc];
        }
    double zz=-1;;
    IloEnv env;
    
    IloNumVar zzz(env, 0.0, IloInfinity);
    stringstream logfile;
    IloInt i,j,t;


    
    try {
        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray con(env);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
      
        cplex.setOut(logfile);
        
        typedef IloArray<IloNumVarArray> NumVarMatrix;
        typedef IloArray<NumVarMatrix> NumVarMatrix3;
        
        NumVarMatrix x(env,nn);
        for(i=0; i< nn; i++)
        {
            x[i] =  IloNumVarArray(env, M);
            for( j=0; j< M; j++)
                x[i][j]= IloNumVar(env, 0.0, 1.0, ILOINT);
              
                    
        }

       
        
         for(i=0; i< nn; i++)
            for( j=0; j< M; j++)
                model.add(x[i][j]<=lambda[i][j]); 
                    
            

    for(i=0;i<nn;i++)
    {
           IloExpr sumx(env);
            for(j=0;j<M;j++)
                sumx+=x[i][j];
            model.add(sumx<=1);
            
    }
      
         for(j=0;j<M;j++)
        {
                IloExpr sumxi(env);
                for(i=0; i< nn; i++)
                    sumxi+=(x[i][j]*jobs[0][i].r);
                
                 model.add(sumxi<=centers[j].C); 

        }
        IloExpr sumxijt(env);
        for(i=0; i< nn; i++) 
            for(j=0;j<M;j++)
               sumxijt+=(lambda[i][j]* x[i][j]);
               
       
        model.add(zzz<=sumxijt);
      
        model.add(IloMaximize(env, zzz));
        
        cplex.solve();

        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
         
            zz= (double)cplex.getValue(zzz);
           
           
                for(i=0; i< nn; i++)
                {
                    
                     for(j=0;j<M;j++)
                     {
                         bar_x[i][j][tt]=(cplex.getValue(x[i][j]));
                          
                     
                     
                  //  if(tt<10 && bar_x[i][j][tt]>0)
                    // cout<<i<<"\t"<<j<<"\t"<<remain[0][i]<<" "<<jobs[0][i].d-tt<<endl;
                    }
                }
          
        }
  }
    catch (IloException& e)
    {
        cerr << "C-Exp: " << e << endl;
    } catch (...)
    {
        cerr << "Unknown Exception" << endl;
    }
    
    
    
    
    env.end();
    return zz;

}
/*
*/
//------------------------------------------------
//------------------------------------------------
double MILP_det(int tt,int sc)// solve the online problem in time slot tt by CPLEX
{
   
   
    double zz=-1;
    IloEnv env;
    
    IloNumVar zzz(env, 0.0, IloInfinity);
    stringstream logfile;
    IloInt i,j,t;
   int nn=Ns[sc];

   
    //cout<<Ind[tt][sc]<<endl;
  //cout<<nn<<endl;
    try {
        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray con(env);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
        cplex.setOut(logfile);
        
        typedef IloArray<IloNumVarArray> NumVarMatrix;
        typedef IloArray<NumVarMatrix> NumVarMatrix3;
        
        NumVarMatrix3 x(env,nn);
        NumVarMatrix y(env,nn);
        for(i=0; i< nn; i++)
        {
            x[i] =  NumVarMatrix(env, M);
            y[i]=IloNumVarArray(env, M);
           
            for( j=0; j< M; j++)
            {
                  x[i][j]=IloNumVarArray(env,T);
                  y[i][j]=IloNumVar(env, 0.0, 1.0, ILOINT);
                  for( t=0;t<T;t++)
                    x[i][j][t]= IloNumVar(env, 0.0, 1.0, ILOINT);
            }
                    
        }
     for(i=0;i<nn;i++)
    {
             IloExpr sumy(env);
            for(j=0;j<M;j++)
                sumy+=y[i][j];
            model.add(sumy==1);
    }


    for(i=0;i<nn;i++)
        for(j=0;j<M;j++)
          for( t=0;t<tt;t++)
           model.add(x[i][j][t]==bar_x[i][j][t]);
  
    for(i=0;i<nn;i++)
        for(j=0;j<M;j++)
        {
             IloExpr sumxi(env);
     
            for( t=jobs[sc][i].a;t<jobs[sc][i].d;t++)
                sumxi+=(x[i][j][t]);
              
            
           model.add(sumxi==y[i][j]*jobs[sc][i].l); 
           
             for( t=0;t<jobs[sc][i].a;t++)
                model.add(x[i][j][t]==0);


              for( t=jobs[sc][i].d;t<T;t++)
               model.add(x[i][j][t]==0);

        }






    for(j=0;j<M;j++)
            for( t=0;t<T;t++)
            {
                IloExpr sumxi(env);
                for(i=0; i< nn; i++)
                    sumxi+=(x[i][j][t]*jobs[sc][i].r);
                
                 model.add(sumxi<=centers[j].C); 

            }
    IloExpr sumxijt(env);
        for(i=0; i< nn; i++) 
        
            for(j=0;j<M;j++)
                for(int t=0;t<T;t++) 
                    sumxijt+=((centers[j].F[t]*centers[j].sigma*jobs[sc][i].r*x[i][j][t])/centers[j].C);
               
   
        model.add(zzz>=sumxijt);
      
        model.add(IloMinimize(env, zzz));
        
        cplex.solve();
       for(j=0;j<M;j++)
                for(i=0; i< nn; i++)
                     alloc[i][j][sc]=0;
        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
           
            zz= (double)cplex.getValue(zzz);
            cout<<"yes"<<zz<<endl;
       
            
            for(j=0;j<M;j++)
                for(i=0; i< nn; i++)
                     alloc[i][j][sc]=(cplex.getValue(x[i][j][tt]));
                    
                    
                
                  
       }
       else
       cout<<"no "<<tt<<endl;
             
  
    


 }
       
    catch (IloException& e)
    {
        cerr << "C-Exp: " << e << endl;
    } catch (...)
    {
        cerr << "Unknown Exception" << endl;
    }
    
    
    
    
    env.end();
    return zz;
}



//------------------------------------------------
//------------------------------------------------
double MILP_single(int tt,int i,int j,int sc)
/* find the best schedule for the remining work of job i on data center j*/
{
    double zz=-1;
    IloEnv env;
    IloNumVar zzz(env, 0.0, IloInfinity);
    stringstream logfile;
    IloInt t;
   
    try {
        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray con(env);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
        cplex.setOut(logfile);
        
        typedef IloArray<IloNumVarArray> NumVarMatrix;
        
       IloNumVarArray x(env,T);
       
           
            for( t=0;t<T;t++)
                x[t]=IloNumVar(env, 0.0, 1.0, ILOINT);
                     
            IloExpr sumxi(env);
     
            for( t=max(tt,jobs[sc][i].a);t<jobs[sc][i].d;t++)
                sumxi+=(x[t]);

            
            model.add(sumxi==remain[sc][i]); 
           
             for( t=0;t<max(tt,jobs[sc][i].a);t++)
                model.add(x[t]==0);


              for( t=jobs[sc][i].d;t<T;t++)
               model.add(x[t]==0);

        
            for( t=max(tt,jobs[sc][i].a);t<jobs[sc][i].d;t++)
            {
                 IloExpr sumxijt(env);
                 
                model.add(x[t]*jobs[sc][i].r<=bar_c[sc][j][t]); 
     }

    IloExpr sumxijt(env);
       for( t=max(tt,jobs[sc][i].a);t<jobs[sc][i].d;t++)
                    sumxijt+=((centers[j].F[t]*centers[j].sigma*jobs[sc][i].r*x[t])/centers[j].C);
               
   
        model.add(zzz>=sumxijt);
      
        model.add(IloMinimize(env, zzz));
        
        cplex.solve();
       
        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
           
            zz= (double)cplex.getValue(zzz);
           // cout<<"ok"<<zz<<endl;
       
            
            
                //for( t=jobs[i].a;t<jobs[sc][i].d;t++)
                    alloc[i][j][sc]=cplex.getValue(x[tt]);
                    if(alloc[i][j][sc]==1)
                        remain[sc][i]-=1; 
                    for( t=tt;t<T;t++)
                        bar_c[sc][j][t]-=(jobs[sc][i].r*cplex.getValue(x[t]));
           //cout<<"yes";             
                  
       }
       else
       ;//cout<<"no "<<endl;
 }
       
    catch (IloException& e)
    {
        cerr << "C-Exp: " << e << endl;
    } catch (...)
    {
        cerr << "Unknown Exception" << endl;
    }
    
    
    
    
    env.end();
    return zz;
}


//------------------------------------------------
//------------------------------------------------
/*bool schedule(int i,int j,int tt, int sc)
{
    //if(G_y[sc][i]>-1 && G_y[sc][i]!=j)//has been assigned to another
      //  return false;
    

    
    int d=jobs[sc][i].d;
    bool result=true;
    int remains=remain[sc][i];
    int CC[T];
    for(int t=0;t<T;t++)
        CC[t]=bar_c[sc][j][t];
    
    
    slot slots[T] ; 
    int cnt=0;
    for(int t=tt;t<d;t++)
    {
        slots[cnt].F=centers[j].F[t];
        slots[cnt].t=t;
        cnt++;
    }
     sort( slots, slots+d-tt, &customer_sorter2);
   
    for(int ts=0;ts<d-tt;ts++)
    {
        int t=slots[ts].t;
        //cout<<t<<endl;
       
        if(CC[t]>=jobs[sc][i].r && jobs[sc][i].a<=t )
        {
            
            CC[t]-=jobs[sc][i].r;
            remains-=1;
            if(t==tt)
                alloc[i][j][sc]=1;
        }
    }
    if(remains>0)
    { 
        alloc[i][j][sc]=0;
       // cout<<"oops";
        return false;
    }
    else
    {
        //cout<<"yay"<<endl;
        for(int t=tt;t<d;t++)
             bar_c[sc][j][t]=CC[t];
             
             G_y[sc][i]=j;
             if(alloc[i][j][sc]==1)
                remain[sc][i]-=1;
              
             
          
    }

    



    return true;


}
*/

//------------------------------------------------
//------------------------------------------------

double greedy(int tt,int sc,int y[N],int rnd)
{
    
    center_av tmp[M];
    int n=Ns[sc];
    double result=-1;
    for(int j=0;j<M;j++)
        for(int t=tt;t<T;t++)
            bar_c[sc][j][t]=centers[j].C;
    for(int j=0;j<M;j++)
        for(int i=0; i< n; i++)
            alloc[i][j][sc]=0;
            int j=0;
    for(int i=0;i<n;i++)
    {
        int d=jobs[sc][i].d;
        
        
        if(remain[sc][i]>0 &&tt<d )
        {
            int d=jobs[sc][i].d-1;
            int a=max(0,jobs[sc][i].a-1);
          for(int j=0;j<M;j++)
            {
                tmp[j].id=j;
                tmp[j].av_F=centers[j].av_F[d]-centers[j].av_F[a];
            }
            int temp[M];
            
            if(rnd==0)
                sort( tmp, tmp+M, &customer_sorter1);
                
        for( j=0;j<M;j++)
        {
            //if(y[i]==tmp[j].id || y[i]<0)
            result= MILP_single( tt,i,tmp[j].id,sc) ;  
            if(result>0)
            {
                
                y[i]=tmp[j].id;
                break; 
            }
        }
        if(j==M)
            cout<<"inf"<<tt<<" "<<i<<" "<<remain[0][i]<<endl;
       
         
        }

    }
    

     

    return 0;
}


//------------------------------------------------
//------------------------------------------------
double MILP_QoS()//our benchmark with average completion time as objcetive
{
    int sc=0;
    double zz=-1;
    IloEnv env;
    
    IloNumVar zzz(env, 0.0, IloInfinity);
    stringstream logfile;
    IloInt i,j,t;
   int nn=Ns[sc];
  //cout<<nn<<endl;
    try {
        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray con(env);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.95);
        cplex.setOut(logfile);
        
        typedef IloArray<IloNumVarArray> NumVarMatrix;
        typedef IloArray<NumVarMatrix> NumVarMatrix3;
        
        NumVarMatrix3 x(env,nn);
        NumVarMatrix y(env,nn);
         NumVarMatrix v(env,nn);
        for(i=0; i< nn; i++)
        {
            x[i] =  NumVarMatrix(env, M);
            y[i]=IloNumVarArray(env, M);
            v[i]=IloNumVarArray(env, T);
           
            for( j=0; j< M; j++)
            {
                  x[i][j]=IloNumVarArray(env,T);
                  y[i][j]=IloNumVar(env, 0.0, 1.0, ILOINT);
                  for( t=0;t<T;t++)
                    x[i][j][t]= IloNumVar(env, 0.0, 1.0, ILOINT);
            }
             for( t=0;t<T;t++)
                v[i][t]=IloNumVar(env, 0.0, 1.0, ILOINT);
                    
        }
     for(i=0;i<nn;i++)
    {
             IloExpr sumy(env);
            for(j=0;j<M;j++)
                sumy+=y[i][j];
            model.add(sumy==1);
    }

    
  
    for(i=0;i<nn;i++)
        for(j=0;j<M;j++)
        {
             IloExpr sumxi(env);
     
            for( t=jobs[sc][i].a;t<jobs[sc][i].d;t++)
                sumxi+=(x[i][j][t]);
              
            
            model.add(sumxi==y[i][j]*jobs[sc][i].l); 

             for( t=0;t<jobs[sc][i].a;t++)
                model.add(x[i][j][t]==0);


              for( t=jobs[sc][i].d;t<T;t++)
               model.add(x[i][j][t]==0);

        }

   


    for(j=0;j<M;j++)
            for( t=0;t<T;t++)
            {
                IloExpr sumxi(env);
                for(i=0; i< nn; i++)
                    sumxi+=(x[i][j][t]*jobs[sc][i].r);
                
                 model.add(sumxi<=centers[j].C); 

            }
              IloExpr sumxijt(env);
    for(i=0; i< nn; i++) 
            for(j=0;j<M;j++)
                for(int t=0;t<T;t++) 
                   model.add(v[i][t]>=x[i][j][t]);;
       
      
        for(i=0; i< nn; i++)
            for(int t=0;t<T;t++)  
                sumxijt+=(t*v[i][t]);
        model.add(zzz>=sumxijt);
      
        model.add(IloMinimize(env, zzz));
        
        cplex.solve();
      double u[M][T];
        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
           
            for(i=0; i< nn; i++) 
                for(j=0;j<M;j++)
                    for(int t=0;t<T;t++) 
                    zz+=centers[j].F[t]*centers[j].sigma*jobs[sc][i].r*(double)cplex.getValue(x[i][j][t])/centers[j].C;

            double load=0;
             for( t=0;t<T;t++)
            {
              for(j=0;j<M;j++)
                {
                    double load=0;
                    for(i=0; i< nn; i++)
                         load+=((double)cplex.getValue(x[i][j][t])*jobs[sc][i].r);
                     u[j][t]=load/centers[j].C;
                }
            }
            ofstream myfile;
    


    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/sto/u_black.txt");
   
    for( t=0;t<T;t++)
    {
        myfile<<endl<<t+1<<",";
        for( j=0;j<M;j++)
             myfile<<u[j][t]<<",";
                        
    }
 myfile.close();
            
           
                  
       }
       else
       cout<<"no "<<endl;
             
  
    


 }
       
    catch (IloException& e)
    {
        cerr << "C-Exp: " << e << endl;
    } catch (...)
    {
        cerr << "Unknown Exception" << endl;
    }
    
    
    
    
    env.end();
    return zz;
}

//------------------------------------------------
//------------------------------------------------
double MILP_EJD()
{
   int sc=0;
    double zz=-1;
    IloEnv env;
     double u[M][T];

    IloNumVar zzz(env, 0.0, IloInfinity);
    stringstream logfile;
    IloInt i,j,t;
   int nn=Ns[sc];
  //cout<<nn<<endl;
    try {
        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray con(env);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
        cplex.setOut(logfile);
        
        typedef IloArray<IloNumVarArray> NumVarMatrix;
        typedef IloArray<NumVarMatrix> NumVarMatrix3;
        
        NumVarMatrix3 x(env,nn);
       // NumVarMatrix y(env,nn);
         
        for(i=0; i< nn; i++)
        {
            x[i] =  NumVarMatrix(env, M);
          //  y[i]=IloNumVarArray(env, M);
            
           
            for( j=0; j< M; j++)
            {
                  x[i][j]=IloNumVarArray(env,T);
                //  y[i][j]=IloNumVar(env, 0.0, 1.0, ILOFLOAT);
                  for( t=0;t<T;t++)
                    x[i][j][t]= IloNumVar(env, 0.0, 1.0, ILOINT);
            }
            
                    
        }
    /* for(i=0;i<nn;i++)
    {
             IloExpr sumy(env);
            for(j=0;j<M;j++)
                sumy+=y[i][j];
            model.add(sumy==1);//relax
    }
*/
    
  
    for(i=0;i<nn;i++)
    {IloExpr sumxi(env);
        for(j=0;j<M;j++)
        {
            // IloExpr sumxi(env);
     
            for( t=jobs[sc][i].a;t<jobs[sc][i].d;t++)
                sumxi+=(x[i][j][t]);
              
            
          //  model.add(sumxi==y[i][j]*jobs[sc][i].l); 
         

             for( t=0;t<jobs[sc][i].a;t++)
                model.add(x[i][j][t]==0);


              for( t=jobs[sc][i].d;t<T;t++)
               model.add(x[i][j][t]==0);

        }
    model.add(sumxi==jobs[sc][i].l); 
    }
   


    for(j=0;j<M;j++)
            for( t=0;t<T;t++)
            {
                IloExpr sumxi(env);
                for(i=0; i< nn; i++)
                    sumxi+=(x[i][j][t]*jobs[sc][i].r);
                
                 model.add(sumxi<=centers[j].C); 

            }
              IloExpr sumxijt(env);
        for(i=0; i< nn; i++) 
        
            for(j=0;j<M;j++)
                for(int t=0;t<T;t++) 
                    sumxijt+=((centers[j].F[t]*centers[j].sigma*jobs[sc][i].r*x[i][j][t])/centers[j].C);
               
   
        model.add(zzz>=sumxijt);
      
        model.add(IloMinimize(env, zzz));
      
      
        cplex.solve();
      
        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
            zz=(double)cplex.getValue(zzz);
           
          
            double load=0;
             for( t=0;t<T;t++)
            {
              for(j=0;j<M;j++)
                {
                    double load=0;
                    for(i=0; i< nn; i++)
                         load+=((double)cplex.getValue(x[i][j][t])*jobs[sc][i].r);
                     u[j][t]=load/centers[j].C;
                }
            }
            ofstream myfile;
    


    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/sto/u_cplex.txt");
   
    for( t=0;t<T;t++)
    {
        myfile<<endl<<t+1<<",";
        for( j=0;j<M;j++)
             myfile<<u[j][t]<<",";
                        
    }
 myfile.close();
            
           
                  
       }
       else
       cout<<"no "<<endl;
  
 }
       
    catch (IloException& e)
    {
        cerr << "C-Exp: " << e << endl;
    } catch (...)
    {
        cerr << "Unknown Exception" << endl;
    }
    
    
    
    
    env.end();
    return zz;

}
//------------------------------------------------
//------------------------------------------------

double online(int mode)
{
    double zz=0;
    double u[M][T];
   int y[N];
    for(int i=0;i<N;i++)
    y[i]=-1;

    for(int sc=0; sc< num_sc; sc++)
    for(int i=0;i<Ns[sc];i++)
            remain[sc][i]=jobs[sc][i].l;
        
        
    
    if(mode==0)//our heuristic approach: sort jobs in increasing order of their deadline
        for(int sc=0; sc< num_sc; sc++)
            sort( jobs[sc], jobs[sc]+Ns[sc], &customer_sorter);
    
    

    for(int i=0;i<N;i++)//reset allocation variables for the current time slot
        for(int j=0;j<M;j++)
            for(int sc=0; sc< num_sc; sc++)
            {
                for(int t=0;t<T;t++)
                    bar_x[i][j][t]=0;
           
                alloc[i][j][sc]=0;
             }
         
    int sum=0;
    job tmp[N];
    
    for(int t=0;t<T;t++) //update workloads of scenarios by realization of sc=0 as the real scenario
    {
       for (int sc=0;sc<num_sc;sc++)
        {
            if(sc>0)
            {
                int cnt=0;
          
                for(int i=Ind[t][sc];i<Ns[sc];i++)//ind[t] index of the first jobs arrive after t
                    {
                        tmp[cnt]=jobs[sc][i];
                        cnt++;
                    }
          
                for(int i=0;i<Ind[t][0];i++)
                    jobs[sc][i]=jobs[0][i];
                
                int ii=Ind[t][0];
                for(int i=0;i<Ns[sc]-Ind[t][sc];i++)
                {
                    jobs[sc][ii]=tmp[i];
                    ii++;
                }

                 int offset=Ind[t][0]-Ind[t][sc];
                 Ns[sc]+=offset;
            
                int tmp=Ind[t][sc];
                Ind[t][sc]=Ind[t][0];;
                
                for(int k=t+1;k<TTT;k++)
                {
                    int diff= Ind[k][sc]-tmp;
                    tmp=Ind[k][sc];
                    Ind[k][sc]=Ind[k-1][sc]+diff;
                }
                 Ns[sc]=Ind[TTT-1][sc];
            }
        }
      for (int sc=0;sc<num_sc;sc++)
        {
            if(mode==0)
                greedy(t,sc,y,0);
            
            else
            {
                if(mode==1) //cplex
                 MILP_det(t,sc);
                 else
                 greedy(t,sc,y,1);//an online version for  MILP-QoS
            }

             
        }
      
        MILP_GAP(t);
    }


int n=Ns[0];
    for(int i=0; i< Ns[0]; i++) 
        for(int j=0;j<M;j++)
            for(int t=0;t<T;t++) 
                zz+=centers[j].F[t]*centers[j].sigma*jobs[0][i].r*bar_x[i][j][t]/centers[j].C;

     for(int t=0;t<T;t++)
        for(int j=0;j<M;j++)
        {
            double load=0;
            for(int i=0; i< n; i++)
                load+=(bar_x[i][j][t]*jobs[0][i].r);
                
            
            u[j][t]=load/centers[j].C;
         }
         
               
    ofstream myfile;
    
 
if(mode==1)
    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/sto/u_online.txt");
    else
    {
         if(mode==0)
             myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/sto/u_greedy.txt");
            else
            myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/sto/u_black.txt");
    }
    for( int t=0;t<T;t++)
    {
        myfile<<endl<<t+1<<",";
        for(int j=0;j<M;j++)
             myfile<<u[j][t]<<",";
                        
    }
    myfile.close();
    

return zz;

}
//------------------------------------------------
//------------------------------------------------

int main()
{
    clock_t begin, end;
    double elapsed_secs,elapsed_secs1,elapsed_secs2;
     
              
   generate_sample();
   
   show();

    begin = clock();
    double opt=MILP_EJD();
    end = clock();
    elapsed_secs1= 1000*double(end - begin) / CLOCKS_PER_SEC;
     
  
    begin = clock();
    double on_b=online(2);
    end = clock();
    elapsed_secs= 1000*double(end - begin) / CLOCKS_PER_SEC;
    
    begin = clock();
    double on=online(0);
    end = clock();
    elapsed_secs2= 1000*double(end - begin) / CLOCKS_PER_SEC;  

    std::cout<<endl<<100*(on-opt)/on<<" "<<100*(on_b-on)/on_b<<" "<<100*(on_b-opt)/on_b<<" "<<on<<" "<<on_b<<endl;
    std::cout<<elapsed_secs<<" "<<elapsed_secs1<<" "<<elapsed_secs2<<endl;

 
}