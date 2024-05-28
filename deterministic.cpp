#include <iostream>

#include <ilcplex/ilocplex.h>
#include <algorithm>
#include<iostream>
#include<list>
#include <random>
#include <omp.h>
using namespace std;
const int hours =4; //duration hours
const double rate =1;
const int slot_duration =1;//10 minutes 
const int T=hours*60/slot_duration;//48
int TTT=T-15;//schedul till a_i<100 (actual T is 100)
const int N=rate*15000;//maximum number of jobs in a scenario within T
 int Ns;
 int Ind[T];//ind[t]:index of the first job arrive after t
const int M=4;
const int C_range=rate*50;//low=104
const  double avalablae_rate= .0005;
const int sigma_range=20;

double alpha=10;
double ff[M];
struct center_av
{
    int id;
    double av_F;
};
bool customer_sorter1 (const center_av & s1, const center_av & s2)
    {
        return(s1.av_F <s2.av_F);
           
    }
struct job
{
    int id;
    int d;
    double l;
    double r;
    int a;
}jobs[N];
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
double bar_y[N][M];
double bar_u[N][M];
double bar_v[M][T];
double bar_w[N][M][T];
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
void generate_sample()/*workload generation*/
{
     ifstream myfile;
    ofstream myfile2;
    double fake=0;
   
    int N_sc[T];
    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/workload/det-Ns.txt");
    /*det-Ns gives the average number of submitted jobs per time slot*/
   
  
    for( int ii=0;ii<TTT;ii++)  
    {
             myfile>>N_sc[ii]; 
            N_sc[ii]*= rate;
    }
         
      for( int ii=TTT;ii<T;ii++)  
         N_sc[ii]=0;
        
    
    myfile.close();
  
    double ls[N];
    int rs[N];
    int tmp=0;
    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/workload/l-r.txt");
    /* l-r gives a random sample for the values of l and r*/
    
    for( int ii=0;ii<N;ii++)  
        myfile>>rs[ii]>>tmp>>ls[ii];  
       
    
       
    
    myfile.close();
//for each time slot define the charctristics of arriving workloads
        int i=0;
        Ns=0;
  
        for(int t=0;t<T;t++)
        {
         
           for( int ii=0;ii<N_sc[t];ii++)
            {
                int in=rand()%N;
                 jobs[i].l=ceil(ls[in]/(60*slot_duration))+rand()%3;
                jobs[i].l=min( jobs[i].l,10);
              
                 jobs[i].l=min(T-t, jobs[i].l);
                jobs[i].a=t;
                jobs[i].r=min(5,rs[in]);//1+rand()%7;
                jobs[i].d=min(T,jobs[i].a +jobs[i].l+alpha);
               i++; 
            }
           
            Ind[t]=i;
            Ns+=N_sc[t]; 
           
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
    ///*
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
    //*/
   //using random CI
   /*
   for(int t=0;t<60*hours;t++)//repeat each of them for 30 slots 
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
       
    }
  
}

//------------------------------------------------
//------------------------------------------------
double MILP_single(int i,int j) //find the optimal schedule for job i on data center j
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
     
            for( t=jobs[i].a;t<jobs[i].d;t++)
                sumxi+=(x[t]);
              
            
            
            model.add(sumxi==jobs[i].l); 
           
             for( t=0;t<jobs[i].a;t++)
                model.add(x[t]==0);


              for( t=jobs[i].d;t<T;t++)
               model.add(x[t]==0);

        
            for( t=jobs[i].a;t<jobs[i].d;t++)
            {
                 IloExpr sumxijt(env);
                 for(int ii=0;ii<i;ii++)
                    sumxijt+=(bar_x[ii][j][t]*jobs[ii].r);
            
                model.add(sumxijt+x[t]*jobs[i].r<=centers[j].C); 
     }

    IloExpr sumxijt(env);
       for(int t=0;t<T;t++) 
                    sumxijt+=((centers[j].F[t]*centers[j].sigma*jobs[i].r*x[t])/centers[j].C);
               
   
        model.add(zzz>=sumxijt);
      
        model.add(IloMinimize(env, zzz));
        
        cplex.solve();
       
        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
           
            zz= (double)cplex.getValue(zzz);
            for( t=jobs[i].a;t<jobs[i].d;t++)
                    bar_x[i][j][t]=cplex.getValue(x[t]);
                                
                  
       }
       else
       cout<<"infeasibke "<<endl;
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


double greedy_ILP()
{
    
    double result=0;
    center_av tmp[M];
     sort( jobs, jobs+Ns, &customer_sorter); //sort jobs in inreasing order of their deadline
    for(int i=0;i<Ns;i++)
    {
        for(int j=0;j<M;j++)
            for(int t=0;t<T;t++)
                bar_x[i][j][t]=0;
        
        int d=jobs[i].d-1;
        int a=max(0,jobs[i].a-1);
          for(int j=0;j<M;j++)
            {
                tmp[j].id=j;
                tmp[j].av_F=centers[j].av_F[d]-centers[j].av_F[a];
            }
                sort( tmp, tmp+M, &customer_sorter1);// sort datacenters based on their average carbon intensity.
               
        for(int j=0;j<M;j++)
        {
            result= MILP_single( i,tmp[j].id) ;  
        if(result>0)
        break; 
        }
            
    }


    result=0;
    for(int i=0; i< Ns; i++) 
        for(int j=0;j<M;j++)
            for(int t=0;t<T;t++) 
                result+=centers[j].F[t]*centers[j].sigma*jobs[i].r*bar_x[i][j][t]/centers[j].C;

    double u[M][T];
    for(int t=0;t<T;t++)
        for(int j=0;j<M;j++)
        {
            double load=0;
            for(int i=0; i< Ns; i++)
            {
                load+=(bar_x[i][j][t]*jobs[i].r);
                
            }
            u[j][t]=load/centers[j].C;
         }
    ofstream myfile;
    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/det/u_greedy.txt");
   
   
    for( int t=0;t<T;t++)
    {
        myfile<<endl<<t+1<<",";
        for(int j=0;j<M;j++)
             myfile<<u[j][t]<<",";
                        
    }
    myfile.close();
    
    return result;
}
//------------------------------------------------
//------------------------------------------------
//------------------------------------------------
//------------------------------------------------


int subproblem(double &obj,int nn)
{
    obj=-1;
   double zz=0;;
    IloEnv env;
    IloNumVar zzz(env,0.0, IloInfinity);
    stringstream logfile;
    IloInt i,j,t;
   
    
    try {
        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray con(env);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0);
        cplex.setParam(IloCplex::EpAGap, 0.0);
        cplex.setOut(logfile);
        cplex.setOut(env.getNullStream());

        typedef IloArray<IloNumVarArray> NumVarMatrix;
        typedef IloArray<NumVarMatrix> NumVarMatrix3;
        
        NumVarMatrix u(env,nn);
        NumVarMatrix v(env,M);
       NumVarMatrix3 w(env,nn);
        for(i=0; i<nn; i++)
        {
            u[i] =  IloNumVarArray(env,M);
            w[i]=NumVarMatrix(env, M);
           
            
            for( j=0; j< M; j++)
            {
                
                string cn="u"+to_string (i)+to_string (j);
                u[i][j]= IloNumVar(env,  -IloInfinity, IloInfinity, ILOFLOAT,cn.c_str());
                
                w[i][j]=IloNumVarArray(env,T);
               
                 
                for( t=0;t<T;t++)
                {
                    
                string cn="w"+to_string (i)+to_string (j)+to_string (t);
                
                    w[i][j][t]= IloNumVar(env,  -IloInfinity, 0, ILOFLOAT,cn.c_str());
                
                }
             }
        }

        for( j=0; j< M; j++)
        {
            v[j] =  IloNumVarArray(env,T);
            for( t=0;t<T;t++)
            {
                string cn="v"+to_string (j)+to_string (t);
                    v[j][t]= IloNumVar(env,  -IloInfinity,0, ILOFLOAT,cn.c_str());
            }
        }
        for(i=0; i< nn; i++)
           for( t=0;t<T;t++)
            for(j=0;j<M;j++)
            {
                if((t<jobs[i].d && t>=jobs[i].a))
                    model.add(u[i][j]+v[j][t]*jobs[i].r+w[i][j][t]<=(centers[j].sigma*centers[j].F[t]*jobs[i].r)/centers[j].C);
              // if(!(t<jobs[i].d && t>=jobs[i].a))
             else
             {
                     model.add(v[j][t]*jobs[i].r+w[i][j][t]<=(centers[j].sigma*centers[j].F[t]*jobs[i].r)/centers[j].C);
            model.add(u[i][j]==u[i][j]);
           
             } }   
 
                 
double rr=0;
    IloExpr sumxijt(env);
        for(i=0; i< nn; i++) 
            for(j=0;j<M;j++)
                sumxijt+=(bar_y[i][j]*jobs[i].l*u[i][j]);
            
         for(j=0;j<M;j++)
            for( t=0;t<T;t++) 
                    sumxijt+=(centers[j].C*v[j][t]);

        for(i=0; i< nn; i++) 
            for(j=0;j<M;j++) 
                for( t=0;t<T;t++) 
                    sumxijt+=  (w[i][j][t]);        
        model.add(zzz<=sumxijt);
      
        model.add(IloMaximize(env, zzz));
      //cplex.setParam(cplex.PreInd, 0); 
   
        cplex.solve();
    
         if (cplex.getStatus() == IloAlgorithm::Unbounded)
        {
            
            for(i=0; i< nn; i++)
                for(j=0;j<M;j++)
                    bar_u[i][j]=0;

            for(i=0; i< nn; i++)
                for(j=0;j<M;j++)
                    for( t=0;t<T;t++)        
                         bar_w[i][j][t]=0;

            for(j=0;j<M;j++)
                for( t=0;t<T;t++)        
                         bar_v[j][t]=0;             
            IloNumVarArray var(env);
            IloNumArray val(env);
           cplex.getRay(val, var);
          
           int cnt=0;
           for(i=0; i< nn; i++)
                for(j=0;j<M;j++)
                    for(int ii=0;ii<var.getSize();ii++)
                        if(var[ii].getId()==u[i][j].getId())
                            bar_u[i][j]=val[ii];
             
            for(i=0; i< nn; i++)
                for(j=0;j<M;j++)
                    for( t=0;t<T;t++)
                        for(int ii=0;ii<var.getSize();ii++)
                            if(var[ii].getId()==w[i][j][t].getId())
                                bar_w[i][j][t]=val[ii];
                             
            for(j=0;j<M;j++)
                for(t=0; t< T; t++)
                for(int ii=0;ii<var.getSize();ii++)
                    if(var[ii].getId()==v[j][t].getId())
                            bar_v[j][t]=val[ii];
           return 1;
        }
            
       else if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
          
            
            obj=cplex.getValue(zzz);
            for(i=0; i< nn; i++)
                for(j=0;j<M;j++)
                    for( t=0;t<T;t++)
                        bar_w[i][j][t]=cplex.getValue(w[i][j][t]);
                      
            for(i=0; i< nn; i++)
                for(j=0;j<M;j++)
                        bar_u[i][j]=cplex.getValue(u[i][j]);
                       
            for(j=0;j<M;j++)
             for(t=0; t< T; t++)
                    bar_v[j][t]=cplex.getValue(v[j][t]);
            
           zz= (double)cplex.getValue(zzz);
           
            return 2;
        }
         else
         {
            
             return 3;
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
   
    return 3;
}


//------------------------------------------------
//------------------------------------------------
double Bend()//Bender's decomposition method
{
    int nn=Ns;
    double LB=-INFINITY;
    double UB=INFINITY;
    double epsilon=1;//LP(nn);//UB/1000;
    
    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<M;j++)
           bar_y[i][j]=0;
       bar_y[i][rand()%M]=1;
    }
  
    IloInt i,j,t;
    double u[M][T];
    
    IloEnv env,env2;
    IloNumVar z(env, 0.0, IloInfinity);
    IloNumVar z2(env2, 0.0, IloInfinity);
     try {
        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray con(env);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
        cplex.setParam(IloCplex::EpAGap, 0.0);
        cplex.setOut(env.getNullStream());
        typedef IloArray<IloNumVarArray> NumVarMatrix;
       
        NumVarMatrix y(env,nn);
    
        for(i=0; i< nn; i++)
        {
            y[i]=IloNumVarArray(env, M);
          
            for( j=0; j< M; j++)
               y[i][j]=IloNumVar(env, 0.0, 1.0, ILOINT);
             
                 
       }
   
        double obj=0;int cnt=0;
     
        for(i=0;i<nn;i++)
        {
             IloExpr sumy(env);
            for(j=0;j<M;j++)
                sumy+=y[i][j];
            model.add(sumy==1);
        }

        model.add(z>=0);
        model.add(IloMinimize(env, z));
        
      
        while (UB-LB>.50)
        {
          
            int result=subproblem(obj,nn);
           
            if(result==1 ) //feasiblity cut
             {
                 
                 IloExpr sumyij(env);
                for(i=0;i<nn;i++)
                  for(j=0;j<M;j++)
                    sumyij+=(y[i][j]*jobs[i].l*bar_u[i][j]);

            
       
                 for(j=0;j<M;j++)
                    for(t=0;t<T;t++)
                        sumyij+=(centers[j].C*bar_v[j][t]);
                
            
                for(i=0;i<nn;i++)
                  for(j=0;j<M;j++)
                    for(t=0;t<T;t++)
                     sumyij+=(bar_w[i][j][t]);

                model.add(sumyij<=0);
                
            }
            else if(result==2)
          {
              
 
                IloExpr sumyij(env);
                for(i=0;i<nn;i++)
                    for(j=0;j<M;j++)
                     sumyij+=y[i][j]*jobs[i].l*bar_u[i][j];
            
              
                for(j=0;j<M;j++)
                     for(t=0;t<T;t++)
                         sumyij+=centers[j].C*bar_v[j][t];
 
           
            for(i=0;i<nn;i++)
                  for(j=0;j<M;j++)
                    for(t=0;t<T;t++)
                     sumyij+=bar_w[i][j][t];
             
           
            model.add(z>= sumyij);
            if(obj<UB)
                 UB=obj;

        }
        else
             return -1;
          
      
       
        cplex.solve();
        
        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {
           

            LB= (double)cplex.getValue(z);
           
            for( i=0;i<nn;i++)
                for( j=0;j<M;j++)
                {
                    bar_y[i][j]=cplex.getValue(y[i][j]);
                    if(cplex.getValue(y[i][j])>.85)
                        bar_y[i][j]=1;
                    else
                        bar_y[i][j]=0;
                }
        }
  // */
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

try {
        IloModel model(env2);
        IloNumVarArray vars(env2);
        IloRangeArray con(env2);
        IloCplex cplex(model);
        
        
        cplex.setParam(IloCplex::TiLim, 10000000.0);
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::EpGap, 0.0);
        cplex.setParam(IloCplex::EpAGap, 0.0);
        cplex.setOut(env2.getNullStream());
        typedef IloArray<IloNumVarArray> NumVarMatrix;
        typedef IloArray<NumVarMatrix> NumVarMatrix3;
     
       NumVarMatrix3 x(env2,nn);
        for(i=0; i< nn; i++)
        {
            x[i]=NumVarMatrix(env2, M);
            for( j=0; j< M; j++)
               { 
                   
                   x[i][j]=IloNumVarArray(env2, T);
                   for(t=0;t<T;t++)
                    x[i][j][t]=IloNumVar(env2, 0.0, 1.0, ILOFLOAT);
               }
       }
        for(i=0;i<nn;i++)
        for(j=0;j<M;j++)
        {
            IloExpr sumyij(env2);
            for(t=jobs[i].a;t<jobs[i].d;t++)
                sumyij+=x[i][j][t];
                model.add(sumyij==bar_y[i][j]*jobs[i].l);
       }
      
         for(j=0;j<M;j++)
            for(t=0;t<T;t++)
            {
                IloExpr sumyij(env2);
                for(i=0;i<nn;i++)
                  sumyij+=x[i][j][t]*jobs[i].r;
                model.add(sumyij<=centers[j].C);
            }
  
    
        IloExpr sumxijt(env2);
        for(i=0; i< nn; i++) 
            for(j=0;j<M;j++)
                for(int t=0;t<T;t++) 
                    sumxijt+=((centers[j].F[t]*centers[j].sigma*jobs[i].r*x[i][j][t])/centers[j].C);
               
     
           
          model.add(z2>=sumxijt);
        
        model.add(IloMinimize(env2, z2));
        
        cplex.solve();
 
        if (cplex.getStatus() == IloAlgorithm::Optimal)
        {

            
             for( t=0;t<T;t++)
            {
              for(j=0;j<M;j++)
                {

                    double load=0;
                    for(i=0; i< nn; i++)
                         load+=((double)cplex.getValue(x[i][j][t])*jobs[i].r);
                     u[j][t]=load/centers[j].C;
                }
            }
               
            ofstream myfile;
            myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/u_bend.txt");
            myfile<<endl<<"0"<<",";
           /* for( j=0;j<M;j++)
                myfile<<u[j][0]<<",";
     */
            for( t=0;t<T;t++)
            {
                myfile<<endl<<t+1<<",";
                for( j=0;j<M;j++)
                    myfile<<u[j][t]<<",";
            }
            myfile.close(); 
        }  
     

   }
   
    catch (IloException& e)
    {
        cerr << "C-Exp: " << e << endl;
    } catch (...)
    {
        cerr << "Unknown Exception" << endl;
    }
  env2.end();

 
 return  UB;
}


//------------------------------------------------
//------------------------------------------------
double MILP_QoS()
{
    int sc=0;
    double zz=-1;
    IloEnv env;
    
    IloNumVar zzz(env, 0.0, IloInfinity);
    stringstream logfile;
    IloInt i,j,t;
   int nn=Ns;
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
                  y[i][j]=IloNumVar(env, 0.0, 1.0, ILOFLOAT);
                  for( t=0;t<T;t++)
                    x[i][j][t]= IloNumVar(env, 0.0, 1.0, ILOFLOAT);
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
     
            for( t=jobs[i].a;t<jobs[i].d;t++)
                sumxi+=(x[i][j][t]);
              
            
            model.add(sumxi==y[i][j]*jobs[i].l); 

             for( t=0;t<jobs[i].a;t++)
                model.add(x[i][j][t]==0);


              for( t=jobs[i].d;t<T;t++)
               model.add(x[i][j][t]==0);

        }
    for(j=0;j<M;j++)
            for( t=0;t<T;t++)
            {
                IloExpr sumxi(env);
                for(i=0; i< nn; i++)
                    sumxi+=(x[i][j][t]*jobs[i].r);
                
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
                    zz+=centers[j].F[t]*centers[j].sigma*jobs[i].r*(double)cplex.getValue(x[i][j][t])/centers[j].C;

            double load=0;
             for( t=0;t<T;t++)
            {
              for(j=0;j<M;j++)
                {
                    double load=0;
                    for(i=0; i< nn; i++)
                         load+=((double)cplex.getValue(x[i][j][t])*jobs[i].r);
                     u[j][t]=load/centers[j].C;
                }
            }
            ofstream myfile;
    


    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/det/u_black.txt");
   
    for( t=0;t<T;t++)
    {
        myfile<<endl<<t+1<<",";
        for( j=0;j<M;j++)
             myfile<<u[j][t]<<",";
                        
    }
 myfile.close();
            
           
                  
       }
       else
       cout<<"infeasible "<<endl;
  
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
     double zz=-1;
    IloEnv env;
     double u[M][T];

    IloNumVar zzz(env, 0.0, IloInfinity);
    stringstream logfile;
    IloInt i,j,t;
   int nn=Ns;
  
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
                    x[i][j][t]= IloNumVar(env, 0.0, 1.0, ILOFLOAT);
            }
            
                    
        }
     for(i=0;i<nn;i++)
    {
             IloExpr sumy(env);
            for(j=0;j<M;j++)
                sumy+=y[i][j];
            model.add(sumy==1);//relax
    }
    for(i=0;i<nn;i++)
        for(j=0;j<M;j++)
        {
             IloExpr sumxi(env);
     
            for( t=jobs[i].a;t<jobs[i].d;t++)
                sumxi+=(x[i][j][t]);
              
            
            model.add(sumxi==y[i][j]*jobs[i].l); 

             for( t=0;t<jobs[i].a;t++)
                model.add(x[i][j][t]==0);


              for( t=jobs[i].d;t<T;t++)
               model.add(x[i][j][t]==0);

        }
    for(j=0;j<M;j++)
            for( t=0;t<T;t++)
            {
                IloExpr sumxi(env);
                for(i=0; i< nn; i++)
                    sumxi+=(x[i][j][t]*jobs[i].r);
                
                 model.add(sumxi<=centers[j].C); 

            }
              IloExpr sumxijt(env);
        for(i=0; i< nn; i++) 
        
            for(j=0;j<M;j++)
                for(int t=0;t<T;t++) 
                    sumxijt+=((centers[j].F[t]*centers[j].sigma*jobs[i].r*x[i][j][t])/centers[j].C);
               
   
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
                         load+=((double)cplex.getValue(x[i][j][t])*jobs[i].r);
                     u[j][t]=load/centers[j].C;
                }
            }
            ofstream myfile;
    


    myfile.open ("/Users/tayebehbahreini/Desktop/cpp/output/det/u_cplex.txt");
   
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

int main()
{
    clock_t begin, end;
    double elapsed_secs,elapsed_secs1,elapsed_secs2;
     
    generate_sample();
   /*compare greedy with MILP-EJD and MILP-QoS*/
    begin = clock();
    double gr=greedy_ILP();
    end = clock();
    elapsed_secs2= 1000*double(end - begin) / CLOCKS_PER_SEC;  

    
    begin = clock();
    double opt=MILP_EJD();
    end = clock();
    elapsed_secs1= 1000*double(end - begin) / CLOCKS_PER_SEC;
     
  
    begin = clock();
    double opt_q=MILP_QoS();
    end = clock();
    elapsed_secs= 1000*double(end - begin) / CLOCKS_PER_SEC;

    std::cout<<endl<<100*(opt_q-opt)/opt_q<<" "<<100*(opt_q-gr)/opt_q<<" "<<gr<<" "<<opt<<" "<<opt_q<<endl;
    std::cout<<elapsed_secs<<" "<<elapsed_secs1<<" "<<elapsed_secs2<<endl;

 
}
