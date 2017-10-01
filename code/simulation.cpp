#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

std::ofstream output;
using namespace std;
int iter = 100000;

double ran()
{
    return 1.0*rand()/RAND_MAX;
}

double* randn(double m, double n)
{
    //Given standard mean m, standard deviation n,
    //calculating 50000 value

    double *a;
    a = (double*)malloc(iter*sizeof(double));
    if(a == NULL)
	{
		printf("No enough memory\n");
		exit(0);
	}

    for(int i=0; i<iter; i++){
        //calculate normal distribution value based on mean m, and stdvar n.
        a[i] = m + n*(ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()-6);
    }
    return a;
}

void outResult(double **s){
    for(int i=0;i<13;i++)
    {
        sort(s[i], s[i]+iter);
        output<<s[i][iter/20]<<","<<s[i][iter/4]<<","<<s[i][iter/2]<<","<<s[i][iter*3/4]<<","<<s[i][iter*19/20]<<"\n";
    }
    output<<"\n";
}

double **storage()
{
    //Input 12 months mean and std deviation for W
    double m[12] = {0.9, 2.3, 5.0, 12.5, 9.0, -2.0, -4.5, -2.6, -1.5, 0.0, 4.0, 4.5};
    double n[12] = {3.6, 4.6, 5.4, 6.0, 5.0, 4.0, 3.5, 3.0, 3.0, 3.5, 6.0, 5.5};
    //Define mean and std deviation for D
    double f[iter], g[12];


    //Calculating Storage of 12 month
    double *w;
    double **s, **q;
    s = (double**)malloc(13*sizeof(double*));
    q = (double**)malloc(13*sizeof(double*));

    int k;
    for(k=0;k<13;++k){
        s[k] = (double*)malloc(iter*sizeof(double));
        q[k] = (double*)malloc(iter*sizeof(double));
    }

    //initiate first month beginning storage
    int i;
    for(i=0;i<iter;i++)
        s[0][i] = 3003.76;

    double d[iter];

    //calculate 12 months distribution
    for(i=0; i<12; i++)
    {
        //calculate std deviation of D
        g[i] = sqrt(pow(0.1,2.0)*(1-pow(-0.75,2.0)));

        //Calculate 50000 trails of W
        w = randn(m[i], n[i]);

        for(int j=0; j<iter; j++)
        {
            f[j] = 0.5+((-0.75)*0.1/n[i])*(w[j]-m[i]);
            //calculate D's value
            //mean + stdvar*(12*ran - 6)
            d[j] = f[j] + g[i]*(ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()-6);

            q[i][j] = 135.147*pow((733.73-0.00155*s[i][j]+50.63*log(s[i][j])-1131.388) ,1.679)*(365*24*3600/(12*pow(10, 9)));
            s[i+1][j] = w[j] - d[j] + s[i][j] - q[i][j];
        }
        free(w);
    }
    //free q's memory
    for(int r=0;r<13;++r){
        free(q[r]);
    }
    free(q);
    return s;
}

void outNew(double **s)
{
    for(int i=0;i<12;i++)
    {
        sort(s[i], s[i]+iter);
        output<<s[i][iter/20]<<","<<s[i][iter/4]<<","<<s[i][iter/2]<<","<<s[i][iter*3/4]<<","<<s[i][iter*19/20]<<"\n";
    }
    output<<"\n";
}

void freeNew(double **s)
{
    for(int i=0;i<12;++i){
        free(s[i]);
    }
    free(s);
}

double **lakeLevel(double **s)
{
    double **h;
    h = (double**)malloc(13*sizeof(double*));

    for(int k=0;k<13;++k){
        h[k] = (double*)malloc(iter*sizeof(double));
    }
    for(int i=0;i<iter;i++){
        h[0][i] = 1134.5;
    }
    for(int i=1;i<13;i++){
        for(int j=0;j<iter;j++){
            h[i][j] = 733.73-0.00155*s[i][j]+50.63*log(s[i][j]);
        }
    }
    return h;
}

double **outflow(double **h)
{
    double **q;
    q = (double**)malloc(12*sizeof(double*));

    for(int k=0;k<12;++k){
        q[k] = (double*)malloc(iter*sizeof(double));
    }

    for(int i=0;i<12;i++){
        for(int j=0;j<iter;j++){
            q[i][j] = 135.147*pow((h[i][j]-1131.388) ,1.679)*(365*24*3600/(12*pow(10, 9)));
        }
    }
    return q;
}

double **energy(double **h, double **q)
{
    double **e;
    e = (double**)malloc(12*sizeof(double*));

    for(int k=0;k<12;++k){
        e[k] = (double*)malloc(iter*sizeof(double));
    }
    for(int i=0;i<12;i++){
        for(int j=0;j<iter;j++){
            e[i][j] = 0.00706*0.94*(h[i][j]-1115)*q[i][j]*(pow(10,9)/(24*(365/12)*3600));
        }
    }
    return e;
}

void freeResult(double **s)
{
    for(int i=0;i<13;++i){
        free(s[i]);
    }
    free(s);
}

void outBelow(double **h)
{
    for(int i=0;i<13;i++)
    {
        int c=0;
        for(int j=0;j<iter;j++){
            if(h[i][j]<1133.5){
                c++;
            }
        }
        output<<c/iter<<"\n";
    }
    output<<"\n";
}

double **futureStorage(double **p, int year)
{
    double stdvar=0;
    if(year==0)
    {
        //Calculate the distribution of year 1
        output<<year+1<<" year Storage Distribution\n";
        outResult(p);
        for(int i=0;i<12;i++){
            double sum=0;
            for(int j=0;j<iter;j++){
                sum = sum + pow((p[i+1][j]-p[i+1][iter/2]),2);
            }
            stdvar = sqrt(sum/iter);
            output<<"Mean,"<<p[i+1][iter/2]<<",stdvar,"<<stdvar<<"\n";
        }
        output<<"\n";
        return p;
    }
    else
    {
        double** pre;
        pre = futureStorage(p, year-1);
        //Input 12 months mean and std deviation for W
        double m[12] = {0.9, 2.3, 5.0, 12.5, 9.0, -2.0, -4.5, -2.6, -1.5, 0.0, 4.0, 4.5};
        double n[12] = {3.6, 4.6, 5.4, 6.0, 5.0, 4.0, 3.5, 3.0, 3.0, 3.5, 6.0, 5.5};
        //Define mean and std deviation for D
        double f[iter], g[12];

        //Calculating Storage of 12 month
        double *w;
        double **s, **q;
        s = (double**)malloc(13*sizeof(double*));
        q = (double**)malloc(13*sizeof(double*));

        for(int k=0;k<13;++k){
            s[k] = (double*)malloc(iter*sizeof(double));
            q[k] = (double*)malloc(iter*sizeof(double));
        }
        //initiate first month beginning storage
        for(int i=0;i<iter;i++)
            s[0][i] = pre[12][i];

        double d[iter];

        //calculate 12 months distribution
        for(int i=0; i<12; i++)
        {
            //calculate std deviation of D
            g[i] = sqrt(pow(0.1,2.0)*(1-pow(-0.75,2.0)));

            //Calculate 50000 trails of W
            w = randn(m[i], n[i]);

            for(int j=0; j<iter; j++)
            {
                f[j] = 0.5+((-0.75)*0.1/n[i])*(w[j]-m[i]);

                d[j] = f[j] + g[i]*(ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()-6);

                q[i][j] = 135.147*pow((733.73-0.00155*s[i][j]+50.63*log(s[i][j])-1131.388) ,1.679)*(365*24*3600/(12*pow(10, 9)));
                s[i+1][j] = w[j] - d[j] + s[i][j] - q[i][j];
            }
            free(w);
        }

        //Calculate the distribution of year n
        output<<year+1<<" year Storage Distribution\n";
        outResult(s);
        for(int i=0;i<12;i++){
            double sum=0;
            for(int j=0;j<iter;j++){
                sum = sum + pow((s[i+1][j]-s[i+1][iter/2]),2);
            }
            stdvar = sqrt(sum/iter);
            output<<"Mean,"<<s[i+1][iter/2]<<",stdvar,"<<stdvar<<"\n";
        }
        output<<"\n";

        //free q's memory
        freeResult(q);
        //free pre's memory
        freeResult(pre);

        return s;
    }
}

double **netRev(double **e)
{
    double **r;
    r = (double**)malloc(12*sizeof(double*));

    for(int k=0;k<12;++k){
        r[k] = (double*)malloc(iter*sizeof(double));
    }
    for(int i=0;i<12;i++){
        for(int j=0;j<iter;j++){
            if(e[i][j]<220){
                r[i][j]=100000*220-250000*(220-e[i][j]);
            }
            else{
                r[i][j]=100000*220;
            }
        }
    }
    return r;
}

void outAnual(double *a)
{
    std::sort(a,a+iter);
    output<<a[iter/20]<<","<<a[iter/4]<<","<<a[iter/2]<<","<<a[iter*3/4]<<","<<a[iter*19/20]<<"\n\n";
}

double *anualRev(double **r)
{
    double *a;
    a = (double*)malloc(iter*sizeof(double));
    for(int i=0;i<iter;i++){
        a[i]=0;
        for(int j=0;j<12;j++){
            a[i] = a[i] + r[j][i];
        }
    }
    return a;
}

double **hInter(double **h)
{
    //Shoreline interest related to H
    double **t;
    t = (double**)malloc(12*sizeof(double*));

    for(int k=0;k<12;++k){
        t[k] = (double*)malloc(iter*sizeof(double));
    }
    for(int i=0;i<12;i++){
        for(int j=0;j<iter;j++){
            if(h[i][j]<1134){
               t[i][j] = 10*h[i][j]-11330;
            }
            if(h[i][j]>=1134 && h[i][j]<=1135){
                t[i][j] = 10;
            }
            else{
                t[i][j] = -20*h[i][j]+22710;
            }
        }
    }
    return t;
}

double **reStorage()
{
    double **m;
    m = (double**)malloc(12*sizeof(double*));

    for(int k=0;k<12;++k){
        m[k] = (double*)malloc(iter*sizeof(double));
        m[k] = (double*)malloc(iter*sizeof(double));
    }
    //Input 12 months mean and std deviation for W
    for(int j=0;j<iter;j++){
        m[0][j]=0.9; m[1][j]=2.3; m[2][j]=5.0; m[3][j]=12.5; m[4][j]=9.0; m[5][j]=-2.0;
        m[6][j]=-4.5; m[7][j]=-2.6; m[8][j]=-1.5; m[9][j]=0; m[10][j]=4; m[11][j]=4.5;
    }
    double n[12] = {3.6, 4.6, 5.4, 6.0, 5.0, 4.0, 3.5, 3.0, 3.0, 3.5, 6.0, 5.5};
    //Define mean and std deviation for D
    double f[iter], g[12];

    //Declare Storage and outflow of 12 month
    double **s, **q;
    s = (double**)malloc(13*sizeof(double*));
    q = (double**)malloc(13*sizeof(double*));
    for(int k=0;k<13;++k){
        s[k] = (double*)malloc(iter*sizeof(double));
        q[k] = (double*)malloc(iter*sizeof(double));
    }

    //initiate first month beginning storage
    int i;
    for(i=0;i<iter;i++)
        s[0][i] = 3003.76;

    //declare mean of D
    double d[iter];

    //calculate 12 months distribution
    for(i=0; i<12; i++)
    {
        //printf("hello ");
        //calculate std deviation of D
        g[i] = sqrt(pow(0.1,2.0)*(1-pow(-0.75,2.0)));
        n[i] = sqrt(pow(n[i],2.0)*(1-pow(0.75,2.0)));

        //Declare net inflow W for a month
        double *w;
        w = (double*)malloc(iter*sizeof(double));

        for(int j=0; j<iter; j++)
        {
            //Calculate iter trails of W
            w[j] = m[i][j] + n[i]*(ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()-6);
            if(i<11){
                m[i+1][j] = m[i+1][j]+0.75*(n[i+1]/n[i])*(w[j]-m[i][j]);
            }
            //Calculate mean of j th D
            f[j] = 0.5+((-0.75)*0.1/n[i])*(w[j]-m[i][j]);
            //calculate D's value
            //mean + stdvar*(12*ran - 6)
            d[j] = f[j] + g[i]*(ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()-6);

            q[i][j] = 135.147*pow((733.73-0.00155*s[i][j]+50.63*log(s[i][j])-1131.388) ,1.679)*(365*24*3600/(12*pow(10, 9)));
            s[i+1][j] = w[j] - d[j] + s[i][j] - q[i][j];
        }
        free(w);

    }
    //free q's memory
    freeResult(q);
    freeNew(m);
    return s;
}

double **refutureStorage(double **p, int year)
{
    double stdvar=0;
    if(year==0)
    {
        //Calculate the distribution of year 1
        output<<year+1<<" year Storage Distribution\n";
        outResult(p);
        for(int i=0;i<12;i++){
            double sum=0;
            for(int j=0;j<iter;j++){
                sum = sum + pow((p[i+1][j]-p[i+1][iter/2]),2);
            }
            stdvar = sqrt(sum/iter);
            output<<"Mean,"<<p[i+1][iter/2]<<",stdvar,"<<stdvar<<"\n";
        }
        output<<"\n";
        return p;
    }
    else
    {
        double** pre;
        pre = futureStorage(p, year-1);

        double **m;
        m = (double**)malloc(12*sizeof(double*));

        for(int k=0;k<12;++k){
            m[k] = (double*)malloc(iter*sizeof(double));
            m[k] = (double*)malloc(iter*sizeof(double));
        }
        //Input 12 months mean and std deviation for W
        for(int j=0;j<iter;j++){
            m[0][j]=0.9; m[1][j]=2.3; m[2][j]=5.0; m[3][j]=12.5; m[4][j]=9.0; m[5][j]=-2.0;
            m[6][j]=-4.5; m[7][j]=-2.6; m[8][j]=-1.5; m[9][j]=0; m[10][j]=4; m[11][j]=4.5;
        }
        double n[12] = {3.6, 4.6, 5.4, 6.0, 5.0, 4.0, 3.5, 3.0, 3.0, 3.5, 6.0, 5.5};
        //Define mean and std deviation for D
        double f[iter], g[12];

        //Declare Storage and outflow of 12 month
        double **s, **q;
        s = (double**)malloc(13*sizeof(double*));
        q = (double**)malloc(13*sizeof(double*));
        for(int k=0;k<13;++k){
            s[k] = (double*)malloc(iter*sizeof(double));
            q[k] = (double*)malloc(iter*sizeof(double));
        }

        //initiate first month beginning storage
        //initiate first month beginning storage
        for(int i=0;i<iter;i++)
            s[0][i] = pre[12][i];

        //declare mean of D
        double d[iter];

        //calculate 12 months distribution
        for(int i=0; i<12; i++)
        {
            //calculate std deviation of D
            g[i] = sqrt(pow(0.1,2.0)*(1-pow(-0.75,2.0)));
            n[i] = sqrt(pow(n[i],2.0)*(1-pow(0.75,2.0)));

            //Declare net inflow W for a month
            double *w;
            w = (double*)malloc(iter*sizeof(double));

            for(int j=0; j<iter; j++)
            {
                //Calculate iter trails of W
                w[j] = m[i][j] + n[i]*(ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()-6);
                if(i<11){
                    m[i+1][j] = m[i+1][j]+0.75*(n[i+1]/n[i])*(w[j]-m[i][j]);
                }
                //Calculate mean of j th D
                f[j] = 0.5+((-0.75)*0.1/n[i])*(w[j]-m[i][j]);
                //calculate D's value
                //mean + stdvar*(12*ran - 6)
                d[j] = f[j] + g[i]*(ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()+ran()-6);

                q[i][j] = 135.147*pow((733.73-0.00155*s[i][j]+50.63*log(s[i][j])-1131.388) ,1.679)*(365*24*3600/(12*pow(10, 9)));
                s[i+1][j] = w[j] - d[j] + s[i][j] - q[i][j];
            }
            free(w);
        }

        //Calculate the distribution of year n
        output<<year+1<<" year Storage Distribution\n";
        outResult(s);
        for(int i=0;i<12;i++){
            double sum=0;
            for(int j=0;j<iter;j++){
                sum = sum + pow((s[i+1][j]-s[i+1][iter/2]),2);
            }
            stdvar = sqrt(sum/iter);
            output<<"Mean,"<<s[i+1][iter/2]<<",stdvar,"<<stdvar<<"\n";
        }
        output<<"\n";

        //free q's memory
        freeResult(q);
        freeNew(m);
        //free pre's memory
        freeResult(pre);

        return s;
    }
}

int main()
{
    output.open("distribution.txt");
    srand(time(NULL));

    //calculate
    double **s=storage();
    double **h=lakeLevel(s);
    double **q=outflow(h);
    double **e=energy(h, q);
    double **r=netRev(e);
    double **t=hInter(h);
    double *a = anualRev(r);
    double *b = anualRev(t);
    //recalculate
    double **rs = reStorage();
    double **rh = lakeLevel(rs);
    double **rq=outflow(rh);
    double **re=energy(rh, rq);
    double **rr=netRev(re);
    double **rt=hInter(rh);
    double *ra= anualRev(rr);
    double *rb= anualRev(rt);

    //output
    output<<"S"<<"\n";
    outResult(s);
    output<<"H"<<"\n";
    outResult(h);
    output<<"Below1133.5"<<"\n";
    outBelow(h);
    output<<"Q"<<"\n";
    outNew(q);
    output<<"E"<<"\n";
    outNew(e);
    output<<"Rev"<<"\n";
    outNew(r);
    output<<"Interest"<<"\n";
    outNew(t);
    output<<"Annual Rev\n";
    outAnual(a);
    output<<"Annual Interest\n";
    outAnual(b);
    //Recalculate output
    output<<"Re storage\n";
    outResult(rs);
    output<<"ReH"<<"\n";
    outResult(rh);
    output<<"ReBelow1133.5"<<"\n";
    outBelow(rh);
    output<<"ReQ"<<"\n";
    outNew(rq);
    output<<"ReE"<<"\n";
    outNew(re);
    output<<"ReRev"<<"\n";
    outNew(rr);
    output<<"ReInterest"<<"\n";
    outNew(rt);
    output<<"ReAnnual Rev\n";
    outAnual(ra);
    output<<"ReAnnual Interest\n";
    outAnual(rb);

    //free memory 1st calculation
    free(a);
    free(b);
    freeNew(t);
    freeNew(r);
    freeNew(q);
    freeNew(e);
    freeResult(h);
    //free memory recalculation
    free(ra);
    free(rb);
    freeNew(rt);
    freeNew(rr);
    freeNew(rq);
    freeNew(re);
    freeResult(rh);

    output<<"Future Storage\n";
    double **futures = futureStorage(s, 20);
    freeResult(futures);

    output<<"Re Future Storage\n";
    double **refutures = refutureStorage(rs, 20);
    freeResult(refutures);

    output.close();
    return 0;
}
