//
//  main.cpp
//  Test Project
//
//  Created by Zixiang  Pei on 5/5/19.
//  Copyright © 2019 Zixiang  Pei. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

float lumda = 0;
float mpc =0.46;
double e=2.7182818284590452353602875;
double a1=-35.0/6;
double a2=29.0/6;
double a3=1.0;
double k1=-35.0/6;
double k2=41.0/6;
double k3=0.0;
const int rowcount=6;
const int colcount=6;
double SP_inteval =1/5.0;
double SC_inteval =1/5.0;
double delta=0.002;
rowvec P_strategy(rowcount);
colvec C_strategy(colcount);

struct entry
{
public:
    double UP;
    double UC;
    double K;
    int year;
    //entry son[6][6];
};

typedef Mat<entry> mat_entry;

double calculate_UP(double p, double c)
{
    double F=c*(a1*pow(p, 2)+a2*p+a3);
    //double F=(a1*pow(p, 2)+a2*p+a3);
    double G=pow(e, p*c-0.5)/((pow(e, p*c-0.5))+1);
    return (0.5+lumda)*F+(0.5-lumda)*G;
}

double calculate_UC(double p,double c)
{
    double H=p*(k1*pow(c, 2)+k2*c+k3);
    //double H=(k1*pow(c, 2)+k2*c+k3);
    //cout << H << "hhh"<< endl;
    double I;
    double temp=(1-c)/(1-mpc);
    if(1<temp)
        I=1;
    else
        I=temp;
    //cout << I << "ppp"<< endl;
    return 0.7*H+0.3*I;
}

struct coord
{
    int x;
    int y;
};

struct utilities
{
    double x;
    double y;
};

double calculate_collaboration_factor(double p,double c,double year_from_start)
{
    double max,min;
    if(p>c)
    {
        max=p;
        min=c;
    }
    else
    {
        max=c;
        min=p;
    }
    if((min==0)&&(max==0))
    {
        //cout<<"hhh"<<endl;
        return 0;
    }
    double temp = cos(atan(1)-atan(min/max));
    // double temp2 = double((5.0/year_from_start)+1)*2*sqrt(pow(p, 2)+pow(c, 2))/sqrt(2);
    double temp2 = double(2)*2*sqrt(pow(p, 2)+pow(c, 2))/sqrt(2);
    return temp*temp2;
}

//double P_calculate_maxmin(entry** matrix)
//{
//    coord temp[rowcount];
//    for(int i=0; i<rowcount;i++)
//    {
//        int min =100000;
//        for (int j=0; j<colcount;j++)
//        {
//            if(min>matrix[i][j].UP)
//            {
//                min=matrix[i][j].UP;
//                temp[i].x=i*SP_inteval;
//                temp[i].y=j*SC_inteval;
//            }
//        }
//    }
//    int max =0;
//    coord track;
//    for(int k=0;k<rowcount;k++)
//    {
//        if(max<matrix[temp[k].x][temp[k].y].UP)
//        {
//            max=matrix[temp[k].x][temp[k].y].UP;
//            track = temp[k];
//        }
//    }
//    return track.x;
//}
//
//double C_calculate_maxmin(entry** matrix)
//{
//    coord temp[colcount];
//    for(int i=0; i<colcount;i++)
//    {
//        int min =5;
//        for (int j=0; j<rowcount;j++)
//        {
//            if(min>matrix[j][i].UC)
//            {
//                min=matrix[j][i].UC;
//                temp[i].x=j*SP_inteval;
//                temp[i].y=i*SC_inteval;
//            }
//        }
//    }
//    int max =0;
//    coord track;
//    for(int k=0;k<colcount;k++)
//    {
//        if(max<matrix[temp[k].x][temp[k].y].UC)
//        {
//            max=matrix[temp[k].x][temp[k].y].UC;
//            track = temp[k];
//        }
//    }
//    return track.y;
//}



utilities calculate_optimal(mat_entry matrix)
{
    mat UP(rowcount,colcount);
    mat UC(rowcount,colcount);
    for(int i =0; i< rowcount;i++)
        for(int j = 0; j< colcount;j++)
        {
            UP(i,j)=matrix(i,j).UP;
            UC(i,j)=matrix(i,j).UC;
        }
    const double INITIAL = 1.0/rowcount;
    rowvec NE_P(rowcount);
    colvec NE_C(colcount);
    for(int i=0;i<rowcount;i++)
        NE_P(i)=INITIAL;
    for(int i=0;i<colcount;i++)
        NE_C(i)=INITIAL;

    while(true)
    {
        rowvec temp1(NE_P);
        colvec temp2(NE_C);
        double CC[rowcount],DD[colcount];
        for(int i =0; i< rowcount;i++)
        {
            double x= as_scalar(UP.row(i)*NE_C-NE_P*UP*NE_C);
            CC[i]=0;
            if(x>0)
                CC[i]=x;
        }
        for(int j=0; j<colcount;j++)
        {
            double y = as_scalar(NE_P*UC.col(j)-NE_P*UC*NE_C);
            DD[j]=0;
            if(y>0)
                DD[j]=y;
        }
        double sum_C=0, sum_D=0;
        for(int i=0;i< rowcount;i++)
            sum_C+=CC[i];
        for(int j=0;j< colcount;j++)
            sum_D+=DD[j];
        for(int i=0;i<rowcount;i++)
            NE_P(i)=(NE_P(i)+CC[i])/(1+sum_C);
        for(int j=0;j<colcount;j++)
            NE_C(j)=(NE_C(j)+DD[j])/(1+sum_D);
        double rate_of_change1=0,rate_of_change2=0;
        for(int i =0; i< rowcount;i++)
            rate_of_change1+=pow((NE_P(i)-temp1(i)),2);
        for(int j =0; j< colcount;j++)
            rate_of_change2+=pow((NE_C(j)-temp2(j)),2);
        if(rate_of_change1<0.001 && rate_of_change2<0.001)
            break;
    }

    utilities UP_UC;
    UP_UC.x = double(as_scalar(NE_P*UP*NE_C));
    UP_UC.y = double(as_scalar(NE_P*UC*NE_C));
    P_strategy = NE_P;
    C_strategy = NE_C;
    //P_strategy.print("PPPP: ");
    //C_strategy.print("CCCC: ");
    return UP_UC;
}

//{
    //    PMatrix Pmat;
    //    CMatrix Cmat;
    //    VectorXf NE_P(rowcount);
    //    VectorXf NE_C(colcount);
    //    VectorXf temp(rowcount+colcount),result(rowcount+colcount);
    //    for(int i=0;i<rowcount;i++)
    //    {
    //        NE_P(i)=1/6;
    //        NE_C(i)=1/6;
    //    }
    //    for(int i=0;i<rowcount;i++)
    //    {
    //        for(int j=0;j<colcount;j++)
    //        {
    //            Pmat(i,j)=matrix[i][j].UP;
    //            Cmat(i,j)=matrix[i][j].UC;
    //        }
    //    }
    //    VectorXf CC(rowcount), DD(colcount);
    //    double c,c1,c2,d,d1,d2;
    //    while(true)
    //    {
    //
    ////        c=0; d=0;
    ////        c1= Pmat.row(Jerry)*NE_C;
    ////        c2= NE_P.transpose()*Pmat*NE_C;
    ////        if((c1-c2)>0)
    ////            c=c1-c2;
    ////        CC(Jerry) =c;
    //        double sumC=0;
    ////        for(int i=0; i<rowcount;i++)
    ////            sumC+=CC[i];
    ////        NE_P(Jerry)=(NE_P(Jerry)+c)/(1+sumC);
    ////
    ////        d1= NE_P.transpose()*Cmat.row(Joe);
    ////        d2= NE_P.transpose()*Cmat*NE_C;
    ////        if((d1-d2)>0)
    ////            d=d1-d2;
    ////        DD(Joe)=d;
    ////        double sumD=0;
    ////        for(int j=0; j<colcount;j++)
    ////            sumD+=DD[j];
    ////        NE_C(Jerry)=(NE_C(Jerry)+d)/(1+sumD);
    ////
    ////        for(int i=0;i< rowcount;i++)
    ////            result(i)=NE_P(i);
    ////
    ////        for(int j=rowcount;j<(rowcount+colcount);j++)
    ////            result[j]=NE_P(j);
    ////
    //        double sum=0;
    ////        for (int k=0;k< (rowcount+colcount);k++)
    ////            sum+=pow((result[k]-temp[k]), 2);
    ////        if (sum<=0.01)
    ////            break;
    ////        for (int l=0; l<(rowcount+colcount);l++)
    ////            temp(l)=result(l);
    ////        Jerry++;
    ////    }
    //        for(int k=0; k<rowcount;c++)
    //        {
    //            c=0;
    //            c1= Pmat.row(k)*NE_C;
    //            c2= NE_P.transpose()*Pmat*NE_C;
    //            if((c1-c2)>0)
    //                c=c1-c2;
    //            CC[k]=c;
    //        }
    //
    //        for(int i=0; i<rowcount;i++)
    //            sumC+=CC[i];
    //        for(int j=0;j<rowcount;j++)
    //            NE_P(j)=(NE_P(j)+c)/(1+sumC);
    //
    //        for(int k=0; k<rowcount;c++)
    //        {
    //            d=0;
    //            d1= NE_P.transpose()*Cmat.row(k);
    //            d2= NE_P.transpose()*Cmat*NE_C;
    //            if((d1-d2)>0)
    //                d=d1-d2;
    //            DD[k]=d;
    //        }
    //
    //        for(int i=0; i<rowcount;i++)
    //            sumD+=DD[i];
    //        for(int j=0;j<rowcount;j++)
    //            NE_C(j)=(NE_C(j)+c)/(1+sumD);
    //
    //        for(int i=0;i< rowcount;i++)
    //            result(i)=NE_P(i);
    //        for(int j=rowcount;j<(rowcount+colcount);j++)
    //            result[j]=NE_P(j);
    //
    //        double sum=0;
    //        for (int k=0;k< (rowcount+colcount);k++)
    //            sum+=pow((result[k]-temp[k]), 2);
    //        if (sum<=0.01)
    //            break;
    //        for (int l=0; l<(rowcount+colcount);l++)
    //            temp(l)=result(l);
    //    }
    //
    //    utilities UP_UC;
    //    UP_UC.x=NE_P.transpose()*Pmat*NE_C;
    //    UP_UC.y=NE_P.transpose()*Cmat*NE_C;
    //
    //    return UP_UC;
//}

//Cvector C_calculate_optimal(entry** matrix)
//{
//
//}

int main(int argc, const char * argv[]) {

//    entry** initial =new entry*[rowcount];
//    for(int i=0; i<rowcount;i++)
//        initial[i]=new entry[colcount];

    mat initial_UP(rowcount,colcount);
    mat initial_UC(rowcount,colcount);
    mat initial_K(rowcount,colcount);

    mat_entry initial(rowcount,colcount);
    entry**first_handler[rowcount*colcount];
    for(int i=0; i<rowcount*colcount;i++)
    {
        first_handler[i]=new entry*[rowcount];
        for (int j=0;j<rowcount;j++)
            first_handler[i][j]=new entry[colcount];
    }

//    cube first_handler_UP(rowcount*colcount,rowcount,colcount);
//    cube first_handler_UC(rowcount*colcount,rowcount,colcount);
//    cube first_handler_K(rowcount*colcount,rowcount,colcount);

    //cube first_handler(rowcount*colcount,rowcount,colcount,entry);
    mat_entry first_handler[rowcount*colcount];
    entry**second_handler[int(pow(rowcount*colcount,2))];
    for(int i=0; i<int(pow(rowcount*colcount,2));i++)
    {
        second_handler[i]=new entry*[rowcount];
        for (int j=0;j<rowcount;j++)
            second_handler[i][j]=new entry[colcount];
    }
//    cube second_handler_UP(int(pow(rowcount*colcount,2)),rowcount,colcount);
//    cube second_handler_UC(int(pow(rowcount*colcount,2)),rowcount,colcount);
//    cube second_handler_K(int(pow(rowcount*colcount,2)),rowcount,colcount);
//    cube second_handler(int(pow(rowcount*colcount,2)),rowcount,colcount,entry);

    mat_entry second_handler[int(pow(rowcount*colcount,2))];

//    entry**third_handler[int(pow(rowcount*colcount,3))];
//    for(int i=0; i<int(pow(rowcount*colcount,3));i++)
//    {
//        third_handler[i]=new entry*[rowcount];
//        for (int j=0;j<rowcount;j++)
//            third_handler[i][j]=new entry[colcount];
//    }
//initial creation
    //entry***first_handler=new entry**[rowcount*colcount];
//    cube third_handler_UP(int(pow(rowcount*colcount,3)),rowcount,colcount);
//    cube third_handler_UC(int(pow(rowcount*colcount,3)),rowcount,colcount);
//    cube third_handler_K(int(pow(rowcount*colcount,3)),rowcount,colcount);
//    cube third_handler(int(pow(rowcount*colcount,2)),rowcount,colcount,entry);

    mat_entry third_handler[int(pow(rowcount*colcount,3))];

    for(int i=0;i<rowcount;i++)
        for(int j=0;j<colcount;j++)
        {
            initial(i,j).UP=calculate_UP(i*SP_inteval, j*SC_inteval);
            initial(i,j).UC=calculate_UC(i*SP_inteval, j*SC_inteval);
            initial(i,j).year=2009;
            initial(i,j).K=calculate_collaboration_factor(double(i*SP_inteval), double(j*SC_inteval), initial(i,j).year-2008);
        }

//first_handler creation
    int counter=0;
    for(int i=0; i<rowcount;i++)
        for(int j=0;j<colcount;j++)
        {
            for(int l=0;l<rowcount;l++)
            {
                for(int m=0;m<colcount;m++)
                {
                    first_handler[counter](l,m).UP=initial(i,j).K*initial(l,m).UP;
                    first_handler[counter](l,m).UC=initial(i,j).K*initial(l,m).UP;
                    first_handler[counter](l,m).year=initial(i,j).year+1;
                    first_handler[counter](l,m).K=calculate_collaboration_factor(l*SP_inteval, m*SC_inteval, first_handler[counter](l,m).year-2008);
                }
            }
            counter++;
        }

    counter =0;
    //for(int i=0;i<int(pow(rowcount*colcount,2));i++)
//second_handler creation
    for(int count=0; count<rowcount*colcount;count++)
        for(int j=0; j<rowcount;j++)
            for(int k=0; k<colcount; k++)
            {
                for(int l=0;l<rowcount;l++)
                    for(int m=0;m<colcount;m++)
                    {
                        second_handler[counter](l,m).UP=first_handler[count](j,k).K*first_handler[count](l,m).UP;
                        second_handler[counter](l,m).UC=first_handler[count](j,k).K*first_handler[count](l,m).UC;
                        second_handler[counter](l,m).year=first_handler[count](l,m).year+1;
                        second_handler[counter](l,m).K=calculate_collaboration_factor(l*SP_inteval, m*SC_inteval,second_handler[counter](l,m).year-2008);
                    }
                counter++;
            }
//third_handler creation
    counter=0;
    for(int count=0; count<int(pow(rowcount*colcount,2));count++)
        for(int j=0; j<rowcount;j++)
            for(int k=0; k<colcount; k++)
            {
                for(int l=0;l<rowcount;l++)
                    for(int m=0;m<colcount;m++)
                    {
                        third_handler[counter](l,m).UP=second_handler[count](j,k).K*second_handler[count](l,m).UP;
                        third_handler[counter](l,m).UC=second_handler[count](j,k).K*second_handler[count](l,m).UC;
                        third_handler[counter](l,m).year=second_handler[count](l,m).year+1;
                        third_handler[counter](l,m).K=calculate_collaboration_factor(l*SP_inteval, m*SC_inteval, third_handler[counter](l,m).year-2008);
                    }
                counter++;
            }

//    double Third_Umax_X[(int(pow(rowcount*colcount,3)))];
//    double Third_Umax_Y[(int(pow(rowcount*colcount,3)))];

    utilities Third_Umax[(int(pow(rowcount*colcount,3)))];

    //vector<double> Third_Umax_C;
    for(int i=0; i< int(pow(rowcount*colcount,3));i++)
    {
        //        int p= P_calculate_maxmin(third_handler[i]);
        //        int c= C_calculate_maxmin(third_handler[i]);
        //        Third_Umax_P[i]=third_handler[i][p][c].UP;
        //        Third_Umax_C[i]=third_handler[i][p][c].UC;
        mat_entry temp(rowcount,colcount);
        for(int m=0; m<rowcount;m++)
            for(int n=0;n<colcount;n++)
                temp(m,n) = third_handler[i](m,n);

        Third_Umax[i]=calculate_optimal(temp);
        //Third_Umax_P.push_back(third_handler[i][p][c].UP);
        //Third_Umax_C.push_back(third_handler[i][p][c].UC);
    }

    //    vector<double> Second_Umax_P;
    //    vector<double> Second_Umax_C;
//    double Second_Umax_P[(int(pow(rowcount*colcount,2)))];
//    double Second_Umax_C[(int(pow(rowcount*colcount,2)))];
    utilities Second_Umax[(int(pow(rowcount*colcount,2)))];
    counter=0;
    for(int i=0; i<int(pow(rowcount*colcount, 2));i++)
    {
        for(int j=0; j<rowcount;j++)
            for(int k=0;k<colcount;k++)
            {
                second_handler[i](j,k).UP=second_handler[i](j,k).UP+delta*Third_Umax[counter].x;
                second_handler[i](j,k).UC=second_handler[i](j,k).UC+delta*Third_Umax[counter].y;
                counter++;
            }
    }

    for(int i=0; i< int(pow(rowcount*colcount,2));i++)
    {
        //        int p= P_calculate_maxmin(second_handler[i]);
        //        int c= C_calculate_maxmin(second_handler[i]);
        //        Second_Umax_P[i]=second_handler[i][p][c].UP;
        //        Second_Umax_C[i]=second_handler[i][p][c].UC;
        mat_entry temp(rowcount,colcount);
        for(int m=0; m<rowcount;m++)
            for(int n=0;n<colcount;n++)
                temp(m,n) = third_handler[i](m,n);
        Second_Umax[i]=calculate_optimal(temp);
    }

//    double First_Umax_P[(int(pow(rowcount*colcount,1)))];
//    double First_Umax_C[(int(pow(rowcount*colcount,1)))];
    utilities First_Umax[(int(pow(rowcount*colcount,1)))];
    double final_c,final_p;
    counter=0;
    for(int i=0; i<int(pow(rowcount*colcount,1));i++)
    {
        for(int j=0; j<rowcount;j++)
            for(int k=0;k<colcount;k++)
            {
                first_handler[i](j,k).UP=first_handler[i](j,k).UP+delta*Second_Umax[counter].x;
                first_handler[i](j,k).UC=first_handler[i](j,k).UC+delta*Second_Umax[counter].y;
                counter++;
            }
    }

    for(int i=0; i< int(pow(rowcount*colcount,1));i++)
    {
        //        int p= P_calculate_maxmin(first_handler[i]);
        //        int c= C_calculate_maxmin(first_handler[i]);
        //        First_Umax_P[i]=first_handler[i][p][c].UP;
        //        First_Umax_C[i]=first_handler[i][p][c].UC;
        mat_entry temp(rowcount,colcount);
        for(int m=0; m<rowcount;m++)
            for(int n=0;n<colcount;n++)
                temp(m,n) = third_handler[i](m,n);
        First_Umax[i]=calculate_optimal(temp);
    }

    counter=0;
    for(int i=0; i<rowcount;i++)
        for(int j=0; j<colcount;j++)
        {
            initial(i,j).UP=initial(i,j).UP+delta*First_Umax[counter].x;
            initial(i,j).UC=initial(i,j).UC+delta*First_Umax[counter].y;
            counter++;
        }

    calculate_optimal(initial);

    //    int final_p=P_calculate_maxmin(initial);
    //    int final_c=C_calculate_maxmin(initial);
//    for(int i=0; i<rowcount;i++)
//        for(int j=0;j<colcount;j++)
//        {
//            initial[i][j].K=calculate_collaboration_factor(final_p, final_c, initial[i][j].year-2008);
//        }

    double P_findK=0,C_findK=0;
    for(int i=0; i< rowcount; i++)
    {
        P_findK+= P_strategy(i)*i*SP_inteval;
    }
    for(int j=0; j<colcount; j++)
    {
        C_findK+= C_strategy(j)*j*SC_inteval;
    }

    double the_K = calculate_collaboration_factor(P_findK,C_findK,initial(0,0).year-2008);

    //entry next_year[rowcount][colcount];

    mat_entry next_year(rowcount,colcount);

    for(int i=0;i<rowcount;i++)
    {
        for(int j=0;j<rowcount;j++)
        {
            next_year(i,j).UP=initial(i,j).K*initial(i,j).UP;
            next_year(i,j).UC=initial(i,j).K*initial(i,j).UC;
            next_year(i,j).year=initial(i,j).year+1;
            next_year(i,j).K=calculate_collaboration_factor(i*SP_inteval, j*SC_inteval, next_year(i,j).year-2008);
        }
    }

    mat P_initial(rowcount,colcount);
    mat C_initial(rowcount,colcount);
    for(int i=0; i< rowcount; i++)
        for(int j=0; j<colcount;j++)
        {
            P_initial(i,j) = initial(i,j).UP;
            C_initial(i,j) = initial(i,j).UC;
        }

    mat P_next(rowcount,colcount), C_next(rowcount,colcount);
    P_next = the_K*P_initial;
    C_next= the_K*C_initial;

    P_initial.print("P_initial ");
    C_initial.print("C_initial ");

    P_next.print("P_next: " );
    C_next.print("C_next: " );

    //    cout << final_p << " " << final_c<< endl;
    //
    //    for(int i=0;i<rowcount;i++)
    //    {
    //        for(int j=0;j<colcount;j++)
    //        {
    //            cout<<initial[i][j].UC<<" ";
    //        }
    //        cout <<endl;
    //    }
    cout <<"---------------------------"<<endl;
    //    for(int what=0; what <46656; what++)
    //    {
    //        for(int i=0;i<rowcount;i++)
    //            {
    //                for(int j=0;j<colcount;j++)
    //                {
    //                    cout<<third_handler[what][i][j].UC<<" ";
    //                }
    //                cout <<endl;
    //            }
    //        cout<<endl;
    //    }
    //    for(int i=0; i<1296;i++)
    //    {
    //        cout<< Second_Umax_P[i]<<endl;
    //    }
    P_strategy.print("P: ");
    C_strategy.print("C: ");

//    for(int i = 0; i < colcount; i++) {
//        delete [] initial[i];
//    }
//    delete [] initial;
//
//    for(int i=0;i<rowcount;i++)
//    {
//        for(int j=0; j<colcount;j++)
//        {
//            delete[] first_handler[i][j];
//        }
//        delete [] first_handler[i];
//    }
    //delete first_handler;

}



  main.cpp
  Test Project

  Created by Zixiang  Pei on 5/5/19.
  Copyright © 2019 Zixiang  Pei. All rights reserved.





