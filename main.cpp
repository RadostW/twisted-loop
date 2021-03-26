#include <iostream>
#include <vector>
#include <cstdio>
#include <stdlib.h> 

#include "simpsons.h"
#include "writhe.h"
#include "splines.h"
#include "point.h"
#include "randompts.h"

// Global to retain spline fit.
sp::spline gammaX, gammaY, gammaZ;
wr::Curve wrSplineCurve;
const int wrpoints = 50; // Number of subsampling points for writhe calc.

wr::Curve nearSplineCurve;
const int nearSamplingN = 500; // Number of subsampling points for steric interactions

int splinePoints; // Stores number of points of initCurve

double dl(double s);

void curveInit(const wr::Curve& initCurve)
{
    int n = initCurve.size();
    splinePoints = n;
    std::vector<double> xi(n),yi(n),zi(n);
    for(int i=0;i<n;i++)
    {
        xi[i] = initCurve.at(i).x;
        yi[i] = initCurve.at(i).y;
        zi[i] = initCurve.at(i).z;
    }
    gammaX.setPoints(xi);
    gammaY.setPoints(yi);
    gammaZ.setPoints(zi);

    wrSplineCurve = wr::Curve(wrpoints);
    double wrdx = 1.0/wrpoints;
    for(int i=0; i < wrpoints; i++)
    {
        wrSplineCurve.at(i) = wr::Point(gammaX(wrdx*i),gammaY(wrdx*i),gammaZ(wrdx*i));
    }

    nearSplineCurve = wr::Curve(nearSamplingN);
    double nedx = 1.0/nearSamplingN;
    for(int i=0; i < nearSamplingN; i++)
    {
        nearSplineCurve.at(i) = wr::Point(gammaX(nedx*i),gammaY(nedx*i),gammaZ(nedx*i));
    }

    return;

}

// Remember to run curveInit for each curve.
double curvature(double s)
{
    wr::Point gammaPrime(gammaX.prime(s),gammaY.prime(s),gammaZ.prime(s));
    wr::Point gammaBis(gammaX.bis(s),gammaY.bis(s),gammaZ.bis(s));
    return (gammaPrime^gammaBis).Norm()/pow(gammaPrime.Norm(),3);
}

// Remember to run curveInit for each curve.
double writhe()
{
    wr::writhe(wrSplineCurve);
}

double nearNeighbours(double near)
{
    std::vector<double> arccoords(nearSamplingN);
    for(int i=1; i < nearSamplingN; i++)
    {
        arccoords[i] = arccoords[i-1] + (nearSplineCurve[i] - nearSplineCurve[i-1]).Norm();
    }

    double accum = 0;
    for(int i=0; i<nearSamplingN; i++)
    {
        wr::Point r2 = nearSplineCurve[(nearSamplingN+i)%nearSamplingN];
        for(int j=0; j<nearSamplingN; j++)
        {
            wr::Point r5 = nearSplineCurve[(nearSamplingN+j)%nearSamplingN];

            if((r5-r2).Norm() < near && 
                    fmod(arccoords[i] - arccoords[j]+arccoords[nearSamplingN-1],arccoords[nearSamplingN-1]) > 0.501*M_PI * near &&
                    fmod(arccoords[j] - arccoords[i]+arccoords[nearSamplingN-1],arccoords[nearSamplingN-1]) > 0.501*M_PI * near
            )
            {
                accum += 0.5*fmod(2*arccoords[nearSamplingN-1] + arccoords[(nearSamplingN+i+1)%nearSamplingN] - arccoords[(nearSamplingN+i-1)%nearSamplingN],arccoords[nearSamplingN-1]);
                //printf("near");
                break;
            }
        }
    }
    return accum;

}

double squaredCurvature(double s)
{
    double tmp = curvature(s);
    return tmp*tmp;
}

double dl(double s)
{
    wr::Point gammaPrime(gammaX.prime(s),gammaY.prime(s),gammaZ.prime(s));
    return gammaPrime.Norm();
}

double squaredCurvatureIntegrand(double s)
{
    return squaredCurvature(s)*dl(s);
}

// Remember to run curveInit for each curve
double totalSquaredCurvature()
{
    double accum = 0;
    int n = splinePoints;
    double dx = 1.0/n;
    for(int i=0;i<n;i++)
    {
        accum += si::adaptiveSimpsons(squaredCurvatureIntegrand,i*dx,(i+1)*dx, 0.0001, 10);
    }
    return accum;
}

// Remember to run curveInit for each curve
double totalLength()
{
    return si::adaptiveSimpsons(dl,0,1, 0.0001, 10);
}

double inverseDistances()
{
    double accum = 0;
    int n = wrpoints;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)
    {
        accum += 1.0 / (0.1+(wrSplineCurve[i]-wrSplineCurve[j]).Norm());
    }
    return accum;
}

const double LkInit = 2.0;
double LkTarget = 1.5;
const double targetLength = 1.0;
double Lk = LkInit;
const double omega = 1.5;
double thickness = 0.17;
double repulsion = 0.;

double energy(const wr::Curve& curve)
{
    curveInit(curve);
    double length = totalLength();
    double squaredCurvature = totalSquaredCurvature();
    double Tw = Lk - writhe();
    double Omega = 2*M_PI*Tw / length;
    double contact = nearNeighbours(thickness);
    double inverseDist = inverseDistances();

    double lengthPenalty = 8000;
    double contactPenalty = 10000;
    return lengthPenalty*(length-targetLength)*(length-targetLength) + squaredCurvature + length*Omega*Omega*omega + contactPenalty*contact + repulsion * inverseDist;
}

void printCurve(const wr::Curve& c,FILE* out)
{
    for(int i=0;i<c.size();i++)
    {
        fprintf(out,"%lf %lf %lf ",c[i].x,c[i].y,c[i].z);
    }
    fprintf(out,"\n");
}

void printCurveJSON(const wr::Curve& c,double beadRadius,FILE* out)
{
    double scale = beadRadius / thickness;
    fprintf(out,"[\n");
    for(int i=0;i<c.size();i++)
    {
        fprintf(out,"[[%lf,%lf,%lf],%lf]",scale*c[i].x,scale*c[i].y,scale*c[i].z,beadRadius);
        if((i+1)<c.size())
        {
            fprintf(out,",\n");
        }
    }
    fprintf(out,"]\n");
}

wr::Curve resample(wr::Curve c, int pts)
{
    curveInit(c);
    auto ret = wr::Curve(pts);
    double dx = 1.0/pts;
    for(int i=0; i < pts; i++)
    {
        ret.at(i) = wr::Point(gammaX(dx*i),gammaY(dx*i),gammaZ(dx*i));
    }
    for(int i=0;i<pts;i++)
    {
        //printf("%lf %lf %lf\n",ret[i].x,ret[i].y,ret[i].z);
    }
    return ret;
}

wr::Curve symmetrized(wr::Curve c)
{
    wr::Curve phase1(c);
    int n = phase1.size();
    for(int i=0; i<n; i++)
    {
        wr::Point l = c[i];
        wr::Point r = c[(n-i)%n];
        phase1[i] = wr::Point( 0.5*(l.x-r.x), 0.5*(l.y-r.y), 0.5*(l.z+r.z) );
    }
    if(n%2!=0)
    {
        return phase1;
    }
    else
    {
        wr::Curve phase2(phase1);
        for(int i=0; i<n; i++)
        {
            wr::Point l = phase1[i];
            wr::Point r = phase1[(n-i+(n/2))%n];
            phase2[i] = wr::Point( 0.5*(l.x+r.x), 0.5*(l.y-r.y), 0.5*(l.z-r.z) );
        }
        return phase2;
    }
}


// =========================================================== main ============
// =========================================================== main ============
// =========================================================== main ============

int main(int argc,char *argv[])
{

    if(argc==3)
    {
        LkTarget = atof(argv[1]);
        thickness = atof(argv[2]);
        printf("Args given LkTarget=%lf thickness=%lf",LkTarget, thickness);
        //return 0;
    }
    else
    {
        printf("Not enough or too many args. (argc==%d)\n",argc);
        printf("Usage:\nprog LkTarget thickness\n");
        printf("Length is kept constant and equal to 1.\n");
        return 0;
    }

    wr::Curve c;
    double step = 2.0*M_PI/8.0;
    Lk=LkInit;


    //Initialize to a decent guess.
    for(double s=0;s<2*M_PI;s+=step)
    {
        wr::Point tmp = wr::Point(0.1*1.5*sin(s),0.1*0.75*sin(2*s),0.1*0.05*cos(s));
        //wr::Point tmp = wr::Point(sin(s),1.5*cos(s),0.2*cos(2*s));
        c.push_back( tmp );
    }

    /*
    c = wr::Curve(12);
    c[0]=wr::Point(0.102523,-0.216333,0.081163);
    c[1]=wr::Point(0.836539,0.288106,-0.214138);
    c[2]=wr::Point(1.646469,0.508622,-0.320211);
    c[3]=wr::Point(2.122765,0.198421,0.062402);
    c[4]=wr::Point(1.793967,-0.252365,0.438142);
    c[5]=wr::Point(0.927355,-0.412181,0.280933);
    c[6]=wr::Point(-0.142164,-0.259649,-0.033488);
    c[7]=wr::Point(-1.162099,-0.070698,0.308664);
    c[8]=wr::Point(-1.740711,0.318208,0.323718);
    c[9]=wr::Point(-1.681976,0.591317,-0.203766);
    c[10]=wr::Point(-1.108678,0.161709,-0.426858);
    c[11]=wr::Point(-0.473384,-0.287895,-0.178630);
    for(int i=0;i<12;i++)c[i]=c[i]*0.1;
    */
    c = wr::Curve(14);
    c[0]=wr::Point(0.       ,	0.        ,	0.008104);
    c[1]=wr::Point(0.060091 ,	0.003151  ,	-0.009004);
    c[2]=wr::Point(0.137288 ,	-0.023273 ,	0.020285);
    c[3]=wr::Point(0.20945  ,	-0.018543 ,	0.019411);
    c[4]=wr::Point(0.20945  ,	0.018543  ,	-0.019411);
    c[5]=wr::Point(0.137288 ,	0.023273  ,	-0.020285);
    c[6]=wr::Point(0.060091 ,	-0.003151 ,	0.009004);
    c[7]=wr::Point(0.       ,	0.	      ,-0.008104);
    c[8]=wr::Point(-0.060091,	0.003151  ,	0.009004);
    c[9]=wr::Point(-0.137288,	-0.023273 ,	-0.020285);
    c[10]=wr::Point(-0.20945 ,	-0.018543 ,	-0.019411);
    c[11]=wr::Point(-0.20945 ,	0.018543  ,	0.019411);
    c[12]=wr::Point(-0.137288,	0.023273  ,	0.020285);
    c[13]=wr::Point(-0.060091,	-0.003151 ,	-0.009004);
    //c = resample(c,14);
    c = symmetrized(c);

    double energyMin = energy(c);
    int n = c.size();

    double globalMin = 9999999999;
    wr::Curve globalMinCurve;

    FILE* out;
    out = fopen("out.txt","w");

    curveInit(c);
    printCurve(nearSplineCurve,out);
    printCurve(c,out);

    int tmax=8000;
    int tfin=4000;
    for(int t=0;t<tmax&&c.size()<15;t++)
    {

        if(t > 40) // Dont look during burn-in
        {
            if(energyMin < globalMin && repulsion==0.0 && Lk == LkTarget)
            {
                globalMin = energyMin;
                globalMinCurve = wr::Curve(resample(c,1.0*targetLength/thickness)); // i'th bead touches  (i+1)^th bead centre
            }
        }

        if(t==10) // End burn-in period
        {
            Lk=LkTarget;
            energyMin=energy(c);
        }
        if(t>10 && t%3000 == 0 && Lk<4) // Adding additional points
        {
            //c = resample(c,c.size()+2);
            //energyMin=energy(c);
            //printf("Poof: more pts %d\n",(int)c.size());
        }

        if(t<tfin && t>20 && t%20 < 2) // Annealing
        {
            repulsion=0.1*(tfin-t)/tfin;
            Lk=LkTarget;
            energyMin=energy(c);
        }
        else if(t<tfin && t>20 && t%80 < 3) // Annealing
        {
            repulsion=0.1*(tfin-t)/tfin;
            Lk=LkTarget;
            energyMin=energy(c);
        }
        else 
        {
            Lk=LkTarget;
            repulsion=0.;
        }


        wr::Curve tmp(c);
        wr::Curve beforeBigStep(c);
        for(int reps=0;reps<30;reps++)
        {
            for(int i=0;i<n;i++)
            {
                double tmpEnergy;

                tmp[i] = tmp[i]+wr::GetSphericalNormal()*0.007;
                //tmp[i] = tmp[i]+wr::Point(0,0,1)*0.01;
                //tmpEnergy = energy(tmp);
                //if( tmpEnergy < energyMin)
                //{
                //    energyMin = tmpEnergy;
                //    c = tmp; // Copy contents
                //}
                //else
                //{
                //    tmp = c; // Copy contents
                //}
             
            }
            tmp = symmetrized(tmp);
            double tmpEnergy = energy(tmp);
            if( tmpEnergy < energyMin)
            {
                energyMin = tmpEnergy;
                c = tmp; // Copy contents
            }
            else
            {
                tmp = c; // Copy contents
            }
        }
        wr::Curve afterBigStep(c);

        wr::Curve linearExtrapolation(c);
        for(int i=0;i<n;i++)
        {
            linearExtrapolation[i] = linearExtrapolation[i] + (afterBigStep[i]-beforeBigStep[i])*0.05;
        }
        double linearExtrapolationEnergy = energy(linearExtrapolation);
        while(linearExtrapolationEnergy < energyMin)
        {
            //printf(":)");
            c = linearExtrapolation;
            energyMin = linearExtrapolationEnergy;
            for(int i=0;i<n;i++) //Extrapolate one step more
            {
                linearExtrapolation[i] = linearExtrapolation[i] + (afterBigStep[i]-beforeBigStep[i])*0.05;
            }
            linearExtrapolationEnergy = energy(linearExtrapolation);
        }

        curveInit(c);
        printCurve(nearSplineCurve,out);
        printCurve(c,out);

        printf("%d curv:%lf len:%lf writhe:%lf contact: %lf invdist %lf; ene:%4.3e global:%4.3e\n",t,totalSquaredCurvature(),totalLength(),writhe(),nearNeighbours(thickness),repulsion*inverseDistances(),energyMin,globalMin);
        if(totalSquaredCurvature() < 0.01)
        {
            break;
        }
    }

    fclose(out);     

    char buffer [200];
    sprintf (buffer, "minimalEnergy_Lk=%s_thickness=%s.beam.json", argv[1], argv[2]);    
    FILE* globalMinFile = fopen(buffer,"w");
    double beadSize = 10.0; // Angstroms
    printCurveJSON(globalMinCurve,beadSize,globalMinFile);
    fclose(globalMinFile);

}
