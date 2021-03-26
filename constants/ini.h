#include"../libs/writhe.h"
#include"math.h"

namespace wr
{

wr::Curve getInitialCondition(double link)
{
    wr::Curve c;
    double step = 2*M_PI/14.9;
    for(double s=0;s<2*M_PI;s+=step)
    {
        wr::Point tmp = wr::Point(0.1*1.5*sin(s),0.1*0.75*sin(2*s),0.1*0.05*cos(s));
        //wr::Point tmp = wr::Point(sin(s),1.5*cos(s),0.2*cos(2*s));
        c.push_back( tmp );
    }
    return c;



    //auto c = wr::Curve(14);
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
    return c;
}

}
