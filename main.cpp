#include<stdio.h> // printf
#include"external/args_parser.h" // parsing commandline args
#include"external/json.h" // parsing jsons

#include"libs/points.h" // Point, Curve
#include"libs/splines.h" // Spline
#include"libs/energy.h" // Energy of loop

#include"constants/ini.h" // getInitialCondition()

using json = nlohmann::json; // convenience
void to_json(json& j, const Point& p)
{
    std::vector<double> tmp {p.x,p.y,p.z};
    j =  json(tmp);
}
double OptionElseDouble(InputParser input, const std::string &option,double def)
{
    return input.cmdOptionExists(option) ? atof(input.getCmdOption(option).c_str()) : def;
}
std::string OptionElseString(InputParser input, const std::string &option,const std::string &option def)
{
    return input.cmdOptionExists(option) ? input.getCmdOption(option) : def;
}

int main(int argc, char **argv)
{
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help"))
    {
        printf("Twisted loops simulator!\n\n");
        printf("Usage:\n");
        printf(" *Find minimum energy of given configuration:\n");
        printf(" main.exe [-L LINK] [-r THICKNESS] [-E ELAELENUM] [-K DEBLEN] [-bp BASEPAIRS] [-omega OMEGA] [-pts POINTS] [-tmax TMAX] [-o OUTPUT]\n");
        printf("\n");
        printf("  -L LINK: linking number\n");
        printf("  -d THICKESS: cross-sectional diameter of the beam [multiples of beam lengh]\n");
        printf("  -E ELAELENUM: elasto-electrostaic number [dimensionless]\n");
        printf("  -K DEBLEN: Debeye length for electrostatics [multiples of beam length]\n");
        printf("  -B BASEPARIS: number of basepairs\n");
        printf("  -omega OMEGA: dimensionless ratio of bending to twisting (for full, circular beam 0.66(6))\n");
        printf("  -pts POINTS: number of nodal points of cubic splines approximating the curve\n");
        printf("  -tmax TMAX: number of iterations in monte carlo\n");
        printf("  -o OUTPUT: path to write output to\n");
        
        return 0;
    }
    else
    {      
        double ParamLink = OptionElseDouble(input, "-L", 3.5);
        double ParamRadius = OptionElseDouble(input, "-d", 20.0 / (3.34*336.0));
        double ParamElaEleNum = OptionElseDouble(input, "-E", 300.0);
        double ParamDebeyeLength = OptionElseDouble(input, "-K", 30.0 / (3.34*336.0));
        int ParamBasePairs = OptionElseDouble(input, "-B", 336);
        double ParamOmega = OptionElseDouble(input, "-omega", 2.0/3.0);
        int ParamPoints = OptionElseDouble(input, "-pts", 14);
        int ParamTMax = OptionElseDouble(input, "-tmax", 5000);
        std::string ParamOutpath = OptionElseString(input, "-o", "out.json");
    }

}
