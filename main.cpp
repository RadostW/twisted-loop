#include<stdio.h> // printf
#include<string> // string
#include<fstream> // std::ofstrem

#include"external/args_parser.h" // parsing commandline args
#include"external/json.h" // parsing jsons

#include"libs/points.h" // Point, Curve
#include"libs/splines.h" // Spline
#include"libs/energy.h" // Energy of loop

#include"constants/ini.h" // GetInitialCondition()

using json = nlohmann::json; // convenience
namespace pt
{
void to_json(json& j, const pt::Point& p)
{
    std::vector<double> tmp {p.x,p.y,p.z};
    j =  json(tmp);
}
}

double OptionElseDouble(InputParser input, const std::string &option,double def)
{
    return input.cmdOptionExists(option) ? atof(input.getCmdOption(option).c_str()) : def;
}
std::string OptionElseString(InputParser input, const std::string &option,const std::string &def)
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
        double ParamLk = OptionElseDouble(input, "-L", 3.5);
        double ParamDiameter = OptionElseDouble(input, "-d", 20.0 / (3.34*336.0));
        double ParamElaEleNum = OptionElseDouble(input, "-E", 50.0);
        double ParamDebeyeLength = OptionElseDouble(input, "-K", 5.0 / (3.34*336.0));
        int ParamBasePairs = OptionElseDouble(input, "-B", 336);
        double ParamOmega = OptionElseDouble(input, "-omega", 2.0/3.0);
        int ParamPoints = OptionElseDouble(input, "-pts", 14);
        int ParamTMax = OptionElseDouble(input, "-tmax", 5000);
        std::string ParamOutpath = OptionElseString(input, "-o", "out.json");
        double ParamStericPenalty = 100000;
        double ParamLengthPenalty = 500000;

        std::vector<pt::Curve> MovieFrames;

        pt::Curve c;
        // Start with sensible guess
        for(double s=-M_PI;s<M_PI;s+=0.01)
        {
            // pt::Point tmp = pt::Point(0.2*sin(s),-0.05*sin(2*s),-0.02*cos(s)); //Old initial guess
            // pt::Point tmp = pt::Point(
            //                    0.194*sin(s) - 0.022*sin(3*s) - 0.002*sin(5*s),
            //                    0.0146*sin(2*s) - 0.013*sin(4*s),
            //                    0.0186*cos(s) - 0.0268*cos(3*s));
            //pt::Point tmp = pt::Point(
            //                    0.2*sin(s) - 0.0177*sin(3*s) - 0.0008*sin(5*s),
            //                    0.025*sin(2*s) - 0.01873*sin(4*s) - 0.0028*sin(6*s) + 
            //                    0.0009*sin(8*s),-0.003486*cos(s) + 0.0066566*cos(3*s) - 
            //                    0.0057575*cos(5*s) - 0.0026587*cos(7*s));
            double q = -6; // has to be even int
            pt::Point tmp = pt::Point( 
                                0.2*sin(s) - 0.02*sin(3*s),
                                -ParamDiameter*sin(q*s),
                                -ParamDiameter*cos(q*s)*tanh(q*(-M_PI/2. + s))*tanh(q*(M_PI/2. + s))
                                );
            c.push_back( tmp );
        }
        //c = pt::GetInitialCondition(3.5);
        c = sp::Resample(c, ParamPoints);
        c = pt::Symmetrized(c);

        lp::Loop loop = lp::Loop(c);

        MovieFrames.push_back(loop.GetShape());

        double MinimalEnergy = lp::LoopEnergy(loop,ParamLk,ParamOmega,
                                          ParamElaEleNum,ParamBasePairs,ParamDebeyeLength,
                                          ParamStericPenalty,ParamDiameter,
                                          ParamLengthPenalty);

        for(int t=0;t<ParamTMax;t++)
        {
            auto tmp = loop.Copy();
            //double TimeScale = cbrt(1.0*(ParamTMax-t)/(1.0*ParamTMax));
            double TimeScale = 1.0;
            tmp.Nudge(0.05*TimeScale*ParamDiameter);
            auto tmpenergy = lp::LoopEnergy(tmp,ParamLk,ParamOmega,
                                          ParamElaEleNum,ParamBasePairs,ParamDebeyeLength,
                                          ParamStericPenalty,ParamDiameter,
                                          ParamLengthPenalty);
            if(t%2000==0)
            {
                printf("%6d: %2.3E %3.2lf %3.2lf %3.2lf\n",
                           t,
                           MinimalEnergy,
                           sp::Writhe(loop.sx,loop.sy,loop.sz),
                           sp::NearNeighbours(loop.sx,loop.sy,loop.sz,ParamDiameter),
                           sp::Length(loop.sx,loop.sy,loop.sz)
                        );
            }
            if(tmpenergy < MinimalEnergy)
            {
                printf("%6d: %2.3E %3.2lf %3.2lf %3.2lf\n",
                           t,
                           MinimalEnergy,
                           sp::Writhe(loop.sx,loop.sy,loop.sz),
                           sp::NearNeighbours(loop.sx,loop.sy,loop.sz,ParamDiameter),
                           sp::Length(loop.sx,loop.sy,loop.sz)
                        );
                MovieFrames.push_back(loop.GetShape());

                loop = tmp.Copy();
                MinimalEnergy = tmpenergy;

            }
            else
            {
                //printf("%d: %lf %lf\n",t,MinimalEnergy,tmpenergy);
            }
        }

        json j;
        j["Lk"] = ParamLk;
        j["Diameter"] = ParamDiameter;
        j["ElastoElectroNumber"] = ParamElaEleNum;
        j["DebyeLength"] = ParamDebeyeLength;
        j["BasePairs"] = ParamBasePairs;
        j["Omega"] = ParamOmega;
        j["TMax"] = ParamTMax;
        j["StericPenalty"] = ParamStericPenalty;
        j["LengthPenalty"] = ParamLengthPenalty;

        j["movie"] = MovieFrames;

        std::ofstream myfile;
        myfile.open (ParamOutpath);
        myfile << j;
        myfile.close();

    }

}
