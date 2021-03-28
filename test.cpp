#include<stdio.h> // printf
#include<iostream>
#include<vector> 
#include<string>
#include<fstream>

#include"external/args_parser.h" // parsing commandline args
#include"external/json.h" // parsing jsons

#include"libs/writhe.h" // wr::curve, wr::writhe, wr::resample
#include"libs/energy.h" // wr::energy, wr::nudge
#include"constants/ini.h" // getInitialCondition()

using json = nlohmann::json; // convenience
namespace wr {
    void to_json(json& j, const Point& p) {
        std::vector<double> tmp {p.x,p.y,p.z};
        j =  json(tmp);
    }
}

int main(int argc, char **argv)
{
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help"))
    {
        printf("Twisted loops simulator!\n\n");
        printf("Usage:\n");
        printf(" *Find minimum energy of given configuration:\n");
        printf(" main.exe -L LINK -r THICKNESS [-l LENGTH] [-omega OMEGA] [-pts POINTS] [-tmax TMAX] [-o OUTPUT] [-v]\n");
        printf("\n");
        printf("  Defaults: LENGTH=1, OMEGA=0.66, POINTS=14, TMAX=5000, OUTPUT=out.json\n");
        printf("  -L LINK: linking number\n");
        printf("  -r THICKESS: cross-sectional radius of the beam\n");
        printf("  -l LENGTH: total length of the centreline\n");
        printf("  -omega OMEGA: dimensionless ratio of bending to twisting (for full, circular beam 0.66(6))\n");
        printf("  -pts POINTS: number of nodal points of cubic spline approximating the curve\n");
        printf("  -tmax TMAX: number of iterations in monte carlo\n");
        printf("  -v: include in output energy and writhe in all iterations\n");
        return 0;
    }
    else
    {
        if(input.cmdOptionExists("-L") && input.cmdOptionExists("-r"))
        {
            double ParamLink = atof(input.getCmdOption("-L").c_str());
            double ParamRadius = atof(input.getCmdOption("-r").c_str());
            double ParamLength = 1.0;
            double ParamOmega = 0.66;
            int ParamPoints = 14;
            int ParamTMax = 5000;
            std::string ParamOutpath = "out.json";
            if(input.cmdOptionExists("-l"))ParamLength = atof(input.getCmdOption("-l").c_str());
            if(input.cmdOptionExists("-omega"))ParamOmega = atof(input.getCmdOption("-omega").c_str());
            if(input.cmdOptionExists("-pts"))ParamPoints = atoi(input.getCmdOption("-pts").c_str());
            if(input.cmdOptionExists("-tmax"))ParamTMax = atoi(input.getCmdOption("-tmax").c_str());
            if(input.cmdOptionExists("-o"))ParamOutpath = input.getCmdOption("-o");
        
            auto c = wr::getInitialCondition(ParamLink);
            auto GlobalBestC = c;
            auto LocalBestC = c;
            double GlobalBestScore = wr::energy(c,ParamLink,ParamRadius,ParamLength,ParamOmega,0);
            double LocalBestScore = wr::energy(c,ParamLink,ParamRadius,ParamLength,ParamOmega,0);
            double Score = 0;

            std::vector<double> WinnerWr;
            std::vector<double> WinnerScores;
            std::vector<wr::Curve> MovieFrames;
            WinnerWr.push_back(wr::writhe(c));
            WinnerScores.push_back(GlobalBestScore);

            for(int t=0;t<ParamTMax;t++)
            {
                wr::Curve ProposeC;
                for(int rep=0;rep<30;rep++)
                {
                    double step=0.05+0.4*(ParamTMax-t)/ParamTMax;
                    ProposeC = wr::nudge(c,ParamRadius*step); // change c a little
                    if(t%20 > 2)
                    {
                        Score = wr::energy(ProposeC,ParamLink,ParamRadius,ParamLength,ParamOmega,0);
                    }
                    else // Annealing
                    {
                        double repulsion=0.5*(ParamTMax-t)/ParamTMax;
                        LocalBestScore = wr::energy(LocalBestC,ParamLink,ParamRadius,ParamLength,ParamOmega,repulsion);
                        Score = wr::energy(ProposeC,ParamLink,ParamRadius,ParamLength,ParamOmega,repulsion);
                    }    

                    if(Score < LocalBestScore)
                    {
                        if(fabs(wr::writhe(resample(c,50)) - wr::writhe(resample(ProposeC,50)) ) > 0.5 ) // Prevent tunelling
                        {
                            continue;
                        }


                        LocalBestC = ProposeC;
                        LocalBestScore = Score;
                        c = ProposeC;
                    }
                }
                if(t%10==0)
                {
                    MovieFrames.push_back(resample(c,50));
                }

                printf("%d writhe:%lf ene:%4.3e global:%4.3e\n",t,wr::writhe(resample(c,50)),Score,GlobalBestScore);

                if(Score < GlobalBestScore)
                {
                    GlobalBestC = ProposeC;
                    GlobalBestScore = Score;
                    WinnerWr.push_back(wr::writhe(resample(GlobalBestC,50)));
                    WinnerScores.push_back(GlobalBestScore);
                    c = ProposeC;
                }
            }

            json j;
            j["link"] = ParamLink;
            j["radius"] = ParamRadius;
            j["length"] = ParamLength;
            j["omega"] = ParamOmega;
            j["points"] = ParamPoints;
            j["tmax"] = ParamTMax;

            j["writhes"] = WinnerWr;
            j["scores"]  = WinnerScores;
            j["winner"] = GlobalBestC;
            j["movie"] = MovieFrames;

            std::cout<<j["writhes"]<<std::endl;

            std::ofstream myfile;
            myfile.open (ParamOutpath);
            myfile << j;
            myfile.close();

        }
        else
        {
            printf("Missing arguments -L LINK -r THICKNESS\n");
            printf("Specify linking number and thickness\n");
            return -1;
        }
    }

   // std::vector<double> q;
   // q.push_back(0.5);
   // q.push_back(0.7);
   // json j;
  //  j["vec"] = c;
 //   std::cout << j;
}
