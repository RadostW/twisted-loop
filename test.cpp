#include<stdio.h> // printf
#include<iostream>
#include<vector> 
#include<utility> //pair
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

bool CompareSecond(std::pair<wr::Curve,double> l,std::pair<wr::Curve,double> r)
{
    return l.second < r.second;
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
            int ParamIntersectPoints = 1000;
            int ParamTMax = 5000;
            std::string ParamOutpath = "out.json";
            if(input.cmdOptionExists("-l"))ParamLength = atof(input.getCmdOption("-l").c_str());
            if(input.cmdOptionExists("-omega"))ParamOmega = atof(input.getCmdOption("-omega").c_str());
            if(input.cmdOptionExists("-pts"))ParamPoints = atoi(input.getCmdOption("-pts").c_str());
            if(input.cmdOptionExists("-ipts"))ParamIntersectPoints = atoi(input.getCmdOption("-ipts").c_str());
            if(input.cmdOptionExists("-tmax"))ParamTMax = atoi(input.getCmdOption("-tmax").c_str());
            if(input.cmdOptionExists("-o"))ParamOutpath = input.getCmdOption("-o");
        
            auto InitialCurve = resample(wr::getInitialCondition(ParamLink),ParamPoints);
            double BestScore = wr::energy(InitialCurve,ParamLink,ParamRadius,ParamLength,ParamOmega,0,0);
            int SinceBestUpdate = 0;
            wr::nearSamplingN = ParamIntersectPoints;

            if(ParamPoints%2!=0)
            {
                printf("Odd number of points currently not supported.\n");
                printf("Set -pts POINTS to an even number\n");
            }

            int SimplexSize = 0;
            std::vector< std::pair < wr::Curve , double > > Simplex(SimplexSize);

            // Init locations
            if(ParamPoints%4==0)
            {
                SimplexSize = 1 + 3*(ParamPoints/4 - 1) + 2;
                Simplex = std::vector< std::pair < wr::Curve , double > >(SimplexSize);
                for(int i=0;i<SimplexSize;i++)
                {
                    Simplex[i] = std::make_pair(InitialCurve,BestScore);
                }

                for(int i=0;i< 3*(ParamPoints/4 - 1);i++)
                {
                    int tmp = i%3;
                    switch(tmp)
                    {
                        case 0:
                            Simplex[1+i].first[1 + i/3].x += ParamRadius;                        
                        break;
                        case 1:
                            Simplex[1+i].first[1 + i/3].y += ParamRadius;                        
                        break;
                        case 2:
                            Simplex[1+i].first[1 + i/3].z += ParamRadius;                        
                        break;
                    }
                    Simplex[1+i].first = nudge(Simplex[1+i].first,ParamRadius);
                }
                Simplex[SimplexSize-2].first[0].z += ParamRadius;
                Simplex[SimplexSize-1].first[ParamPoints/4].x += ParamRadius;
            }
            else
            {
                SimplexSize = 1 + 3*(ParamPoints/4) + 1;
                Simplex = std::vector< std::pair < wr::Curve , double > >(SimplexSize);
                for(int i=0;i<SimplexSize;i++)
                {
                    Simplex[i] = std::make_pair(InitialCurve,BestScore);
                }

                for(int i=0;i< 3*(ParamPoints/4);i++)
                {
                    int tmp = i%3;
                    switch(tmp)
                    {
                        case 0:
                            Simplex[1+i].first[1 + i/3].x += ParamRadius;                        
                        break;
                        case 1:
                            Simplex[1+i].first[1 + i/3].y += ParamRadius;                        
                        break;
                        case 2:
                            Simplex[1+i].first[1 + i/3].z += ParamRadius;                        
                        break;
                    }
                }
                Simplex[SimplexSize-1].first[0].z += ParamRadius;
            }

            // Set energies of verticies
            for(int i=0;i<SimplexSize;i++)
            {
                auto tmpc = symmetrized(Simplex[i].first);
                double tmpe = wr::energy(tmpc,ParamLink,ParamRadius,ParamLength,ParamOmega,0,0);
                Simplex[i] = std::make_pair(tmpc,tmpe);
            }

            if(input.cmdOptionExists("-v"))
            {
                json j2;
                j2["simplex"]=Simplex;
                std::ofstream myfile2;
                myfile2.open ("simplex_init.json");
                myfile2 << j2;
                myfile2.close();
            }

            int WrithePoints = 50;
            std::vector<double> WinnerWr;
            std::vector<double> WinnerScores;
            std::vector<wr::Curve> MovieFrames;
            WinnerWr.push_back(wr::writhe(resample(Simplex[0].first,WrithePoints)));
            WinnerScores.push_back(BestScore);

            double alpha=1.0;
            double gamma=2.0;
            double rho=0.5;
            double sigma=0.5;

            for(int t=0;t<ParamTMax*100;t++)
            {
                std::sort( Simplex.begin(), Simplex.end(), CompareSecond );

                if(Simplex[0].second < BestScore)
                {
                    BestScore = Simplex[0].second;
                    SinceBestUpdate = 0;
                }
                else
                {
                    SinceBestUpdate++;
                    if(SinceBestUpdate > 300)
                    {
                        break; //Stop optimizing
                    }
                }

                if(t%10==0)
                {
                    MovieFrames.push_back(resample(Simplex[0].first,100));   
                    double WritheTmp = wr::writhe(resample(Simplex[0].first,WrithePoints));
                    double ScoreTmp = Simplex[0].second;
                    WinnerWr.push_back(WritheTmp);
                    WinnerScores.push_back(ScoreTmp);
                    if(t%100==0)
                    {

                        printf("\n%d writhe:%lf ene:%6.6e sbu: %d\n",t/100,WritheTmp,ScoreTmp,SinceBestUpdate);
                    }
                }

                //Find centroid of all but worst point
                wr::Curve Centroid(ParamPoints);
                for(int i=0;i<SimplexSize-1;i++)
                {
                    Centroid = wr::add( Centroid , wr::scale( Simplex[i].first , 1.0 / (SimplexSize-1) ) );
                }

                // x_r = x_o + \alpha ( x_0 - x_{n+1} )
                wr::Curve Reflection;
                Reflection = wr::add ( Centroid , wr::scale( wr::add( Centroid , wr::scale( Simplex[SimplexSize-1].first , -1) ) , alpha ) );  
                double ScoreReflected = wr::energy(Reflection,ParamLink,ParamRadius,ParamLength,ParamOmega,0,0);

                // f(x_1) <= f(x_r) < f(x_n)
                if(Simplex[0].second <= ScoreReflected && ScoreReflected < Simplex[SimplexSize-2].second)
                {
                    Simplex[SimplexSize-1] = std::make_pair(Reflection,ScoreReflected);
                    if(input.cmdOptionExists("-v"))printf("r");
                    continue; // next iteration of Nelder-Mead
                }

                if(ScoreReflected < Simplex[0].second)
                {
                    //Expansion
                    wr::Curve Expansion;
                    Expansion = wr::add( Centroid , wr::scale( wr::add( Reflection, wr::scale( Centroid , -1)) , gamma) );
                    double ScoreExpanded = wr::energy(Expansion,ParamLink,ParamRadius,ParamLength,ParamOmega,0,0);

                    if(ScoreExpanded < ScoreReflected)
                    {
                        Simplex[SimplexSize-1] = std::make_pair(Expansion,ScoreExpanded);
                        if(input.cmdOptionExists("-v"))printf("e");
                        continue; // next iteration of Nelder-Mead
                    }
                    else
                    {
                        Simplex[SimplexSize-1] = std::make_pair(Reflection,ScoreReflected);
                        if(input.cmdOptionExists("-v"))printf("q");
                        continue; // next iteration of Nelder-Mead
                    }
                }
                else
                {
                    //Contraction
                    wr::Curve Contraction;
                    Contraction = wr::add( Centroid, wr::scale( wr::add( Simplex[SimplexSize-1].first , wr::scale( Centroid , -1) ) , rho));
                    double ScoreContracted = wr::energy(Contraction,ParamLink,ParamRadius,ParamLength,ParamOmega,0,0);

                    if(ScoreContracted < Simplex[SimplexSize-1].second)
                    {
                        Simplex[SimplexSize-1] = std::make_pair( Contraction, ScoreContracted);
                        if(input.cmdOptionExists("-v"))printf("c");
                        continue; // next iteration of Nelder-Mead
                    }
                    else
                    {
                        //Shrink
                        for(int i=1;i<SimplexSize;i++)
                        {
                            wr::Curve tmp;
                            tmp = wr::add( Simplex[0].first , wr::scale( wr::add(Simplex[i].first,wr::scale(Simplex[0].first,-1))  ,sigma));
                            double ShrinkScore = wr::energy(tmp,ParamLink,ParamRadius,ParamLength,ParamOmega,0,0);
                            Simplex[i] = std::make_pair(tmp,ShrinkScore);
                        }
                        printf("s");
                        continue;
                    }

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
            j["winner"] = Simplex[0].first;
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
}
