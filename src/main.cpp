//    This is the source code for Volator version 0.001
//    Volator simulates life-histories and population dynamics in a population.
//
//	  Copyright (C) 2016-2018  Timothee Bonnet - timotheebonnetc@gmail.com
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

//pre-processor instructions
#include <iostream>
#include <sstream> //stringstream
#include <fstream> // ofstream
#include <ctime>//clock
#include <vector>
#include <map>
#include <random> // For Gaussian sampling

using namespace std;

#include "MersenneTwister.h"
#include "Main.h"
#include "param_input.h"
#include "Individual.h"
#include "CLifestage.h"
#include "CSelectionYear.h"
#include "CSelEvol.h"

//
//#include "eigenmvn.h"
#ifndef M_PI
#define M_PI REAL(3.1415926535897932384626433832795029)
#endif

///**
//  Take a pair of un-correlated variances.
//  Create a covariance matrix by correlating
//  them, sandwiching them in a rotation matrix.
//*/
//Eigen::Matrix2d genCovar(double v0,double v1,double theta)
//{
//  Eigen::Matrix2d rot = Eigen::Rotation2Dd(theta).matrix();
//  return rot*Eigen::DiagonalMatrix<double,2,2>(v0,v1)*rot.transpose();
//}

/// GLOBAL VARIABLES ///



// POPULATION PARAMETERS //
bool clonal=false;
bool MendelianSegregation=true;
bool NoExtinction=false;
unsigned int InitialAdultMales=15;
unsigned int InitialAdultFemales=15;
//unsigned int InitialJuvenileMales=30;
//unsigned int InitialJuvenileFemales=30;

vector<unsigned int> JuvenileMales;
vector<unsigned int> JuvenileFemales;

double PropMalesImm=0.5;

///Trait distribution
long double VAz=1.;
long double VEz=1.;
long double DiffImmAz=0.;
long double DiffImmEz=0.;

///Selection parameters
bool SoftViabilitySelection=true;

bool SoftFertilitySelection=true;
long double MeanReproIntercept=2.; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)
long double VarReproLatentIntercept=0.1; //THIS IS A LATENT POISSON VARIANCE

long double MeanReproSlopeA=0.; //Selection gradient on additive genetic variation, on the latent scale
long double VarReproSlopeA=0.; //Variance in Selection gradient on additive genetic variation, on the latent scale
long double QuadraticReproSlopeA=0.; //Quadratic Selection gradient on additive genetic variation, on the latent scale

long double MeanReproSlopeE=0.; //Selection gradient on environmental variation, on the latent scale
long double VarReproSlopeE=0.; //Variance in Selection gradient on environmental variation, on the latent scale
long double QuadraticReproSlopeE=0.; //Quadratic Selection gradient on environmental variation, on the latent scale
long double OptRhoA = 0.;
long double OptRhoE = 0.;


long double SurvivalInterceptAdultMale=0.5;
long double SurvivalInterceptAdultFemale=0.5;
long double SurvivalInterceptJuvenileMale=0.5;
long double SurvivalInterceptJuvenileFemale=0.5;
long double VarSurvivalLatentIntercept=0.;

long double MeanSurvivalSlopeA = 2.;
long double VarSurvivalSlopeA = 0.;
long double QuadraticSurvivalSlopeA = 0.;
long double MeanSurvivalSlopeE = 0.;
long double VarSurvivalSlopeE = 0.;
long double QuadraticSurvivalSlopeE = 0.;
long double MaxAge = 4;
long double OptPhiA = 0.;
long double OptPhiE = 0.;


// MONITORING PARAMETERS //
unsigned int MonitoringDuration=10;

// WRITING PARAMETERS //
bool WriteInput=false;
bool WritePedigree=false;
bool WriteCaptures=false;
bool WriteEvolSel=false;

bool pauseGP=false;
bool cinGetOnError=false;

// SIMULATION PARAMETERS //

MTRand alea;
unsigned long int _ptSamplingSeed=67144630;
unsigned int RunNumber=2;
std::default_random_engine generator(std::random_device{}()); // Random number for Gaussian generator

string cmdlinefilename="cmdlineArguments.txt";
string settingsfilename="InputVolator.txt";//fichiers d'entree avec valeurs des parametres

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{

/////////////////////////////// READ PARAMETERS /////////////////////////////////////
if (argc>1)
     {
        // to give inline the name of the file in which command line is written

        string buf(argv[1]);
        string::size_type pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
        string var=buf.substr(0,pos).c_str();
        if(cmp_nocase(var,"CmdlineFileName")==0) cmdlinefilename=buf.substr(pos+1);
        ofstream cmdline(cmdlinefilename.c_str(),ios::out);
        for (int it=1;it<=argc;it++) cmdline<<argv[it]<<endl;
        cmdline<<endl;
        cmdline.close();
        // seeks optional SettingsFile in cmdline
        seeks_settings_file_name(cmdlinefilename,settingsfilename);
     }
    read_settings_file(settingsfilename); //cf migraine.cpp... READS and SETS...
    if (argc>1) read_settings_file(cmdlinefilename);

    alea.seed(_ptSamplingSeed);

    FClearFiles(); //UNCOMMENT IF YOU HAVE MULTIPLE WRITING ROUNDS
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// NOW DOING THE THINGS ////////////////////////////////////
    unsigned int CurrentYear(0);
    for (unsigned int RUN(1); RUN<=RunNumber; RUN++)
    {
        CurrentYear=0;
        cout << "Run number "<<RUN<<" begins"<<endl;
        cout << "Hello world! This year is "<< CurrentYear << endl;

        FTranslateGeneral(JuvenileMales,JuvenileFemales);//Convert the number of individuals per year into something manageable


        map<char, vector<unsigned int> > Alive;
        map<char, vector<unsigned int> > PreviousAlive;
        map<char, vector<unsigned int> > PreviousRecruits;

        map<unsigned int, CSelectionYear> YearSelection = FVarSelection(InitialAdultMales, InitialAdultFemales,
                                                                        JuvenileMales[0]+JuvenileFemales[0]);

    //Hypothetical survival before!
        map<unsigned int, Cindividual> Population=FInitPop(Alive, CurrentYear, YearSelection[CurrentYear]);

        vector<CSelEvol> StatsSelEvol(MonitoringDuration+1);
        vector<long double> MeanBV(MonitoringDuration+1);
        vector<long double> MeanBVRecruits(MonitoringDuration+1);
        vector<long double> SelReal(MonitoringDuration);
        vector<long double> SelRealRecruits(MonitoringDuration);

        PreviousAlive['M']=Alive['M'];
        PreviousAlive['F']=Alive['F'];

        PreviousRecruits['M']=Alive['M'];
        PreviousRecruits['F']=Alive['F'];
        ///BV0
        StatsSelEvol[CurrentYear].LivingZygoteBreedingValues = FMeanBV(PreviousAlive, Population);//before repro, we exclude juv
        StatsSelEvol[CurrentYear].LivingRecruitsBreedingValues = StatsSelEvol[CurrentYear].LivingZygoteBreedingValues;//before repro, we exclude juv

        FRealSel(PreviousAlive, Population, true, StatsSelEvol[CurrentYear]);///No viability selection, only fertility in adults

        //StatsSelEvol[CurrentYear].SelectionRecruits = StatsSelEvol[CurrentYear].SelectionZygotes;//FRealSel(PreviousAlive, Population, false);///No viability selection, only fertility in adults
///Cannot compute selection on recruits before juvenile recruitment

     //map<char, vector<unsigned int> > Alive=FInitAlive(Population);
        /////////////////////////////// INITIALIZE POPULATION /////////////////////////////////
        if (MonitoringDuration>1)
        {
            for (unsigned int CurrentYear(1); CurrentYear<(MonitoringDuration); CurrentYear++)//Must be <MD-1 because another round after the loop
                {
                    cout << "Hello world! This year is "<< CurrentYear << endl;
                    StatsSelEvol[CurrentYear].LivingZygoteBreedingValues = FMeanBV(Alive, Population);// with juveniles and adults

//                    PreviousRecruits.clear();
//                    PreviousRecruits['M']=Alive['M'];
//                    PreviousRecruits['F']=Alive['F'];


                    PreviousAlive.clear();
                    PreviousAlive = Alive; // all the zygotes, adults and still to mature juveniles
                    YearSelection = FVarSelection(Alive['M'].size(), Alive['F'].size(),
                                                  JuvenileMales[CurrentYear]+JuvenileFemales[CurrentYear]);

                    FMaturation(Alive, Population, CurrentYear); //no more juv in Alive

                    FSurvival(Alive, Population, CurrentYear, YearSelection[CurrentYear]);

                    FRealSel(PreviousRecruits, Population, false, StatsSelEvol[CurrentYear-1]);
                    PreviousRecruits.clear();
                    PreviousRecruits['M']=Alive['M'];
                    PreviousRecruits['F']=Alive['F'];
                    StatsSelEvol[CurrentYear].LivingRecruitsBreedingValues = FMeanBV(PreviousRecruits, Population);//before repro, we exclude juv

    if(((Alive['M'].size()) * (Alive['F'].size()))  ==0)///
    {
        CurrentYear=MonitoringDuration; //ABORT THE LOOP
        cout << "GOING EXTINCT"<<flush<<endl;
    }else{
            FImmigration(Alive, Population, CurrentYear);
                    FContinuePop(Alive, Population, CurrentYear, YearSelection[CurrentYear]);
        }

                    FRealSel(PreviousAlive, Population, true, StatsSelEvol[CurrentYear]);

                }
//            ///Last year ?? maybe change the maturation??
//            MeanBV[CurrentYear] = FMeanBV(Alive, Population);
//            PreviousAlive.clear();
//            PreviousAlive = Alive;
//            YearSelection = FVarSelection();
//            unsigned int CurrentYear(MonitoringDuration-1);
//            FMaturation(Alive, Population, CurrentYear);
//
//            cout << "Hello world! This year is "<< CurrentYear << endl;
//            FSurvival(Alive, Population, CurrentYear, YearSelection[CurrentYear]);
//            RecruitsAlive.clear();
//            RecruitsAlive['M']=Alive['M'];
//            RecruitsAlive['F']=Alive['F'];
//            MeanBVRecruits[CurrentYear] = FMeanBV(RecruitsAlive, Population);//before repro, we exclude juv
//            SelRealRecruits[CurrentYear]=FRealSel(RecruitsAlive, Population, false); ///NOT CORRECT
//
//
//            FContinuePop(Alive, Population, CurrentYear, YearSelection[CurrentYear]);
//            SelReal[CurrentYear] = FRealSel(PreviousAlive, Population, true);

        }//end if (MonitoringDuration>1)

        /////////////////////////////// TIME HAPPENS //////////////////////////////////////////

        FWriteInput();

        FWritePedigree(Population, RUN);

        FWriteCaptures(Population, RUN);

        FWriteEvolSel(StatsSelEvol, RUN);

    }//end for(unsigned int RUN(0); RUN<RunNumber; RUN++}

        return 0;
}// end main()
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
int FTranslateGeneral(vector<unsigned int>& JuvenileMales, vector<unsigned int>& JuvenileFemales)
{
    FTranslatePopsize(JuvenileMales);
    FTranslatePopsize(JuvenileFemales);
    return 0;
}//end FTranslateGeneral

int FTranslatePopsize(vector<unsigned int>& YearPopSizeList)
 {
    if (YearPopSizeList.size()==0)//no input = constant offspring production N=30 per sex
        {
            for (unsigned int i(0);i<MonitoringDuration;i++)
                {
                    YearPopSizeList.push_back(30);
                }
        }
    if (YearPopSizeList.size()==1)//Constant offspring production
        {
            for (unsigned int i(1);i<MonitoringDuration;i++)
                {
                   YearPopSizeList.push_back(YearPopSizeList[0]);
                }
        }
    if (YearPopSizeList.size()==2)//the first year is different from the other ones (or there are 2 years only)
        {
            for (unsigned int i(2);i<MonitoringDuration;i++)
                {
                    YearPopSizeList.push_back(YearPopSizeList[1]);
                }
        }
    if ((YearPopSizeList.size()<MonitoringDuration && YearPopSizeList.size()>2) ||
        (YearPopSizeList.size()>MonitoringDuration))
        {
            cerr<<"ERROR in TranslatePopsize(): YearPopSizeList.size="<<YearPopSizeList.size()<<
                ". It should be equal to 0,1,2 or to MonitoringDuration ("<<MonitoringDuration<<")"<<endl;
            cerr<<"I exit"<<endl;
            if (cinGetOnError)
			cin.get();
            exit(-1);
        }
    //else YearPopSizeList.size()==MonitoringDuration and needs no transformation

    return 0;
 }// end int FTranslatePopsize()

long double FGaussianR(const long double mean, const long double var)// Wrapper for univariate gaussian generator
{
    const long double sd(sqrt(var));// standard deviation
    //std::normal_distribution<long double> distribution(mean,sd);
    //long double number = distribution(generator);
    long double u = alea();
    long double v = alea();
    long double r = sqrt( -2*log(u) );
    long double number = mean + sd*r * sin(2*M_PI*v);
    return number;
}// end FGaussianR

long double FMean(vector<long double>& v)
{
    long double Mean(0.);
    unsigned int n(v.size());
    for (unsigned int i(0); i < n; i++)
    {
        Mean += v[i];
    }

    Mean = Mean/n;
    return Mean;
}// end FMean()

long double FCov(vector<long double>& v1, vector<long double>& v2)
{
    if(v1.size() != v2.size())
    {
        cerr << "Trying to compute the covariance between two variables of different size";
        cerr << "v1.size()="<<v1.size()<<" v2.size()="<<v2.size();
        cerr << "\nI am sad and upset, I exit" << endl;
        if (cinGetOnError==true) cin.get();
        exit(-1);
    }
    unsigned int n(v1.size());
    long double m1(FMean(v1));
    long double m2(FMean(v2));
    long double cov(0.);
    for (unsigned int i(0); i < n; i++)
    {
        cov += ((v1[i]-m1)*(v2[i]-m2));
    }
    long double nn(n);
    cov = cov/(nn);
return cov;
}//end FCov()

long double FVar(vector<long double>& v)
{
    long double var(FCov(v,v));// covariance with itself
    return var;
}//end FVar()


map<unsigned int, Cindividual> FInitPop(map<char, vector<unsigned int> >& Alive,
                                        unsigned int const& CurrentYear, CSelectionYear const& SelYear) // SHOULD I PASS BY POINTER INSTEAD?
{
    map<unsigned int, Cindividual> Population;

    vector<unsigned int> AdultMales(InitialAdultMales); //CAN I INIT WITH A GIVEN SIZE?

    char AM('M');
    vector<unsigned int> AdultFemales(InitialAdultFemales);
    char AF('F');
    unsigned int JuvNumb(JuvenileMales[CurrentYear]+ JuvenileFemales[CurrentYear]);
    vector<unsigned int> Juveniles(JuvNumb);
    char J('J');

    unsigned int count(1);// ID must start at 1, NOT 0

    for(unsigned int i(0); i<InitialAdultMales;i++)// CREATES ADULT MALES
    {
        Population[count]= FAgeSexClasses(1, 0,count, -1, CurrentYear, 1);// You can use count safely because the population is just created
        AdultMales[i] = count;
        count++;

    }
    for(unsigned int i(0); i<(InitialAdultFemales);i++)//
    {
        Population[count]=FAgeSexClasses(0, 0, count, -1, CurrentYear, 1);// You can use count safely because the population is just created
        AdultFemales[i] = count;
        count++;
    }

    for(unsigned int i(0); i<JuvenileMales[CurrentYear];i++)//
    {
        Population[count]=FAgeSexClasses(1, 0, count, CurrentYear, CurrentYear, 0);// You can use count safely because the population is just created
        Juveniles[i] = count;
        count++;
    }
    for(unsigned int i(0); i<JuvenileFemales[CurrentYear];i++)//
    {
        Population[count]=FAgeSexClasses(0, 0, count, CurrentYear, CurrentYear, 0);// You can use count safely because the population is just created
        Juveniles[i+JuvenileMales[CurrentYear]] = count;
        count++;
    }

    Alive[AM] = AdultMales;
    Alive[AF] = AdultFemales;
    Alive[J] = Juveniles;

    FReproSoft(Population, Alive, SelYear);

return Population;
}// end FInitPop()

int FContinuePop(map<char, vector<unsigned int> >& Alive, map<unsigned int, Cindividual>& Population ,
                  unsigned int const& CurrentYear, CSelectionYear const& SelYear)
{

    if(SoftFertilitySelection){

        //to find the last individual Key in Population
        unsigned int count(Population.size()+1);// the first ind is "1", so the last one is "n", we want "n+1"

        unsigned int JuvNumb(JuvenileMales[CurrentYear]+JuvenileFemales[CurrentYear]);
        vector<unsigned int> Juveniles(JuvNumb);

        for(unsigned int i(0); i<JuvenileMales[CurrentYear];i++)//
        {
            Population[count]=FAgeSexClasses(1, 0, count, CurrentYear, CurrentYear, 0);
            Juveniles[i] = count;
            count++;
        }
        for(unsigned int i(0); i<JuvenileFemales[CurrentYear];i++)//
        {
            Population[count]=FAgeSexClasses(0, 0, count, CurrentYear, CurrentYear, 0);
            Juveniles[i+JuvenileMales[CurrentYear]] = count;
            count++;
        }
        Alive['J'] = Juveniles;
        FReproSoft(Population, Alive, SelYear);
    }else{ //hard selection
        FReproHard(Population, Alive, SelYear, CurrentYear);
    }
    return 0;
}// end FContinuePop()


int FImmigration(map<char, vector<unsigned int> >& Alive, map<unsigned int, Cindividual>& Population, unsigned int const& CurrentYear)
{
    //to find the last individual Key in Population
    unsigned int count(Population.size()+1);// the first ind is "1", so the last one is "n", we want "n+1"

    ///Fnumber immigrants
    long double lambda(0);
    unsigned int nbind(0);
    nbind = Alive['M'].size()+Alive['F'].size() + Alive['J'].size();

    lambda = exp(0.2395763  + nbind * 0.0611026 - nbind*nbind * 0.0004846); //empirical estimate
    long double L=exp(-(lambda));//generates Poisson distribution
    long unsigned int k=0;
    double p=1.;
    do{
    k++;
    p*=alea();
    }while(p>L);
    k--;
    if(k<0){k=0;}// should not be necessary, but does not hurt

    unsigned int NbImm(k);
    bool sex(0);
    double sexrand(0.);
    for (unsigned int i(0); i<NbImm; i++)
    {
        sexrand = alea();
        if(sexrand<PropMalesImm)
        {
            sex=1; //males
        }else{
        sex=0;
        }
        Population[count] = FAgeSexClasses(sex, 1, count, CurrentYear+1, CurrentYear, 1);
        Population[count].BreedingValueZ = FGaussianR(DiffImmAz, VAz);
        Population[count].EnvValueZ = FGaussianR(DiffImmEz, VEz);
        if(sex)//males
        {
            Alive['M'].push_back(count);

        }else{//females
            Alive['F'].push_back(count);
        }
        count++;
    }
return 0;
}//end FImmigration()

Cindividual FAgeSexClasses(bool const& Sex, bool const& Immigrant, unsigned int const& count,
                           unsigned int const& Cohort, unsigned int const& CurrentYear,
                           unsigned int const& age)
{

    Cindividual Ind;
    Ind.IndividualKey = count;
    Ind.Sex = Sex;
    Ind.Immigrant = Immigrant;
    Ind.Cohort = Cohort;
    Ind.Mother = 0;
    Ind.Father = 0;
    Ind.BreedingValueZ = FGaussianR(0., VAz);
    Ind.EnvValueZ = FGaussianR(0., VEz);

    Ind.LifeHistory.push_back(FInitLifeStage(CurrentYear, age));

return Ind;
}// end FAgeSexClasses()

Clifestage FInitLifeStage(unsigned int const& CurrentYear, unsigned int const& age)
{
    Clifestage LifeStageInit;
    LifeStageInit.age = age;
    LifeStageInit.year = CurrentYear;
    LifeStageInit.repro = 0;
    LifeStageInit.survival = 1; //VERY IMPORTANT, MUST BE ONE GIVEN OUR LIFE CYCLE§
return LifeStageInit;
}//end FInitLifeStage()
int FReproHard(map<unsigned int, Cindividual>& Population, map<char, vector<unsigned int> >& Alive,
                CSelectionYear const& SelYear, unsigned int const& CurrentYear)
{
    vector<Clifestage> ::iterator cit;
    ///HARD SELECTION ON FEMALES
    Cindividual focalind;
    unsigned int NBAdultFemales(Alive['F'].size());
    //FEMALE ABSOLUTE FITNESS
    vector<long double> reprofemale(NBAdultFemales);
    vector<long double> CumulSumF;
    long unsigned int NbOffspring(0);
    long double lambda(0.);

    for (unsigned int i(0);i<NBAdultFemales;i++)
    {
        focalind = Population[Alive['F'][i]];
        lambda = FFitnessRepro(focalind.Sex, focalind.LifeHistory[focalind.LifeHistory.size()-1].age,
                                      focalind.BreedingValueZ, focalind.EnvValueZ, SelYear); // cumulative vector of fitness // SHOULD BE THE RELATIVE FITNESS FUNCTION

        long double L=exp(-(lambda));//generates Poisson distribution
        long unsigned int k=0;
        double p=1.;
        do{
        k++;
        p*=alea();
        }while(p>L);
        k--;
        if(k<0){k=0;}// should not be necessary, but does not hurt
        reprofemale[i] = k;
        NbOffspring += k;
        CumulSumF.push_back( (long double) NbOffspring );
    }

    long unsigned int NBJuvFemale(0);
    long unsigned int NBJuvMale(0);

    NBJuvFemale=round(NbOffspring/2); ///ASSUMES EVEN SEX-RATIO
    NBJuvMale = NbOffspring - NBJuvFemale;

     //to find the last individual Key in Population
    unsigned int counte(Population.size()+1);// the first ind is "1", so the last one is "n", we want "n+1"
    vector<unsigned int> Juveniles(NbOffspring);

    for(unsigned int i(0); i<NBJuvMale;i++)//
    {
        Population[counte]=FAgeSexClasses(1, 0, counte, CurrentYear, CurrentYear, 0);
        Juveniles[i] = counte;
        counte++;
    }
    for(unsigned int i(0); i<NBJuvFemale;i++)//
    {
        Population[counte]=FAgeSexClasses(0, 0, counte, CurrentYear, CurrentYear, 0);
        Juveniles[i+NBJuvMale] = counte;
        counte++;
    }
    Alive['J'] = Juveniles;
    //WE FIND THE OFFSPRING MOTHER FROM THE MOTHER POINT OF VIEW
    unsigned int JuvCount(0);
    for (unsigned int i(0);i<NBAdultFemales;i++)
    {
        if(reprofemale[i]>0)
        {
            for (unsigned int j(1); j<=reprofemale[i]; j++)
            {
                Population[ Alive['J'][JuvCount] ].Mother = Alive['F'][i];//give mother
                cit = Population[ Alive['F'][i] ].LifeHistory.end();
                cit--;
                (*cit).repro++;// give offspring
                (*cit).offspring[Alive['J'][JuvCount]]=true;///default survival of the offspring as a recruit
                JuvCount++;
            }
        }
    }


    ///NOW SOFT SELECTION ON MALES
    unsigned int NBAdultMales(Alive['M'].size());
    //MALE RELATIVE FITNESS
    vector<long double> fitmale(NBAdultMales);
    long double SumM(0.);
    vector<long double> CumulSumM;
    for (unsigned int i(0);i<NBAdultMales;i++)
    {
        focalind = Population[Alive['M'][i]];
        fitmale[i] = FFitnessRepro(focalind.Sex, focalind.LifeHistory[focalind.LifeHistory.size()-1].age,
                                    focalind.BreedingValueZ, focalind.EnvValueZ, SelYear); // cumulative vector of fitness // SHOULD BE THE RELATIVE FITNESS FUNCTION
        SumM += fitmale[i];
        CumulSumM.push_back(SumM);
    }

    //FILIATION PART
    long double random(0.);
    long double h(0);//can still become zero thanks to h--
    long double w(0.);

    for (unsigned int i(0); i<NbOffspring; i++)//Loop over juveniles who need parents
    {
        // for Males (Fathers)
        random=alea()*SumM;//nombre aleatoire compris entre 0 et Sum maximale
        while (w<random)
        {
            w=CumulSumM[h];
            h++;
        }
        h--;
        if(h<0){h=0;}

        Population[ Alive['J'][i] ].Father = Alive['M'][h];// give parent
        cit=Population[ Alive['M'][h] ].LifeHistory.end();
        cit--;
        (*cit).repro++;// give offspring
        (*cit).offspring[Alive['J'][i]]=true;///default survival of the offspring as a recruit
        h=0;//can still become zero thanks to h--
        w=0.;// ready for next round

        FInheritance(Population, Alive['J'][i]); // transmit breeding values
    }

return 0;
}//end FReproHard()

int FReproSoft(map<unsigned int, Cindividual>& Population, map<char, vector<unsigned int> >& Alive,
                CSelectionYear const& SelYear)
{
    Cindividual focalind;
    vector<Clifestage> ::iterator cit;
    unsigned int NBAdultFemales(Alive['F'].size());
    unsigned int NBAdultMales(Alive['M'].size());
    unsigned int NBJuveniles(Alive['J'].size());
    //FEMALE RELATIVE FITNESS
    vector<long double> fitfemale(NBAdultFemales);
    long double SumF(0.);
    vector<long double> CumulSumF;
    for (unsigned int i(0);i<NBAdultFemales;i++)
    {
        focalind = Population[Alive['F'][i]];
        fitfemale[i] = FFitnessRepro(focalind.Sex, focalind.LifeHistory[focalind.LifeHistory.size()-1].age,
                                      focalind.BreedingValueZ, focalind.EnvValueZ, SelYear); // cumulative vector of fitness // SHOULD BE THE RELATIVE FITNESS FUNCTION
        SumF += fitfemale[i];
        CumulSumF.push_back(SumF);
    }

    //MALE RELATIVE FITNESS
    vector<long double> fitmale(NBAdultMales);
    long double SumM(0.);
    vector<long double> CumulSumM;
    for (unsigned int i(0);i<NBAdultMales;i++)
    {
        focalind = Population[Alive['M'][i]];
        fitmale[i] = FFitnessRepro(focalind.Sex, focalind.LifeHistory[focalind.LifeHistory.size()-1].age,
                                    focalind.BreedingValueZ, focalind.EnvValueZ, SelYear); // cumulative vector of fitness // SHOULD BE THE RELATIVE FITNESS FUNCTION
        SumM += fitmale[i];
        CumulSumM.push_back(SumM);
    }

    //FILIATION PART
    long double random(0.);
    long double h(0);//can still become zero thanks to h--
    long double w(0.);

    for (unsigned int i(0); i<NBJuveniles; i++)//Loop over juveniles who need parents
    {
    // for Females (Mothers)
    random=alea()*SumF;//nombre aleatoire compris entre 0 et Sum maximale
        while (w<random)
        {
            w=CumulSumF[h];
            h++;
        }
    h--;
    if(h<0){h=0;}

    Population[ Alive['J'][i] ].Mother = Alive['F'][h];// give parent
    cit=Population[ Alive['F'][h] ].LifeHistory.end();
    cit--;
    (*cit).repro++;// give offspring
    (*cit).offspring[Alive['J'][i]]=true;///default survival of the offspring as a recruit
    h=0;//can still become zero thanks to h--
    w=0.;

    // for Males (Fathers)
    random=alea()*SumM;//nombre aleatoire compris entre 0 et Sum maximale
        while (w<random)
        {
            w=CumulSumM[h];
            h++;
        }
    h--;
    if(h<0){h=0;}

    Population[ Alive['J'][i] ].Father = Alive['M'][h];// give parent
    cit=Population[ Alive['M'][h] ].LifeHistory.end();
    cit--;
    (*cit).repro++;// give offspring
    (*cit).offspring[Alive['J'][i]]=true;///default survival of the offspring as a recruit
    h=0;//can still become zero thanks to h--
    w=0.;// ready for next round

    FInheritance(Population, Alive['J'][i]); // transmit breeding values
    }

return 0;
}// end FReproFitness

int FInheritance(map<unsigned int, Cindividual>& Population, unsigned int& focaljuv)
{
    long double BV(0);
    long double midP((Population[Population[focaljuv].Mother].BreedingValueZ + //Mother's BV
                      Population[Population[focaljuv].Father].BreedingValueZ)/2);//Father's BV
    BV = midP;
    if(MendelianSegregation)
    {
        BV += FGaussianR(0., VAz/2); //Mendelian segregation
    }

    if (clonal)
        {
            BV=Population[Population[focaljuv].Mother].BreedingValueZ;
        }
    Population[focaljuv].BreedingValueZ=BV;

    return 0;
}

long double FFitnessRepro(bool const& Sex, unsigned int const& age, long double const& a, long double const& e,
                          CSelectionYear const& SelYear)
{
    long double fitness(0.);
    long double intercept(1.);
    if (intercept<0)
        {
        cerr << "The interecpt of the fitness function is "<< intercept << endl;
        cerr << "\nbut we want to take its logarithm..." << endl;
        cerr << "\ndo you see the problem?" << endl;
        if (cinGetOnError==true) cin.get();
        exit(-1);
        }
    long double repintercept(0.);
    if(Sex)//Male
    {
        repintercept=SelYear.ReproInterceptMale;
    }else{
        repintercept=SelYear.ReproInterceptFemale;
    }

    fitness = exp( log(repintercept) + SelYear.ReproSlopeA*(a-OptRhoA) + SelYear.ReproSlopeE*(e-OptRhoE) +
                  SelYear.QuadraticReproSlopeA*(a-OptRhoA)*(a-OptRhoA) + SelYear.QuadraticReproSlopeE*(e-OptRhoE)*(e-OptRhoE));
    return fitness;
};// end FFitnessRepro()

int FSurvival( map<char, vector<unsigned int> >& Alive, map<unsigned int, Cindividual>& Population,
              unsigned int const& CurrentYear, CSelectionYear const& SelYear)
{
    if(SoftViabilitySelection)
    {
        FSurvivalSexAgeSoft(0, Alive['F'], Population, CurrentYear, SelYear);
        FSurvivalSexAgeSoft(1, Alive['M'], Population, CurrentYear, SelYear);
    }else{
        FSurvivalSexAgeHard(Alive['M'], Population, CurrentYear, SelYear);
        FSurvivalSexAgeHard(Alive['F'], Population, CurrentYear, SelYear);
    }
    //FSurvivalSexAge(Alive['J'], Population, CurrentYear, SelYear); /// Juv Alive should always be empty at this point in the code
    if (Alive['J'].size()>0)
    {
        cerr << "Oh wait! Shouldn't this be empty?!"<<flush<<endl;
        exit(-1);
    }
return 0;
}// end FSurvival


int FSurvivalSexAgeSoft(bool const& Sex, vector<unsigned int>& AliveV, map<unsigned int, Cindividual>& Population,
                     unsigned int const& CurrentYear, CSelectionYear const& SelYear)
{
    vector<unsigned int> StillAlive;
    // First count how many survivors
    unsigned int AdSurvivor(0);
    unsigned int JuvSurvivor(0);
    long double GroupSizeJuv(0);
    long double GroupSizeAd(0);
    long double SurvivalNumber(0.5);
    if(Sex)//Males
    {
        GroupSizeJuv = (long double) JuvenileMales[CurrentYear-1];
        SurvivalNumber = GroupSizeJuv*SelYear.SurvivalInterceptJuvenileMale + 0.5;
        JuvSurvivor = (unsigned int) SurvivalNumber;

        GroupSizeAd = (long double) AliveV.size() - GroupSizeJuv;
        SurvivalNumber = GroupSizeAd*SelYear.SurvivalInterceptAdultMale + 0.5 ;
        AdSurvivor = (int) SurvivalNumber;
    }else{//Females
        GroupSizeJuv = (long double) JuvenileMales[CurrentYear-1];
        SurvivalNumber = GroupSizeJuv*SelYear.SurvivalInterceptJuvenileFemale + 0.5;
        JuvSurvivor = (unsigned int) SurvivalNumber;

        GroupSizeAd = (long double) AliveV.size() - GroupSizeJuv;
        SurvivalNumber = GroupSizeAd*SelYear.SurvivalInterceptAdultFemale + 0.5;
        AdSurvivor = (int) SurvivalNumber;
    }

    vector<unsigned int>::iterator it;
    vector<long double>::iterator ldit;
    vector<Clifestage>::iterator cit;
    long double fitness(0.);
    vector<unsigned int> ListJuv;
    map<unsigned int, unsigned int> ListAtRiskJuv; // Keeps a list of individuals who might die, and remove individuals progressively
    vector<long double> ListFitnessJuv;
    long double SumFitnessJuv;
    vector<long double> CumulSumFitnessJuv;
    vector<unsigned int> ListAd;
    map<unsigned int, unsigned int> ListAtRiskAd; // Keeps a list of individuals who might die, and remove individuals progressively
    vector<long double> ListFitnessAd;
    long double SumFitnessAd;
    vector<long double> CumulSumFitnessAd;
// FEMALES
// Now compute "absolute" fitnesses

    for(it=AliveV.begin();it!=AliveV.end();it++)//
    {
        cit = Population[(*it)].LifeHistory.end();
        cit--; //the last element
        fitness = FFitnessSurvival(Population[(*it)].Sex, (*cit).age,
                        Population[(*it)].BreedingValueZ, Population[(*it)].EnvValueZ, SelYear);
        if((*cit).age>1)// Old Adult
        {
            ListAd.push_back((*it));
            ListAtRiskAd[(*it)] = (*it);
            ListFitnessAd.push_back(fitness);
            SumFitnessAd+=fitness;
            CumulSumFitnessAd.push_back(SumFitnessAd);
        }else
        {
            if((*cit).age==1)//Juv newly matured
            {
                ListJuv.push_back((*it));
                ListAtRiskJuv[(*it)] = (*it);
                ListFitnessJuv.push_back(fitness);
                SumFitnessJuv+=fitness;
                CumulSumFitnessJuv.push_back(SumFitnessJuv);
            }else{
cout<<"CRAP"<<flush<<endl; exit(-1);
            }
        }
    }//end females

//Now select who survives

    FAgeSurvivor(Population, StillAlive, AdSurvivor, SumFitnessAd, CumulSumFitnessAd, ListFitnessAd, ListAtRiskAd, ListAd);
    FAgeSurvivor(Population, StillAlive, JuvSurvivor, SumFitnessJuv, CumulSumFitnessJuv, ListFitnessJuv, ListAtRiskJuv, ListJuv);

    AliveV.clear();
    AliveV = StillAlive;

return 0;
}//end FSexSurvivor()

int FAgeSurvivor(map<unsigned int, Cindividual>& Population, vector<unsigned int>& StillAlive, unsigned int & SurvivorNB, long double& SumFitness,
                  vector<long double>& CumulSumFitness, vector<long double>& ListFitness,
                   map<unsigned int, unsigned int>& ListAtRisk, vector<unsigned int> const& ListInd)
 {
    vector<Clifestage>::iterator cit;
    long double random(0.);
    long double h(0.);
    long double w(0.);
if(SurvivorNB>ListInd.size()){SurvivorNB=ListInd.size();}
//Adults
    for (unsigned int i(0); i<SurvivorNB; i++)
    {
        h=0.;
        w=0.;
        random=alea()*SumFitness;
        while(w<random)
        {
            w=CumulSumFitness[h];
            h++;
        }
        h--;
        if(h<0){h=0;}

        ListFitness[h] = 0;//cannot be selected again!
        cit = Population.find(ListInd[h])->second.LifeHistory.end();
        cit--;
        (*cit).survival=true;
        StillAlive.push_back(ListInd[h]);
        ListAtRisk.erase(ListInd[h]);/// Removed from the list of dead

        ///NOW LOOKUP PARENTS AND WRITE WHETHER THE OFFSPRING FROM LAST YEAR RECRUITED
        FLookupRecruitment(Population, (*cit), ListInd[w], true);

        // Ready for next round
        SumFitness=0;
        for (unsigned int j(0); j< CumulSumFitness.size(); j++)
        {
            SumFitness+=ListFitness[j];
            CumulSumFitness[j]= SumFitness;
        }

    }// end(for AdFemaleSurvivor)
    ///Now must say the non-surviving ones died
    map<unsigned int, unsigned int>::iterator uuit;
    for (uuit=ListAtRisk.begin(); uuit!=ListAtRisk.end(); uuit++)
    {
        cit = Population.find((*uuit).second)->second.LifeHistory.end();
        cit--;
        (*cit).survival=false;
        FLookupRecruitment(Population, (*cit), (*uuit).second, false);
    }
     return 0;
 }//end FAgeSurvivor()

int FLookupRecruitment(map<unsigned int, Cindividual>& Population, Clifestage const& FocalLS, unsigned int const& FocalID,
                       bool const& phi)
{
    vector<Clifestage>::iterator cit;
    if(FocalLS.age==1)///We look at recruitment only for newly matured adults
        {
            if(Population[FocalID].Mother!=0)
            {
                cit = Population[Population[FocalID].Mother].LifeHistory.end();
                cit--; cit--;///Two steps backward to find last year offspring
//cerr<<"LHS mother="<<Population[Population[FocalID].Mother].LifeHistory.size()<<endl;
                if((*cit).offspring.find(FocalID) == (*cit).offspring.end()  )///TEST FOR EXISTENCE MAYBE ADD AN EXCEPTION
                {
                    cerr<<"I did not find this offspring "<<FocalID<<" in the mother"<<endl;
                    exit(-1);
                }else{
                        (*cit).offspring.find(FocalID)->second=phi;
                   }
            }
            if(Population[FocalID].Father!=0)
            {
                cit = Population[Population[FocalID].Father].LifeHistory.end();
//cerr<<"LHS father="<<Population[Population[FocalID].Father].LifeHistory.size()<<flush<<endl;
                cit--; cit--;///Two steps backward to find last year offspring
                if((*cit).offspring.find(FocalID) == (*cit).offspring.end()  )///TEST FOR EXISTENCE MAYBE ADD AN EXCEPTION
                {
                    cerr<<"I did not find this offspring "<<FocalID<<" in the father"<<endl; exit(-1);
                }else{
                        (*cit).offspring.find(FocalID)->second=phi;
                    }
            }
        }

    return 0;
}

int FSurvivalSexAgeHard(vector<unsigned int>& AliveV, map<unsigned int, Cindividual>& Population,
                     unsigned int const& CurrentYear, CSelectionYear const& SelYear)
{

    vector<unsigned int>::iterator it;
    vector<Clifestage>::iterator cit;
    vector<unsigned int> StillAlive;
    long double random(0.);
    long double p(0.);//just to initialize
    bool phi(true);
    for(it=AliveV.begin();it!=AliveV.end();it++)//
    {
        cit = Population[(*it)].LifeHistory.end();
        cit--; //the last element
        p=FFitnessSurvival(Population[(*it)].Sex, (*cit).age,
                        Population[(*it)].BreedingValueZ, Population[(*it)].EnvValueZ, SelYear);//Survival proba
        random=alea();
        if (random>=p)
        {
            phi = false;
            Population[(*it)].Death=CurrentYear;
        }else{
            phi=true;
            StillAlive.push_back((*it));
        }
        (*cit).survival=phi;

        ///NOW LOOKUP PARENTS AND WRITE WHETHER THE OFFSPRING FROM LAST YEAR RECRUITED
        FLookupRecruitment( Population, (*cit), (*it), phi);
    }

    if(NoExtinction)
    {
        if(StillAlive.size()==0)//TO be sure at least one survives
        {
            it=AliveV.end();
            it--;
    /// if(AliveV.size()>0)/// NOT SURE IF NICE, MAYBE MEANS SOMETHING IS WRONG IF TRIGGERED
            {
                cerr <<" ej ="<<AliveV.size()<<endl;
                StillAlive.push_back((*it));
                cit = Population[(*it)].LifeHistory.end();
                cerr<< "pop size=" << Population.size() << "LH size ="<< Population[(*it)].LifeHistory.size();
                cit--;
                (*cit).survival=true;
            }
        }
    }


    AliveV.clear();
    AliveV = StillAlive;
    return 0;
}// end FSurvivalSexAge()

long double FFitnessSurvival(bool const& Sex, unsigned int const& Age, long double const& a, long double const& e,
                           CSelectionYear const& SelYear)
{
    long double p(0.);//0 if Age is beyond limit
    if(Age<MaxAge)
    {
        long double SurvivalIntercept(0.5);
        if(Age>1)//ADULTS /// THE AGE LIMIT IS ONE BECAUSE Juveniles do not "exist" any more, but are in the adult vectors!
        {
            if(Sex)//MALES
            {
                SurvivalIntercept = SelYear.SurvivalInterceptAdultMale;
            }else{
                SurvivalIntercept = SelYear.SurvivalInterceptAdultFemale;
            }//FEMALES
        }else{ // Juveniles
            if(Sex)//Males
            {
                SurvivalIntercept = SelYear.SurvivalInterceptJuvenileMale;
            }else{
                SurvivalIntercept = SelYear.SurvivalInterceptJuvenileFemale;
            }
        }
        long double lp(0);
        lp = log(SurvivalIntercept/(1-SurvivalIntercept)) +  SelYear.SurvivalSlopeA*(a-OptPhiA) +  SelYear.QuadraticSurvivalSlopeA*(a-OptPhiA)*(a-OptPhiA) +
          SelYear.SurvivalSlopeE*(e-OptPhiE) +  SelYear.QuadraticSurvivalSlopeE*(e-OptPhiE)*(e-OptPhiE);
        p = 1/(1+exp(-lp));
    }//else p = 0

return p;
}// end FFitnessRepro()

map<unsigned int, CSelectionYear> FVarSelection(unsigned int const& AdultMales, unsigned int const& AdultFemales,
                                                     unsigned int const& Juvenile)
{
    map<unsigned int, CSelectionYear> YearSelection;
    CSelectionYear ThisYear;

    ///CAST SECTION
    long double AMrep(0.);
    long double AFrep(0.);

    long double AMnum(0.);
    long double AFnum(0.);
    long double Jnum(0.);

    if(SoftFertilitySelection)// maybe does not matter, but relative fitness will be not too far from absolute
    {
        AMnum = (long double) AdultMales;
        AFnum = (long double) AdultFemales;
        Jnum = (long double) Juvenile;

        AMrep = Jnum/AMnum;
        AFrep = Jnum/AFnum;
    }else{ // reproductive expectation when hard selection
        AMrep = MeanReproIntercept;
        AFrep = MeanReproIntercept;
    }

    //END cast section

    for (unsigned int i(0); i<MonitoringDuration; i++)
    {
        ///REPRODUCTION //ATTENTION WE USE A INSTEAD OF E!!!!!
        long double ReproSlopeShift(FGaussianR(0., VarReproSlopeA));
        ThisYear.ReproIntercept = exp(log(MeanReproIntercept) + FGaussianR(0., VarReproLatentIntercept) );///NOT USEFULL YET
        ThisYear.ReproInterceptMale = exp(log(AMrep) + FGaussianR(0., VarReproLatentIntercept) );
        ThisYear.ReproInterceptFemale = exp(log(AFrep) + FGaussianR(0., VarReproLatentIntercept) );
        ThisYear.ReproSlopeA = MeanReproSlopeA + ReproSlopeShift;
        ThisYear.ReproSlopeE = MeanReproSlopeA + ReproSlopeShift;//ATTENTION WE USE A INSTEAD OF E!!!!!
        ThisYear.QuadraticReproSlopeA = QuadraticReproSlopeA;
        ThisYear.QuadraticReproSlopeE = QuadraticReproSlopeA;

        ///SURVIVAL //ATTENTION WE USE A INSTEAD OF E!!!!!

        long double PhiShift(FGaussianR(0., VarSurvivalLatentIntercept));//additive shift to apply to all sex-age classes

        ThisYear.SurvivalInterceptAdultMale = FPhiShift(SurvivalInterceptAdultMale,PhiShift);
        ThisYear.SurvivalInterceptAdultFemale = FPhiShift(SurvivalInterceptAdultFemale,PhiShift);
        ThisYear.SurvivalInterceptJuvenileMale = FPhiShift(SurvivalInterceptJuvenileMale,PhiShift);
        ThisYear.SurvivalInterceptJuvenileFemale = FPhiShift(SurvivalInterceptJuvenileFemale,PhiShift);

        long double SurvSlopeShift(FGaussianR(0., VarSurvivalSlopeA));
        ThisYear.SurvivalSlopeA = MeanSurvivalSlopeA + SurvSlopeShift;
        ThisYear.SurvivalSlopeE = MeanSurvivalSlopeA + SurvSlopeShift;
        ThisYear.QuadraticSurvivalSlopeA = QuadraticSurvivalSlopeA;
        ThisYear.QuadraticSurvivalSlopeE = QuadraticSurvivalSlopeA;
        YearSelection[i] = ThisYear;
    }

return YearSelection;
}// end FVarSelection

long double FPhiShift(long double p, long double Deviation)
{
    long double PhiShift =0.;
    PhiShift = 1/(1+exp(-(log(p/(1-p))+Deviation)));
    return PhiShift;
}// end FPhiShift()


int FMaturation(map<char, vector<unsigned int> >& Alive, map<unsigned int, Cindividual>& Population,
                 unsigned int const& CurrentYear)
{
    unsigned int age(0);
    vector<unsigned int>::iterator vit;
        ///Adult part
    //Creates a new life history element to be filled in during the next year

    for(vit=Alive['M'].begin();vit!=Alive['M'].end();vit++)//
    {
        age = CurrentYear-Population[(*vit)].Cohort;
        Clifestage newlifestage;
        newlifestage = FInitLifeStage(CurrentYear, age);
        Population[(*vit)].LifeHistory.push_back(newlifestage);//For Next year
    }

    for(vit=Alive['F'].begin();vit!=Alive['F'].end();vit++)//
    {
        age = CurrentYear-Population[(*vit)].Cohort;
        Clifestage newlifestage;
        newlifestage = FInitLifeStage(CurrentYear, age);
        Population[(*vit)].LifeHistory.push_back(newlifestage);//For Next year
    }

    ///Juvenile part
    for(vit=Alive['J'].begin();vit!=Alive['J'].end();vit++)//
    {
        age = CurrentYear-Population[(*vit)].Cohort;
        Clifestage newlifestage;
        newlifestage = FInitLifeStage(CurrentYear, age);
        Population[(*vit)].LifeHistory[0]=newlifestage;
        if (Population[(*vit)].Sex==0)//FEMALES
        {
            Alive['F'].push_back((*vit));
        }else{
            Alive['M'].push_back((*vit));//MALES
        }
    }
    Alive['J'].clear();

    return 0;
}// end FMaturation()

long double FMeanBV(map<char, vector<unsigned int> > const& Alive, map<unsigned int, Cindividual> const& Population)
{
    long double MeanBV(0.);
    map<char, vector<unsigned int> >::const_iterator mit;
    vector<unsigned int>::const_iterator vit;
    unsigned int nbalive(0);
    for (mit=Alive.begin(); mit!=Alive.end();mit++)
    {
        nbalive += (*mit).second.size();
        for (vit=(*mit).second.begin(); vit!=(*mit).second.end(); vit++)
        {
            MeanBV+= Population.find((*vit))->second.BreedingValueZ;
        }
    }
    MeanBV = MeanBV/nbalive;
    return MeanBV;
}//end FMeanBV

int FRealSel(map<char, vector<unsigned int> > const& PreviousAlive, map<unsigned int, Cindividual> const& Population,
                     bool const& Zygote, CSelEvol& StatsSelEvolYear)
{
    long double RealSel(0.);
    map<char, vector<unsigned int> >::const_iterator mit;
    vector<unsigned int>::const_iterator vit;

    unsigned int n(0);
    for(mit=PreviousAlive.begin(); mit!=PreviousAlive.end();mit++)
    {
            n+=(*mit).second.size();
    //(PreviousAlive.find('M')->second.size() + PreviousAlive.find('F')->second.size() + PreviousAlive.find('J')->second.size());
    }
    vector<long double> phenotype(n);
    vector<long double> fitness(n);
    Clifestage life;
    Clifestage pastlife;
    unsigned int counter(0);
    unsigned int survival(0);
    vector<Clifestage>::const_iterator lifehistory;
    map<unsigned int, bool>::iterator oit;
    for (mit=PreviousAlive.begin(); mit!=PreviousAlive.end();mit++)//loop on males, females, old juveniles
    {
        for(vit=(*mit).second.begin(); vit!=(*mit).second.end();vit++)//loop on individuals
        {
            phenotype[counter] = Population.find((*vit))->second.BreedingValueZ + Population.find((*vit))->second.EnvValueZ;
            lifehistory = Population.find((*vit))->second.LifeHistory.end();
            lifehistory--;
            life = (*lifehistory);
            if(life.survival){survival=1;}else{survival=0;}//to avoid cast problems
            if(Zygote) /// If we measure selection from zygote to living zygote
                {
                    fitness[counter] = life.repro + 2*survival;
                }else{ /// If we measure selection from recruit to living recruit
                    fitness[counter] = 2*survival;
                    if(life.age>1)//One year old don't have recruits
                    {
                        lifehistory--;
                        pastlife = (*lifehistory);

                        if(pastlife.offspring.size()>0) ///If there were offspring born this year
                        {
                            for (oit=pastlife.offspring.begin();oit!=pastlife.offspring.end();oit++)
                            {
                                if((*oit).second)///if juvenile recruited
                                {
                                    fitness[counter]++;
                                }
                            }
                        }
                    }// end if(life.age>1)
                } //End if/else Zygote
            counter++;
        }//end loop on individuals
    }//end loop on age/sex classes

    ///standardizes fitness
    vector<long double>::iterator fit;
    long double meanfit = FMean(fitness);
    for (fit=fitness.begin(); fit!=fitness.end();fit++)
    {
        (*fit) = (*fit)/meanfit;
    }

    ///now can get to selection differentials and gradients
    if(Zygote)
        {
            StatsSelEvolYear.SelectionZygotes=FCov(fitness, phenotype);
            StatsSelEvolYear.PhenoVarZygotes=FVar(phenotype);
            StatsSelEvolYear.BetaZygotes = StatsSelEvolYear.SelectionZygotes / StatsSelEvolYear.PhenoVarZygotes;
        }else{
            StatsSelEvolYear.SelectionRecruits=FCov(fitness, phenotype);
            StatsSelEvolYear.PhenoVarRecruits=FVar(phenotype);
            StatsSelEvolYear.BetaRecruits = StatsSelEvolYear.SelectionRecruits/StatsSelEvolYear.PhenoVarRecruits;
        }
return 0;
}//end FRealSel()


/// WRITING FUNCTIONS ///
int FClearFiles() // Clear all the files that should be written during the run. // USELESS IF ONLY ONE WRITING ROUND
{
    if (WritePedigree==true)
    {
        ofstream Pedigree("Pedigree.txt");
        Pedigree <<"RUN"<<"\t"<< "Id"<<"\t"<< "Dam"<<"\t" << "Sire" << "\t" << "Immigrant" <<endl;
        Pedigree.close();
    }
    if (WriteCaptures==true)
    {
        ofstream Captures("Captures.txt");
        Captures<< "RUN" <<"\t" << "Key"<<"\t" << "ID"<<"\t"<<"Sex"<<"\t"<<"BVZ"<<"\t"<<"EVZ"<<"\t"<<"Age"<<"\t"<<"Year"<<"\t"<<"Cohort"<<"\t"<<"Rho"<<"\t"<<"Phi" "\t" << "Immigrant"<<endl;
        Captures.close();
    }
    if (WriteEvolSel==true)
    {
        ofstream EvolSel("EvolSel.txt");
        EvolSel<< "RUN" <<"\t" << "Year"<<"\t"<<"DeltaBV"<<"\t"<<"Sel"<<"\t"<<"Var"<<"\t"<<"Beta"<<"\t"
        <<"DeltaBVRecruits"<<"\t"<<"SelRecruits"<<"\t"<<"VarRecruits"<<"\t"<<"BetaRecruits"<<endl;
        EvolSel.close();
    }
return 0;
}// end FClearFiles()

int FWriteInput()
{
    if (WriteInput==true)
        {
            //cout << "Kitty!" << endl;
        }
    return 0;
}// end FWriteInput()


int FWritePedigree(map<unsigned int, Cindividual> const& Population, unsigned int const& RUN)
{
    map<unsigned int, Cindividual>::const_iterator it;
    if (WritePedigree==true)
    {
        ofstream Pedigree("Pedigree.txt", ios::app);
        for (it=Population.begin();it!=Population.end();it++)//
            {
                Pedigree << RUN <<"\t" << (*it).second.IndividualKey<<"\t"<< (*it).second.Mother<<"\t" << (*it).second.Father<< "\t" << (*it).second.Immigrant<< endl;
            }
        Pedigree.close();
    }
return 0;
}// end FWritePedigree()

int FWriteCaptures(std::map<unsigned int, Cindividual> const& Population, unsigned int const& RUN)
{
    map<unsigned int, Cindividual>::const_iterator it;
    if (WriteCaptures==true)
    {
        ofstream Captures("Captures.txt", ios::app);
        for (it=Population.begin();it!=Population.end();it++)//
        {
            vector<Clifestage>::const_iterator vit;

            for (vit=(*it).second.LifeHistory.begin(); vit!=(*it).second.LifeHistory.end(); vit++)
                {
                    Captures << RUN <<"\t" << (*it).first<<"\t" << (*it).second.IndividualKey<<"\t"<<(*it).second.Sex
                    <<"\t"<<(*it).second.BreedingValueZ<<"\t"<<(*it).second.EnvValueZ
                    <<"\t"<<(*vit).age<<"\t"<<(*vit).year<<"\t"<<(*it).second.Cohort<<"\t"<<(*vit).repro<<"\t"<<(*vit).survival<< "\t" << (*it).second.Immigrant<<endl;

                }
        }
        Captures.close();
    }
return 0;
}// end FWriteCaptures

int FWriteEvolSel(vector<CSelEvol> const& StatsSelEvol, unsigned int const& RUN)
{
    if(WriteEvolSel==true)
    {
        ofstream EvolSel("EvolSel.txt", ios::app);

        for(unsigned int i(0); i<(StatsSelEvol.size()-2); i++)//minus 2 because we don't care about the last year where we cannot define evolution
        {
            EvolSel << RUN <<"\t"<< i << "\t"
                <<(StatsSelEvol[i+1].LivingZygoteBreedingValues-StatsSelEvol[i].LivingZygoteBreedingValues) << "\t"
                << StatsSelEvol[i].SelectionZygotes << "\t" << StatsSelEvol[i].PhenoVarZygotes<<"\t" << StatsSelEvol[i].BetaZygotes<<"\t"
                << (StatsSelEvol[i+1].LivingRecruitsBreedingValues-StatsSelEvol[i].LivingRecruitsBreedingValues) << "\t"
                << StatsSelEvol[i].SelectionRecruits<< "\t"<< StatsSelEvol[i].PhenoVarRecruits<< "\t" <<StatsSelEvol[i].BetaRecruits << endl;
        }
        EvolSel.close();
    }
return 0;
}// end FWriteEvolSel()
