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


#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED
#include <vector>
#include <map>

#include <Individual.h>
#include "CSelectionYear.h"
#include "CSelEvol.h"

/// TOOL FUCNTIONS
long double FGaussianR(const long double mean, const long double var);

long double FMean(std::vector<long double>& v);

long double FCov(std::vector<long double>& v1, std::vector<long double>& v2);

long double FVar(std::vector<long double>& v);


int FTranslatePopsize(std::vector<unsigned int>& YearPopSizeList);

int FTranslateGeneral(std::vector<unsigned int>& JuvenileMales, std::vector<unsigned int>& JuvenileFemales);

/// CORE FUNCTIONS
std::map<unsigned int, Cindividual> FInitPop(std::map<char, std::vector<unsigned int> >& Alive,
                                             unsigned int const& CurrentYear,  CSelectionYear const& YearSel);

int FContinuePop(std::map<char, std::vector<unsigned int> >& Alive, std::map<unsigned int, Cindividual>& Population ,
                 unsigned int const& CurrentYear, CSelectionYear const& SelYear);

int FImmigration(std::map<char, std::vector<unsigned int> >& Alive, std::map<unsigned int, Cindividual>& Population ,
                 unsigned int const& CurrentYear);

Cindividual FAgeSexClasses(bool const& Sex, bool const& Immigrant, unsigned int const& count, unsigned int const& Cohort,
                            unsigned int const& CurrentYear, unsigned int const& age);

//std::map<char, std::vector<unsigned int> > FInitAlive(std::map<unsigned int, Cindividual> const& Population);
int FReproHard(std::map<unsigned int, Cindividual>& Population, std::map<char, std::vector<unsigned int> >& Alive,
                CSelectionYear const& SelYear, unsigned int const& CurrentYear);

int FReproSoft(std::map<unsigned int, Cindividual>& Population, std::map<char, std::vector<unsigned int> >& Alive,
               CSelectionYear const& SelYear);

int FInheritance(std::map<unsigned int, Cindividual>& Population, unsigned int& focaljuv);

long double FFitnessRepro(bool const& Sex, unsigned int const& Age, long double const& a, long double const& e,
                          CSelectionYear const& SelYear);

int FSurvival( std::map<char, std::vector<unsigned int> >& Alive, std::map<unsigned int,
              Cindividual>& Population, unsigned int const& CurrentYear, CSelectionYear const& SelYear);

int FSurvivalSexAgeSoft(bool const& Sex, std::vector<unsigned int>& AliveV,
                  std::map<unsigned int, Cindividual>& Population,
                     unsigned int const& CurrentYear, CSelectionYear const& SelYear);

int FAgeSurvivor(std::map<unsigned int, Cindividual>& Population, std::vector<unsigned int>& StillAlive, unsigned int & SurvivorNB, long double& SumFitness,
                  std::vector<long double>& CumulSumFitness, std::vector<long double>& ListFitness,
                   std::map<unsigned int, unsigned int>& ListAtRisk, std::vector<unsigned int> const& ListInd);

int FSurvivalSexAgeHard(std::vector<unsigned int>& AliveV, std::map<unsigned int, Cindividual>& Population,
                     unsigned int const& CurrentYear, CSelectionYear const& SelYear);

int FLookupRecruitment(std::map<unsigned int, Cindividual>& Population, Clifestage const& FocalLS, unsigned int const& FocalID,
                       bool const& phi);

long double FFitnessSurvival(bool const& Sex, unsigned int const& Age, long double const& a, long double const& e,
                          CSelectionYear const& SelYear);

Clifestage FInitLifeStage(unsigned int const& CurrentYear, unsigned int const& age);

std::map<unsigned int, CSelectionYear> FVarSelection(unsigned int const& AdultMales, unsigned int const& AdultFemales,
                                                     unsigned int const& Juvenile);

long double FPhiShift(long double p, long double Deviation);

int FMaturation(std::map<char, std::vector<unsigned int> >& Alive, std::map<unsigned int, Cindividual>& Population,
                unsigned int const& CurrentYear);

long double FMeanBV(std::map<char, std::vector<unsigned int> > const& Alive, std::map<unsigned int, Cindividual> const& Population);

int FRealSel(std::map<char, std::vector<unsigned int> > const& PreviousAlive,
                      std::map<unsigned int, Cindividual> const& Population, bool const& Zygote, CSelEvol& StatsSelEvolYear);

/// WRITING FUNCTIONS
int FClearFiles();

int FWriteInput();

int FWritePedigree(std::map<unsigned int, Cindividual> const& Population, unsigned int const& RUN);

int FWriteCaptures(std::map<unsigned int, Cindividual> const& Population, unsigned int const& RUN);

int FWriteEvolSel(std::vector<CSelEvol> const& StatsSelEvol, unsigned int const& RUN);

#endif // MAIN_H_INCLUDED

