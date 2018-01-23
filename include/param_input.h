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

#ifndef PARAM_INPUT_H_INCLUDED
#define PARAM_INPUT_H_INCLUDED
#include <vector>

/// GLOBAL VARIABLES ///
// SIMULATION PARAMETERS //
extern unsigned int RunNumber;
extern unsigned long int _ptSamplingSeed;

// POPULATION PARAMETERS //
extern bool clonal;
extern bool MendelianSegregation;
extern bool NoExtinction;
extern unsigned int InitialAdultMales;
extern unsigned int InitialAdultFemales;
//extern unsigned int InitialJuvenileMales;
//extern unsigned int InitialJuvenileFemales;
extern std::vector<unsigned int> JuvenileMales;
extern std::vector<unsigned int> JuvenileFemales;

extern double PropMalesImm;

extern long double VAz;
extern long double VEz;
extern long double DiffImmAz;
extern long double DiffImmEz;

extern bool SoftViabilitySelection;
extern bool SoftFertilitySelection;

extern long double MeanReproIntercept; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)
extern long double VarReproLatentIntercept; //THIS IS A LATENT POISSON VARIANCE

extern long double MeanReproSlopeA; //Selection gradient on additive genetic variation, on the latent scale
extern long double VarReproSlopeA; //Variance in Selection gradient on additive genetic variation, on the latent scale
extern long double QuadraticReproSlopeA; //Quadratic Selection gradient on additive genetic variation, on the latent scale

extern long double MeanReproSlopeE; //Selection gradient on environmental variation, on the latent scale
extern long double VarReproSlopeE; //Variance in Selection gradient on environmental variation, on the latent scale
extern long double QuadraticReproSlopeE; //Quadratic Selection gradient on environmental variation, on the latent scale

extern long double OptRhoA;
extern long double OptRhoE;

extern long double SurvivalInterceptAdultMale;
extern long double SurvivalInterceptAdultFemale;
extern long double SurvivalInterceptJuvenileMale;
extern long double SurvivalInterceptJuvenileFemale;
extern long double VarSurvivalLatentIntercept;


extern long double MeanSurvivalSlopeA;
extern long double VarSurvivalSlopeA;
extern long double QuadraticSurvivalSlopeA;
extern long double MeanSurvivalSlopeE;
extern long double VarSurvivalSlopeE;
extern long double QuadraticSurvivalSlopeE;
extern long double MaxAge;
extern long double OptPhiA;
extern long double OptPhiE;

// MONITORING PARAMETERS //
extern unsigned int MonitoringDuration;

// WRITING PARAMETERS //
extern bool WriteInput;
extern bool WritePedigree;
extern bool WriteCaptures;
extern bool WriteEvolSel;

extern bool pauseGP;
extern bool cinGetOnError;
////////////////////////////////////////

int cmp_nocase(const std::string& s, const std::string& s2);
void rtrim(std::string *s);
int evaluateBool(bool &boolean, std::string buf);
int seeks_settings_file_name(const std::string cmdlinefilename,std::string& settingsfilename);
int read_settings_file(const std::string filename);

#endif // PARAM_INPUT_H_INCLUDED

