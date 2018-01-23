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

#include "CSelectionYear.h"

CSelectionYear::CSelectionYear()
{
    //ctor
    long double ReproIntercept=2.; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)
    long double ReproInterceptMale=2.; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)
    long double ReproInterceptFemale=2.; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)

    long double ReproSlopeA=0.; //Selection gradient on additive genetic variation, on the latent scale
    long double QuadraticReproSlopeA=0.; //Quadratic Selection gradient on additive genetic variation, on the latent scale

    long double ReproSlopeE=0.; //Selection gradient on environmental variation, on the latent scale
    long double QuadraticReproSlopeE=0.; //Quadratic Selection gradient on environmental variation, on the latent scale

    long double SurvivalInterceptAdultMale=0.5;
    long double SurvivalInterceptAdultFemale=0.5;
    long double SurvivalInterceptJuvenileMale=0.5;
    long double SurvivalInterceptJuvenileFemale=0.5;

    long double SurvivalSlopeA = 0;
    long double QuadraticSurvivalSlopeA = 0;
    long double SurvivalSlopeE = 0;
    long double QuadraticSurvivalSlopeE = 0;
}

CSelectionYear::~CSelectionYear()
{
    //dtor
}
