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

#ifndef CSELECTIONYEAR_H
#define CSELECTIONYEAR_H


class CSelectionYear
{
    public:
        CSelectionYear();
        virtual ~CSelectionYear();
    long double ReproIntercept; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)
    long double ReproInterceptMale; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)
    long double ReproInterceptFemale; //THIS IS AN INTERCEPT ON THE OBSERVED SCALE (we will take its log in the latent Poisson function)

    long double ReproSlopeA; //Selection gradient on additive genetic variation, on the latent scale
    long double QuadraticReproSlopeA; //Quadratic Selection gradient on additive genetic variation, on the latent scale

    long double ReproSlopeE; //Selection gradient on environmental variation, on the latent scale
    long double QuadraticReproSlopeE; //Quadratic Selection gradient on environmental variation, on the latent scale

    long double SurvivalInterceptAdultMale;
    long double SurvivalInterceptAdultFemale;
    long double SurvivalInterceptJuvenileMale;
    long double SurvivalInterceptJuvenileFemale;

    long double SurvivalSlopeA;
    long double QuadraticSurvivalSlopeA;
    long double SurvivalSlopeE;
    long double QuadraticSurvivalSlopeE;

    protected:

    private:
};

#endif // CSELECTIONYEAR_H
