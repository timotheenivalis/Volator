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

#include "Individual.h"

Cindividual::Cindividual()
{
    //ctor
    unsigned int IndividualKey = 0;
    bool Immigrant = 0;
    bool Sex = 0; // female 0, male 1
    int Cohort = 0;
    unsigned int Death = 0;
    long double BreedingValueZ = 0.;
    long double EnvValueZ = 0.;
    unsigned int Mother = 0;
    unsigned int Father = 0;

}

Cindividual::~Cindividual()
{
    //dtor
}
