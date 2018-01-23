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

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include "CLifestage.h"

class Cindividual
{
    public:
        Cindividual();
        virtual ~Cindividual();

        unsigned int IndividualKey;
        bool Sex;
        bool Immigrant;
        int Cohort;
        unsigned int Death;
        long double BreedingValueZ;
        long double EnvValueZ;
        unsigned int Mother;
        unsigned int Father;
        std::vector<Clifestage> LifeHistory;

    protected:

    private:
};

#endif // INDIVIDUAL_H
