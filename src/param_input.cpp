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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream> //stringstream

#include "param_input.h"


//variable globale du fichier
bool inputCheckBool=true;// si true, affiche tout ce qui est lu, ligne par ligne, dans le fichier de parametres.


using namespace std;


// utilitaire de comparaison de chaines insensible [`a] la casse
//attention sortie non intuitive! voir des usages prec['e]dents
int cmp_nocase(const string& s, const string& s2)
{
	string::const_iterator p = s.begin();
	string::const_iterator p2 = s2.begin();

	while(p != s.end() && p2 != s2.end())
        {
            if (toupper(*p) != toupper(*p2)) return((toupper(*p)<toupper(*p2)) ? -1 : 1);
            ++p; ++p2;
        }
	return((s2.size()==s.size()) ? 0 : (s2.size()<s.size()) ? -1 : 1);
}
// vire les blancs [`a] droite
void rtrim(string *s)
{
	while ((s->length()>0)  && (s->substr(s->length()-1,1)) == " ")
        {
			s->erase(s->length()-1,s->length());
        }
}

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

int cmp_nocase_no__(const string& cs, const string& cs2) {

    string s(cs); // makes local copy that one could modify. Aletrnatively, (string& s, string& s2) [no const] would modify its arguments
    string s2(cs2);
    replaceAll(s,"_","");
    replaceAll(s2,"_","");

	return(cmp_nocase(s,s2));
}

/***********************************************************************************************/


int evaluateBool(bool &boolean, string buf) { // safe assignment of value `buf' to 'boolean'
	stringstream strstr(buf);
    string locstring;
    strstr>>locstring;
    if(cmp_nocase(locstring,"")==0 || cmp_nocase(locstring,"T")==0 ||
	   cmp_nocase(locstring,"True")==0 || cmp_nocase(locstring,"Yes")==0 ||
	   cmp_nocase(locstring,"Y")==0)
		boolean=true;
    else if(cmp_nocase(locstring,"F")==0 || cmp_nocase(locstring,"False")==0 ||
			cmp_nocase(locstring,"No")==0 || cmp_nocase(locstring,"N")==0)
		boolean=false;
    else {
        cout<<"(!) Suspicious specification for a boolean: "<< buf <<endl;
		cout<<"(!) Only \"\", \"T\", \"True\", \"Yes\", \"Y\", \"F\", \"False\", \"No\", and \"N\" are allowed"<<endl;
		cout<<"I exit..."<< endl;
		cin.get();
		exit(-1);
    }
	return(0);
}

/*********************************************************/

int seeks_settings_file_name(const string cmdlinefilename,string& settingsfilename) {
// conceived to be executed only if the file has been created.
string buf,var;
string::size_type pos;
ifstream settings(cmdlinefilename.c_str(),ios::in);
    if(!settings.is_open()) {
    	cout << "Unable to open file "<<cmdlinefilename<< endl;
    	cerr << "Unable to open file "<<cmdlinefilename<< endl;
    	cerr << "Possible cause: several processes try to read the same file\n on a client/file server architecture.\n";
    	cerr << "Further execution would ignore the command line option.\n";
    	cerr << "Use the CmdLineFileName setting to avoid this.\n";
    	cerr << "I exit for safety";
    	cin.get();
        exit(-1);
    }
    else
    do {
    		getline(settings,buf);
    		if(buf.length()==0) break;
    		while((buf[0]==' ')||(buf[0]=='\t')) {buf.erase(0,1);}//vire les blancs initiaux
    		pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
    		var=buf.substr(0,pos).c_str();
    		if(cmp_nocase(var,"SettingsFile")==0) {
    			stringstream strstr(buf.substr(pos+1));
    			strstr>>settingsfilename;
                cout<<"Will take settings from file "<<settingsfilename<<endl;
    			settings.close();
    			return 0;
    		}
    } while(!settings.eof());
settings.close();
return 0;
}

/************************************************************/


int read_settings_file(const string filename) {
	string buf,var;
	//stringstring good mostly for string/numbers conversions
	string::size_type pos;
	//int tempo;
	char bidon;
	//double val;
	///RcodeStString.clear();
	int lineCounter=1;// to check whether cmdline is the first argument



	ifstream settings(filename.c_str(),ios::in);
	if(!settings.is_open()) {
		cout << "Unable to open file "<<filename<< endl;
		cerr << "Unable to open file "<<filename<< endl;
		cin.get();
		exit(-1);
	}
	else {
		cout<<"Reading settings file "<<filename<<endl;
		do {
			getline(settings,buf);
			if (inputCheckBool) {cout<<"Read line:\n"<<buf<<endl;}
			if(buf.length()==0) goto nextline;
			while((buf[0]==' ')||(buf[0]=='\t')) {buf.erase(0,1);}//vire les blancs initiaux
			pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
			//		pos=std::min(buf.find('='),buf.length());
			if ((buf[pos])=='=') while (buf[pos-1]==' ') {buf.erase(pos-1,1); pos--;}// vire les blancs avant le =
			var=buf.substr(0,pos).c_str();
			if (inputCheckBool) {cout<<"Parsed: "<<var<<endl;}
			if(var.length()==0) goto nextline;

			if(cmp_nocase(var,"InputCheck")==0) {//tordu
				inputCheckBool=true;
				goto nextline;
			}


			if(cmp_nocase(var,"cmdlinefilename")==0) {
				// meaningful on command line only (and processed at the beginning of main()), but written in <cmdline> text file where it should be ignored
				if (lineCounter>1) {
					cerr<<"(!)(!) Cmdlinefilename setting not first argument. Anything could happen! I exit.";
					cin.get();
					exit(-1);
				}
				goto nextline;
			}

            if(cmp_nocase_no__(var,"Pause")==0) {
    /* Pause determines only two correct contexts for cin.get():
           cerr<< error message + if(cinGetOnError) cin.get() + exit
    and
           cout<< some info + if(pauseGP) cin.get() + execution continues
    cinGetOnError is true at declaration in genepop.cpp but may then be set to false in latin.cpp
    */
                string locstring;
                stringstream strstr(buf.substr(pos+1));
                strstr>>locstring;
                if(cmp_nocase_no__(buf.substr(pos+1),"Final")==0) pauseGP=true;
                if(cmp_nocase_no__(locstring,"OnError")==0) cinGetOnError=true;
                if((cmp_nocase_no__(locstring,"Default")==0)) {cinGetOnError=false; pauseGP=false;}
                goto nextline;}

            if(cmp_nocase(var,"SamplingSeed")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>_ptSamplingSeed;
				goto nextline;}
            if(cmp_nocase(var,"RunNumber")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>RunNumber;
				goto nextline;}

            if(cmp_nocase(var,"InitialAdultMales")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>InitialAdultMales;
				goto nextline;}
			if(cmp_nocase(var,"InitialAdultFemales")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>InitialAdultFemales;
				goto nextline;}
            if(cmp_nocase(var,"JuvenileMales")==0 || cmp_nocase(var,"JuvenileMales")==0) {
                JuvenileMales.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                        while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        JuvenileMales.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"JuvenileFemales")==0 || cmp_nocase(var,"JuvenileFemales")==0) {
                JuvenileFemales.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                        while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        JuvenileFemales.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"PropMalesImm")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>PropMalesImm;
				goto nextline;}
            if(cmp_nocase(var,"VAz")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VAz;
				goto nextline;}
            if(cmp_nocase(var,"VEz")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VEz;
				goto nextline;}
            if(cmp_nocase(var,"DiffImmAz")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>DiffImmAz;
				goto nextline;}
            if(cmp_nocase(var,"DiffImmEz")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>DiffImmEz;
				goto nextline;}
            if(cmp_nocase(var,"MeanReproIntercept")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>MeanReproIntercept;
				goto nextline;}
            if(cmp_nocase(var,"VarReproLatentIntercept")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VarReproLatentIntercept;
				goto nextline;}
            if(cmp_nocase(var,"MeanReproSlopeA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>MeanReproSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"VarReproSlopeA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VarReproSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"QuadraticReproSlopeA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>QuadraticReproSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"MeanReproSlopeE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>MeanReproSlopeE;
				goto nextline;}
            if(cmp_nocase(var,"VarReproSlopeE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VarReproSlopeE;
				goto nextline;}
            if(cmp_nocase(var,"SoftViabilitySelection")==0) {
                evaluateBool(SoftViabilitySelection,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"SoftFertilitySelection")==0) {
                evaluateBool(SoftFertilitySelection,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"QuadraticReproSlopeE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>QuadraticReproSlopeE;
				goto nextline;}
            if(cmp_nocase(var,"OptRhoA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>OptRhoA;
				goto nextline;}
            if(cmp_nocase(var,"OptRhoE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>OptRhoE;
				goto nextline;}
            if(cmp_nocase(var,"SurvivalInterceptAdultMale")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>SurvivalInterceptAdultMale;
				goto nextline;}
            if(cmp_nocase(var,"SurvivalInterceptAdultFemale")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>SurvivalInterceptAdultFemale;
				goto nextline;}
            if(cmp_nocase(var,"SurvivalInterceptJuvenileMale")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>SurvivalInterceptJuvenileMale;
				goto nextline;}
            if(cmp_nocase(var,"SurvivalInterceptJuvenileFemale")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>SurvivalInterceptJuvenileFemale;
				goto nextline;}
            if(cmp_nocase(var,"VarSurvivalLatentIntercept")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VarSurvivalLatentIntercept;
				goto nextline;}
            if(cmp_nocase(var,"MeanSurvivalSlopeA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>MeanSurvivalSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"VarSurvivalSlopeA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VarSurvivalSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"QuadraticSurvivalSlopeA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>QuadraticSurvivalSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"MeanSurvivalSlopeE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>MeanSurvivalSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"VarSurvivalSlopeE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>VarSurvivalSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"QuadraticSurvivalSlopeE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>QuadraticSurvivalSlopeA;
				goto nextline;}
            if(cmp_nocase(var,"MaxAge")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>MaxAge;
				goto nextline;}
            if(cmp_nocase(var,"OptPhiA")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>OptPhiA;
				goto nextline;}
            if(cmp_nocase(var,"OptPhiE")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>OptPhiE;
				goto nextline;}

            if(cmp_nocase(var,"MonitoringDuration")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>MonitoringDuration;
				goto nextline;}

            if(cmp_nocase(var,"WriteInput")==0) {
                evaluateBool(WriteInput,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WritePedigree")==0) {
                evaluateBool(WritePedigree,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteCaptures")==0) {
                evaluateBool(WriteCaptures,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteEvolSel")==0) {
                evaluateBool(WriteEvolSel,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"Clonal")==0) {
                evaluateBool(clonal,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"MendelianSegregation")==0) {
                evaluateBool(MendelianSegregation,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"NoExtinction")==0) {
                evaluateBool(NoExtinction,buf.substr(pos+1));
                goto nextline;}

			if(cmp_nocase(var,"SettingsFile")==0) { //ignored when read from file ! cf seeks_settings_file_name()
				// meaningful on command line only but written in cmdline.txt which is read _after_ settingsfile
				goto nextline;
			}


			if ( ! (var[0]=='%' || var[0]=='#' || var[0]=='/' )) // exclude comments
				cout<<"(!) Unknown keyword \""<<var<<"\" in file "<<filename<<"."<<endl;
			// ici test on cinGetOnError cannot work since this setting may not yet have been set
		nextline: ;	// the pleasure of sin
			lineCounter++;
		} while(!settings.eof());
	}
	settings.close();
	return 0;
}

