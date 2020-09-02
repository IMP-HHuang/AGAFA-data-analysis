
#include "tree.h"
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <TROOT.h>

using namespace std;

string getTime()
{
	time_t timep;
	time (&timep);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
	return tmp;	
}

int main(int argc, char* argv[])
{
		if(argc < 2)
		{
			cout << "eg. ./main [run number] or ./main [start run number] [end run number]" << endl;
			exit(1);
		}
		int run1 = 0;
		int run2 = 0;
		if(argc == 2)
		{
			run1 = atoi(argv[1]);
			run2 = atoi(argv[1]);
		}
		if(argc == 3) 
		{
				run1 = atoi(argv[1]);
				run2 = atoi(argv[2]);
		}

		for(int run = run1; run <= run2; run++)
		{
			auto ipf = new TFile(Form("../data/Bi/ana_3rd/dssd%d.root", run), "READ");
			ipf->cd();
			auto ipt = (TTree*)ipf->Get("opt");
			string time = getTime();
			cout<< "start " << TString::Format("run%d", run) << " at " << time << endl;
			auto opf = new TFile(Form("../data/Bi/correlation/correlation%d.root", run), "RECREATE");
			opf->cd();
			auto opt = new TTree("opt", "opt");
			auto it = new tree(ipt);
			it->Loop(opt);
			it->fill();
			opt->Write();
			opf->Close();
			ipf->Close();
			delete ipf;
			delete opf;
			delete it;
			time = getTime();
			cout<< "Finsihed " << TString::Format("run%d", run) << " at " << time << endl;
		}

		return 0;
}

