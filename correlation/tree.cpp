//tree.cpp


#include "tree.h"
#include <TTree.h>
#include <vector>
#include <iostream>

using namespace std;

void tree::Init()
{
	
	if(ipt==NULL) 
	{
		cout << "Cannot Get input tree!" << endl;
		return;
	}
	mesum = 0;
	xesum = 0;
	yesum = 0;
	desum = 0;
	xvec = 0;
	yvec = 0;
	boxvec = 0;
	pavec = 0;
	mwvec = 0;
	xavec = 0;
	gsvec = 0;
	ipt->SetBranchAddress("mesum", &mesum);
	ipt->SetBranchAddress("xesum", &xesum);
	ipt->SetBranchAddress("yesum", &yesum);
	ipt->SetBranchAddress("desum", &desum);
	ipt->SetBranchAddress("x", &xvec);
	ipt->SetBranchAddress("y", &yvec);
	ipt->SetBranchAddress("box", &boxvec);
	ipt->SetBranchAddress("pa", &pavec);
//	ipt->SetBranchAddress("mw", &mwvec);
	ipt->SetBranchAddress("xa", &xavec);
	ipt->SetBranchAddress("gs", &gsvec);
}

void tree::Loop(TTree *opt_)
{	
	if (opt_==NULL) return;
	opt = opt_;
	BranchOpt();
	Long64_t nentries = ipt->GetEntriesFast();
	for(Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		ipt->GetEntry(jentry);
		if((*xvec).size()==0 || (*yvec).size()==0)
			continue;
		Double_t pasum = 0;
		for(auto pa : (*pavec))
			pasum = pasum + pa.e;
		//Recoil
		if((mesum > 0 || pasum > 0) && desum > 10000 && desum < 50000)
		{
//			cout << "Reciol start" << endl;
			event.x = -1.0;
			event.y = -1.0;
			event.e = -1.0;
			event.emin = -1.0;
			event.ts = -1;
			chain.clear();
			event.emin = desum;
			if((*xvec).size()==1)
			{
				event.x = (*xvec)[0].id;
				event.e = (*xvec)[0].e;
				event.ts = (*xvec)[0].ts;
				if((*yvec).size()==1)
				{
					event.y = (*yvec)[0].id;
				}
				else
				{
					event.y = Float_t((*yvec)[0].id + (*yvec)[1].id)/2.0;
				}
			}
			else
			{
				event.x = Float_t((*xvec)[0].id + (*xvec)[1].id)/2.0;
				event.ts = (*xvec)[0].ts;
				if((*yvec).size()==1)
				{
					event.e = (*yvec)[0].e;
					event.y = (*yvec)[0].id;
					event.ts = (*yvec)[0].ts;
				}
				else
				{
					event.e = xesum;
					event.y = Float_t((*yvec)[0].id + (*yvec)[1].id)/2.0;
				}
			}
			chain.push_back(event);
//			cout << "(x,y) is" << int(event.x*10.0) << "  " <<int(event.y*10.0) << endl;
			multi_chain[int(event.x*10.0)][int(event.y*10.0)].push_back(chain);
//			cout << "Recoil over" <<endl;
		}
		//Decay
		if(mesum==0 && pasum==0 && desum < 10000)
		{
//			cout << "Decay start" << endl;
			event.x = -1.0;
			event.y = -1.0;
			event.e = -1.0;
			event.emin = desum;
			event.ts = -1;
			if((*xvec).size()==1 && (*yvec).size()==1)
			{
				event.x = (*xvec)[0].id;
				event.y = (*yvec)[0].id;
				event.e = (*xvec)[0].e;
				event.ts = (*xvec)[0].ts;
			}
			else if((*xvec).size()!=1 && (*yvec).size()==1)
			{
				event.x = Float_t((*xvec)[0].id + (*xvec)[1].id)/2.0;
				event.y = (*yvec)[0].id;
				event.e = (*yvec)[0].e;
				event.ts = (*yvec)[0].ts;
			}
			else if((*xvec).size()==1 && (*yvec).size()!=1)
			{
				event.x = (*xvec)[0].id;
				event.y = Float_t((*yvec)[0].id + (*yvec)[1].id)/2.0;
				event.e = (*xvec)[0].e;
				event.ts = (*xvec)[0].ts;
			}
			else if((*xvec).size()!=1 && (*yvec).size()!=1)
			{
				event.x = Float_t((*xvec)[0].id + (*xvec)[1].id)/2.0;
				event.y = Float_t((*yvec)[0].id + (*yvec)[1].id)/2.0;
				event.e = xesum;
				event.ts = (*xvec)[0].ts;
			}
			alpha.x = event.x;
			alpha.y = event.y;
			alpha.e = event.e;
			alpha.emin = event.emin;
			alpha_all[int(event.x*10.0)][int(event.y*10)].insert(pair<ULong64_t, Alpha>(event.ts, alpha));
			int i = 0;
			Event eevent;
			for(auto &sub_chain : multi_chain[int(event.x*10.0)][int(event.y*10)])
			{
				if(sub_chain.size() >= chain_long)
					continue;
				else
				sub_chain.push_back(event);
				i++;
			}
			if(int(event.x*10.0) > 0)
			{
				alpha_all[int(event.x*10.0)-5][int(event.y*10)].insert(pair<ULong64_t, Alpha>(event.ts, alpha));
				i = 0;
				for(auto &sub_chain : multi_chain[int(event.x*10.0)-5][int(event.y*10)])
				{
					if(sub_chain.size() >= chain_long)
						continue;
					else
					sub_chain.push_back(event);
					i++;
				}
			}
			if(int(event.x*10.0) < 1590)
			{
				alpha_all[int(event.x*10.0)+5][int(event.y*10)].insert(pair<ULong64_t, Alpha>(event.ts, alpha));
				i = 0;
				for(auto &sub_chain : multi_chain[int(event.x*10.0)+5][int(event.y*10)])
				{
					if(sub_chain.size() >= chain_long)
						continue;
					else
					sub_chain.push_back(event);
					i++;
				}
			}
			if(int(event.y*10) > 0)
			{
				alpha_all[int(event.x*10.0)][int(event.y*10)-5].insert(pair<ULong64_t, Alpha>(event.ts, alpha));
				i = 0;
				for(auto &sub_chain : multi_chain[int(event.x*10.0)][int(event.y*10)-5])
				{
					if(sub_chain.size() >= chain_long)
						continue;
					else
					sub_chain.push_back(event);
					i++;
				}
			}
			if(int(event.y*10) < 1590)
			{
				alpha_all[int(event.x*10.0)][int(event.y*10)+5].insert(pair<ULong64_t, Alpha>(event.ts, alpha));
				i = 0;
				for(auto &sub_chain : multi_chain[int(event.x*10.0)][int(event.y*10)+5])
				{
					if(sub_chain.size() >= chain_long)
						continue;
					else
					sub_chain.push_back(event);
					i++;
				}
			}
//			cout << "Decay over "<<endl;
		}
		if(jentry % 1000 == 0 || jentry == nentries)
		{
			printf("Loop process %.2f %%, %lld  / %lld \r", Double_t(jentry)/(nentries)*100., jentry, nentries);
			fflush(stdout);
		}
	}
	printf("\n");
}

void tree::fill()
{
//	cout << "fill" << endl;
	for(int i=0; i<1600; i = i+5)
	{
		for(int j=0; j<1600; j=j+5)
		{
			for(auto &sub_chain : multi_chain[i][j])
			{
				imp.x = -1;
				imp.y = -1;
				imp.e = -1.0;
				imp.emin = -1.0;
				imp.ts = -1;
				decay.clear();
				chain.clear();
//				for(auto it = sub_chain.begin())
				for(auto &ssub_chain : sub_chain)
				{
					event.x = ssub_chain.x;
					event.y = ssub_chain.y;
					event.e = ssub_chain.e;
					event.emin = ssub_chain.emin;
					event.ts = ssub_chain.ts;
					if(ssub_chain.ts == sub_chain[0].ts)
					{
						imp.x = ssub_chain.x;
						imp.y = ssub_chain.y;
						imp.e = ssub_chain.e;
						imp.emin = ssub_chain.emin; 
						imp.ts = ssub_chain.ts;
						auto itsb = alpha_all[int(imp.x*10.0)][int(imp.y*10.0)].lower_bound(event.ts-front_window);
						auto itse = alpha_all[int(imp.x*10.0)][int(imp.y*10.0)].lower_bound(event.ts);
						for(auto it=itsb; it != alpha_all[int(imp.x*10.0)][int(imp.y*10.0)].end(); it++)
						{
							if(it->first > imp.ts)
									break;
							event.x = (*it).second.x;
							event.y = (*it).second.y;
							event.e = (*it).second.e;
							event.emin = (*it).second.emin;
							event.ts = (*it).first;
							decay.push_back(event);
						}
						continue;
					}
					else
					{
						chain.push_back(event);
					}
					if((event.ts - imp.ts) < back_window)
					{
						decay.push_back(event);
					}
				}
				if(imp.ts > 0 && (chain.size() != 0 || decay.size() != 0))
					opt->Fill();
			}

		}
		if(i % 10 == 0 || i==1599)
		{
			printf("Fill  process %.2f %%, %d  / %d \r", Double_t(i)/Double_t(1599)*100.0, i, 1599);
			fflush(stdout);
		}
	}
	printf("\n");
}

void tree::BranchOpt()
{
	opt->Branch("imp", &imp);
	opt->Branch("chain", &chain);
	opt->Branch("decay", &decay);
}
