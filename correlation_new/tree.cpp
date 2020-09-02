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
	mesum = 0;  //must
	xesum = 0;
	yesum = 0;
	desum = 0;
	xvec = 0;
	yvec = 0;
	boxvec = 0;
	pavec = 0;
	mwvec = 0;
	xavec = 0;
	xavec_addback = 0;
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
	ipt->SetBranchAddress("xa_addback", &xavec_addback);
	ipt->SetBranchAddress("gs", &gsvec);
}

void tree::Loop(TTree *opt_)
{	
	if (opt_==NULL) return;
	opt = opt_;
	BranchOpt();
	imp_all.clear();
	alpha_all.clear();
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
			imp.x = -1.0;
			imp.y = -1.0;
			imp.ex_max = 0.0;
			imp.ey_max = 0.0;
			imp.e_min = 0.0;
			imp.ge.clear();
			if((*xvec).size()==1)
				imp.x = (*xvec)[0].id;
			else
				imp.x = Float_t((*xvec)[0].id + (*xvec)[1].id)/2.0;
			if((*yvec).size()==1)
				imp.y = (*yvec)[0].id;
			else
				imp.y = Float_t((*yvec)[0].id + (*yvec)[1].id)/2.0;
			imp.ex_max = (*xvec)[0].e;
			imp.ey_max = (*yvec)[0].e;
			imp.e_min = desum;
			if((*gsvec).size()>0)
			{
				for(int i=0; i<(*gsvec).size(); i++)
				{
					ge_data.e = (*gsvec)[i].e;
					ge_data.id = (*gsvec)[i].id;
					imp.ge.push_back(ge_data);
				}
			}
			imp_all.insert(pair<ULong64_t, Imp>((*xvec)[0].ts, imp));
//			cout << "Recoil over" <<endl;
		}
		//Decay
		if(mesum==0 && pasum==0 && desum < 10000)
		{
//			cout << "Decay start" << endl;
			alpha.x = -1.0;
			alpha.y = -1.0;
			alpha.ex_max = 0.0;
			alpha.ey_max = 0.0;
			alpha.e_min = 0.0;
			alpha.box_e = 0.0;
			alpha.box_id = -1.0;
			alpha.xa.clear();
			alpha.xa_addback.clear();
			if((*xvec).size()==1)
				alpha.x = (*xvec)[0].id;
			else
				alpha.x = Float_t((*xvec)[0].id + (*xvec)[1].id)/2.0;
			if((*yvec).size()==1)
				alpha.y = (*yvec)[0].id;
			else
				alpha.y = Float_t((*yvec)[0].id + (*yvec)[1].id)/2.0;
			alpha.ex_max = (*xvec)[0].e;
			alpha.ey_max = (*yvec)[0].e;
			alpha.e_min = desum;
			if((*boxvec).size()>0)
			{
				alpha.box_e = (*boxvec)[0].e;
				alpha.box_id = (*boxvec)[0].id;
			}
			if((*xavec).size()>0)
			{
				for(int i = 0; i< (*xavec).size(); i++)
				{
					ge_data.e = (*xavec)[i].e;
					ge_data.id = (*xavec)[i].id;
					alpha.xa.push_back(ge_data);
				}
			}
			if((*xavec_addback).size()>0)
			{
				for(int i = 0; i< (*xavec_addback).size(); i++)
				{
					ge_data.e = (*xavec_addback)[i].e;
					ge_data.id = (*xavec_addback)[i].id;
					alpha.xa_addback.push_back(ge_data);
				}
			}
			alpha_all.insert(pair<ULong64_t, Alpha>((*xvec)[0].ts, alpha));
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
	ULong64_t i = 0;
	for(auto it_imp=imp_all.begin(); it_imp != imp_all.end(); it_imp++)
	{
//		cout << "Pick" << endl;
		chain.clear();
		decay.clear();
		auto it_alpha_begin = alpha_all.lower_bound(it_imp->first - front_window);
		auto it_alpha_end = alpha_all.lower_bound(it_imp->first + back_window);
		for(auto it_alpha=it_alpha_begin; it_alpha!=it_alpha_end; it_alpha++)
		{
			if(abs(it_alpha->second.x-it_imp->second.x)>0.5 || abs(it_alpha->second.y-it_imp->second.y)>0.5)
				continue;
			alpha_data.ts = it_alpha->first;
			alpha_data.x = it_alpha->second.x;
			alpha_data.y = it_alpha->second.y;
			alpha_data.ex_max = it_alpha->second.ex_max;
			alpha_data.ey_max = it_alpha->second.ey_max;
			alpha_data.e_min = it_alpha->second.e_min;
			alpha_data.box_e = it_alpha->second.box_e;
			alpha_data.box_id = it_alpha->second.box_id;
			alpha_data.xa.clear();
			if(it_alpha->second.xa.size() > 0)
			{
				for(int j=0; j<it_alpha->second.xa.size(); j++)
				{
					ge_data.e = it_alpha->second.xa[j].e;
					ge_data.id = it_alpha->second.xa[j].id;
					alpha_data.xa.push_back(ge_data);
				}
			}
			alpha_data.xa_addback.clear();
			if(it_alpha->second.xa_addback.size() > 0)
			{
				for(int j=0; j<it_alpha->second.xa_addback.size(); j++)
				{
					ge_data.e = it_alpha->second.xa_addback[j].e;
					ge_data.id = it_alpha->second.xa_addback[j].id;
					alpha_data.xa_addback.push_back(ge_data);
				}
			}
			decay.push_back(alpha_data);

			bool bfill = true;
			if(alpha_data.ts<=it_imp->first)
				continue;
			auto it_imp1 = it_imp;
			it_imp1++;
			if(chain.size()!=0)
				it_imp1 = imp_all.lower_bound(chain[int(chain.size())-1].ts);
			auto it_imp2 = imp_all.lower_bound(alpha_data.ts);
			for(auto it_imp3=it_imp1; it_imp3!=it_imp2; it_imp3++)
			{
				if(abs(it_imp3->second.x-it_imp->second.x)<=0.5 && abs(it_imp3->second.y-it_imp->second.y)<=0.5)
				{
					bfill = false;
					break;
				}
			}
			if(bfill)
			{
				chain.push_back(alpha_data);
			}
		}
		if(decay.size()>0 || chain.size()>0)
		{
			imp_data.ts = it_imp->first;
			imp_data.x = it_imp->second.x;
			imp_data.y = it_imp->second.y;
			imp_data.ex_max = it_imp->second.ex_max;
			imp_data.ey_max = it_imp->second.ey_max;
			imp_data.e_min = it_imp->second.e_min;
			imp_data.ge.clear();	
			if(it_imp->second.ge.size()>0)
			{
				for(int j=0; j<imp.ge.size(); j++)
				{
					ge_data.e = it_imp->second.ge[j].e;
					ge_data.id = it_imp->second.ge[j].id;
					imp_data.ge.push_back(ge_data);
				}
			}
			opt->Fill();
//			cout << "fill tree over" << endl;
		}
		i++;
		if(i % 10 == 0 || i==1599)
		{
			printf("Fill  process %.2f %%, %lld  / %lld \r", Double_t(i)/Double_t(imp_all.size())*100.0, i, ULong64_t(imp_all.size()));
			fflush(stdout);
		}
	}
	printf("\n");
}

void tree::BranchOpt()
{
	opt->Branch("imp", &imp_data);
	opt->Branch("chain", &chain);
	opt->Branch("decay", &decay);
}
