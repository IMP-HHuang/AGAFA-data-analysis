//tree.h

#include <TROOT.h>
#include <TTree.h>
#include <vector>
#include <map>
#ifndef tree_h
#define tree_h

using namespace std;
typedef struct GE
{
	Double_t e;
	Int_t id;
}GE;
typedef vector<GE> Ge;
typedef struct Imp
{
	Float_t x;
	Float_t y;
	Double_t ex_max;
	Double_t ey_max;
	Double_t e_min;
	Ge ge;
}Imp;
typedef struct IMP
{
	ULong64_t ts;
	Float_t x;
	Float_t y;
	Double_t ex_max;
	Double_t ey_max;
	Double_t e_min;
	Ge ge;
}IMP;
typedef struct Alpha
{
	Float_t x;
	Float_t y;
	Double_t ex_max;
	Double_t ey_max;
	Double_t e_min;
	Double_t box_e;
	Int_t box_id;
	Ge xa;
	Ge xa_addback;
}Alpha;
typedef struct ALPHA
{
	ULong64_t ts;
	Float_t x;
	Float_t y;
	Double_t ex_max;
	Double_t ey_max;
	Double_t e_min;
	Double_t box_e;
	Int_t box_id;
	Ge xa;
	Ge xa_addback;
}ALPHA;
typedef struct dssd
{
	Double_t e;
	Int_t id;
	ULong64_t ts;
} dssd;


class tree 
{
	public:
			tree(TTree *ipt_)
			{
				ipt = ipt_;
				chain_long = 10;
				back_window = 5E10;
				front_window = 1.0E10;
				Init();
			}
			virtual ~tree() {};
			void Init();
			virtual void Loop(TTree *opt_);
			virtual void fill();
			virtual void BranchOpt();

			TTree *ipt;
			TTree *opt;
			vector<dssd>* xvec, *yvec, *boxvec, *pavec, *mwvec;
			vector<dssd> *xavec, *xavec_addback, *gsvec;
	 		double mesum;
			double xesum, yesum, desum;
			multimap<ULong64_t, Imp> imp_all;
			multimap<ULong64_t, Alpha> alpha_all;
			Alpha alpha;
			Imp imp;
			IMP imp_data;
			ALPHA alpha_data;
			GE ge_data;
			vector<ALPHA> chain;
			vector<ALPHA> decay;
			dssd data;
			int chain_long;
			ULong64_t back_window;
			ULong64_t front_window;

};
#endif
