//tree.h

#include <TROOT.h>
#include <TTree.h>
#include <vector>
#include <map>
#ifndef tree_h
#define tree_h

using namespace std;

typedef struct Event
{
	Float_t x;
	Float_t y;
	Double_t e;
	Double_t emin;
	ULong64_t ts;
}Event;

typedef struct Alpha
{
	Float_t x;
	Float_t y;
	Double_t e;
	Double_t emin;
}Alpha;
typedef multimap<ULong64_t, Alpha> Alpha_all;

typedef vector<Event> Chain;
typedef vector<Chain> Multi_chain;

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
				back_window = 1.02E11;
				front_window = 2.0E10;
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
			vector<dssd> *xavec, *gsvec;
			double mesum;
			double xesum, yesum, desum;
			Multi_chain multi_chain[1600][1600];
			Alpha_all alpha_all[1600][1600];
			Alpha alpha;
			Event event;
			Chain chain;
			dssd data;
			Event imp;
			Chain decay;
			int chain_long;
			ULong64_t back_window;
			ULong64_t front_window;

};
#endif
