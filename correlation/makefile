
flag = -std=c++11

correlation: tree.h tree.cpp main.cpp LinkDict.cpp
		g++  $(flag) $^ -o $@ $(shell root-config --libs --cflags)

LinkDict.cpp : tree.h Linkdef.h
	    #root 6
		rootcling -f LinkDict.cpp  $^
		# root 5
		#rootcint -f LinkDict.cpp -c $^

clean:
	rm -rf correlation *Dict*
