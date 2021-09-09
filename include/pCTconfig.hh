#ifndef PCTCONFIG_h
#define PCTCONFIG_h
// Utility class to set the pCT_Preprocessing command-line defaults from a text
// file
// The text file format is a set of lines of format:
// key = value
// where key and value are text strings that can be interpreted as string,
// integer, or float
// Comment lines beginning with # can be embedded in the file

#include "Util.hh"
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <fstream> 
using namespace std;

class pCTconfig {
public:
  struct cfgItm {
    string key;
    string longKey;
    string type;
    int *valInt;
    float *valFloat;
    string *valString;
  };
  static inline pCTconfig* GetInstance() { return theConfig; }

  
  vector<cfgItm> itemList;
  string configFileName;
  map<string, string> item_str;
  map<string, int> item_int;
  map<string, float> item_float;

  // Functions
  pCTconfig(string fileName); // constructor
  void addItem(string longKey, string &value);
  void addItem(string longKey, int &value);
  void addItem(string longKey, float &value);
  void SetStageThreshold(int stage, float value); 

  int Configure();
 private:
  static pCTconfig* theConfig;
};
#endif
