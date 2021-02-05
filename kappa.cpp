#include <Rcpp.h>
#include <sstream>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace std;

typedef vector<string> stringList;
typedef vector<int> numList;

// [[Rcpp::export]]
stringList split_(const string &s, char delim)
{
    stringList result;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim))
    {
        result.push_back(item);
    }
    return result;
}
// [[Rcpp::export]]
SEXP kappa(const string &x,const string &y, StringVector &g) {
    StringVector xx = Rcpp::wrap(split_(x,','));
    StringVector yy = Rcpp::wrap(split_(y,','));
    float kap;
    float Oab;
    float Aab;
    StringVector t1 = intersect(xx,yy);
    StringVector t2 = setdiff(xx,yy);
    StringVector t3 = setdiff(yy,xx);
    StringVector uu = union_(xx,yy);
    StringVector t4 = setdiff(g,uu);
    int l1 = t1.size();
    int l2 = t2.size();
    int l3 = t3.size();
    int l4 = t4.size();
    if(l1 == 0){
        kap = 0.0;
    }else{
        float tol = float(l1+l2+l3+l4);
        Oab = (l1+l4)/tol;
        Aab = ((l1+l2)*(l1+l3)+(l3+l4)*(l2+l4))/(tol*tol);
        if(Aab==1.0){
            kap = 0.0;
        }else{
            kap = (Oab-Aab)/(1-Aab);
        }
    }
    return(Rcpp::wrap(kap));
}

