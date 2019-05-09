namespace RooUnfolding {
  template<class Hist> Hist* createHist(const TVectorD& vec, const char* name, const char* title, const Variable<Hist>& x, bool overflow) { return createHist<Hist>(vec,name,title,std::vector<Variable<Hist>>{x},overflow); };
  template<class Hist> Hist* createHist(const TVectorD& vec, const TVectorD& errvec, const char* name, const char* title, const Variable<Hist>& x, bool overflow) { return createHist<Hist>(vec,errvec,name,title,std::vector<Variable<Hist>>{x},overflow); };
  template<class Hist> std::vector<RooUnfolding::Variable<Hist>> vars(const Hist* h){
    int d = dim(h);
    std::vector<Variable<Hist> > v;
    v.push_back(var(h,X));
    if(d>1) v.push_back(var(h,Y));
    if(d>2) v.push_back(var(h,Z));
    return v;
  }
}
