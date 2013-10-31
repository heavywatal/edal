// -*- mode: c++; coding: utf-8 -*-
#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include <vector>

class Patch {
  public:
    double height() const {return height_;}
    double diameter() const {return diameter_;}

  private:
    double height_;
    double diameter_;
};




class Space {
  public:
    
  private:
    std::vector<std::vector<Patch> > matrix_;
};

#endif /* ENVIRONMENT_H_ */
