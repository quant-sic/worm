/* $Id: worm.hpp,v 1.1 2006/09/09 9:21:44 pollet Exp $ */

#ifndef worm_Element_HPP
#define worm_Element_HPP

#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <stdint.h>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include "model.hpp"

class Element;  // forward declaration

#if defined DSTRUC_AVL
#include "avl_tree.hpp"

class AVL_diagram : public AVL::tree<Element> {
public:

    
    static void print_name(std::ostream& os) {
        os << "# config : DSTRUC                                   : AVL\n";

    }

    using AVL::tree<Element>::insert;
    inline iterator insert(const_iterator pos, const Element& v) {
        return insert(v);
    }
    inline iterator insert(const_iterator pos, Element&& v) {
        return insert(std::move(v));
    }
};

typedef AVL_diagram Diagram_type;


#elif defined DSTRUC_LIST
#include <list>

class list_diagram : public std::list<Element > {
public:

    
    static void print_name(std::ostream& os) {
        os << "# config : DSTRUC                                   : LIST\n";
    }

    iterator at(size_t i) {
        if (i > size())
            return end();
        iterator it = begin();
        for (; i > 0; --i)
            ++it;
        return it;
    }
    const_iterator at(size_t i) const {
        if (i > size())
            return end();
        const_iterator it = begin();
        for (; i > 0; --i)
            ++it;
        return it;
    }
};

typedef list_diagram Diagram_type;

#elif defined DSTRUC_LIST_STACK
#include "list_stack.hpp"

class list_stack_diagram : public list_stack<Element> {
public:

    static void print_name(std::ostream& os) {
        os << "# config : DSTRUC                                   : LIST_STACK\n";
    }
};

typedef list_stack_diagram Diagram_type;
#endif



class Element
{
 public:
    using StateType = model::StateType;
    using SiteIndex = size_t;

  Element() {};

  Element(const StateType n0, const StateType n1, const SiteIndex s0, const double t, const int color, const size_t assoc_size) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    mLink = s0;
    mColor = color;
    mAssoc.resize(assoc_size);
  }

  ~Element() {}
  Element(const Element& src) 
    : mTime(src.mTime)
    , mBefore(src.mBefore)
    , mAfter(src.mAfter)
    , mLink(src.mLink)
    , mColor(src.mColor)
  {
    mAssoc.resize(src.mAssoc.size());
    for (size_t i = 0; i < mAssoc.size(); i++)
      mAssoc[i] = src.mAssoc[i];
  }

  bool operator<(const Element& rhs) const {
      if (time() == rhs.time()) {
        return after() == rhs.before();
      }
      else {
        return time() < rhs.time();
      }
    } 

  Element& operator=(const Element rhs) {
    if (this == &rhs) return (*this);
    mTime = rhs.mTime;
    mBefore = rhs.mBefore;
    mAfter = rhs.mAfter;
    mLink = rhs.mLink;
    mColor = rhs.mColor;
    if (mAssoc.size() != rhs.mAssoc.size())
        mAssoc.resize(rhs.mAssoc.size());
    for (size_t i = 0; i < mAssoc.size(); i++)
        mAssoc[i] = rhs.mAssoc[i];
    return *this;
  }

  friend std::ostream &operator<<( std::ostream &os, const Element& rhs) {
    os  << rhs.color() << "\t"
        << rhs.time() << "\t"
        << rhs.link() << "\t"
        << rhs.before() << "\t"
        << rhs.after();
    return (os);
  }

  friend std::istream& operator>>(std::istream& is, Element& e) {
    is >> e.mColor>> e.mTime >> e.mLink >> e.mBefore >> e.mAfter;
    return (is);
  }

  friend bool operator== (const Element& lhs, const Element& rhs) {
    return  ( (lhs.mColor == rhs.mColor) && (lhs.mLink == rhs.mLink) && (lhs.mTime == rhs.mTime) && (lhs.mBefore == rhs.mBefore) && (lhs.mAfter == rhs.mAfter) );
  }


  double time() const {return mTime;}
  StateType before() const {return mBefore;}
  StateType after() const {return mAfter;}
  int color() const {return mColor;}
  SiteIndex link() const { return mLink;}
  
  void time(const double t) {mTime=t;};
  void before(const StateType n) {mBefore = n;};
  void after(const StateType n) {mAfter = n;};
  void link(const SiteIndex l)  {mLink = l;}
  void color(const int col) {mColor = col;}

  void print() const {
    std::cout << "\n" << this << "\tcolor : " << mColor << "\ttime : " 
              << mTime << "\tlink : " << mLink << "\tbefore : " << mBefore << "\tafter : " << mAfter;
    for (auto const& dit : mAssoc) {
        std::cout << "\t" << dit->time() * dit->color();
    }
  }

  Diagram_type::iterator& get_assoc(SiteIndex s) { return mAssoc[s];}
  double get_assoc_time(SiteIndex s) const {return mAssoc[s]->time();}
  void set_assoc(const SiteIndex j, const Diagram_type::iterator v) {
    mAssoc[j] = v;
  }

 private:
  double mTime;  // time of the interaction
  StateType mBefore;     // occupation on the site before the interaction
  StateType mAfter;      // occupation on the site after the interaction
  SiteIndex mLink;   // site to which element is linked
  std::vector<Diagram_type::iterator> mAssoc;   // associations; i.e. iterator to the kink on nb sites that are equal or just greater in time
  int mColor;
};






#endif
