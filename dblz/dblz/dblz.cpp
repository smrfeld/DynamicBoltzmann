//
//  dblz.cpp
//  dblz
//
//  Created by Oliver Ernst on 1/11/19.
//  Copyright © 2019 oke. All rights reserved.
//

#include <iostream>
#include "dblz.hpp"
#include "dblzPriv.hpp"

void dblz::HelloWorld(const char * s)
{
    dblzPriv *theObj = new dblzPriv;
    theObj->HelloWorldPriv(s);
    delete theObj;
};

void dblzPriv::HelloWorldPriv(const char * s) 
{
    std::cout << s << std::endl;
};

