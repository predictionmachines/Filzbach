#include "stdafx_base.h"

#include "R_ext/Print.h"
#define printf Rprintf
#define CHECK(cond,msg) if (!(cond)) {error(msg);}
#define drand unif_rand
