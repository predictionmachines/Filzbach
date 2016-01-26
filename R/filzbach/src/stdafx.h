#include "stdafx_base.h"

#include "R_ext/Print.h"
#define printf Rprintf
#include "R.h"
#undef CHECK
#define CHECK(cond,msg) if (!(cond)) {error(msg);}
#define drand unif_rand
