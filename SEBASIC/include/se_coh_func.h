#ifndef SE_COH_FUNC_H
#define SE_COH_FUNC_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "se_basic_math.h"

float se_coh_func_semblance(const float *data, int n_tr, int n_samp);
float se_coh_func_variation(const float *data, int n_tr, int n_samp);

#endif