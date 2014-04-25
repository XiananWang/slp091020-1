/*
 * logMMSE_functions.h
 *
 *  Created on: 2014Äê4ÔÂ25ÈÕ
 *      Author: Timberjack
 */

#ifndef LOGMMSE_FUNCTIONS_H_
#define LOGMMSE_FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mel_functions.h"
#include "Transform.h"

void logMMSE(float* inputBuffer,int sr,int winduration,int hoptime);


#endif /* LOGMMSE_FUNCTIONS_H_ */
