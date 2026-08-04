/*
 * Uranus Atmosphere Circulation Model (ATURAN)
 * Dry convective adjustment — ATURAN's binding of the SHARED implementation.
 *
 * The algorithm, the reasoning and the knobs live in ConvectiveAdjustment.h, which is
 * byte-identical to ATSAT's, ATJUP's and ATNEPT's copies and knows nothing about which planet it
 * runs on. This file exists only so the call sites keep their familiar name.
 */

#pragma once

#include "ConvectiveAdjustment.h"

class cUranusModel;

typedef ConvectiveAdjustment<cUranusModel> ConvectiveAdjustmentUran;
