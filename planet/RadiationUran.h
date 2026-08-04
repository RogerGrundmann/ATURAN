/*
 * Uranus Atmosphere Circulation Model (ATURAN)
 * Grey multi-layer radiation — ATURAN's binding of the SHARED implementation.
 *
 * The scheme lives in Radiation.h, byte-identical to ATSAT's, ATJUP's and ATNEPT's copies.
 * Uranus's own numbers are in cUranusModel.h, where measured properties of the planet belong.
 */

#pragma once

#include "Radiation.h"

class cUranusModel;

typedef Radiation<cUranusModel> RadiationUran;
