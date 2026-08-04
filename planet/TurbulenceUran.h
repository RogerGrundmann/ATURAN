/*
 * Uranus Atmosphere Circulation Model (ATURAN)
 * Turbulence closure — ATURAN's binding of the SHARED implementation.
 *
 * The closure lives in Turbulence.h, byte-identical to ATSAT's, ATJUP's and ATNEPT's copies.
 */

#pragma once

#include "Turbulence.h"

class cUranusModel;

typedef Turbulence<cUranusModel> TurbulenceUran;
