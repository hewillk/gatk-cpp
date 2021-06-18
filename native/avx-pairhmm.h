#include <cstdint>
#include "pairhmm_common.h"
#include "Context.h"

#include "avx-types.h"

#undef SIMD_ENGINE
#define SIMD_ENGINE avx

#include "avx-functions-float.h"
#include "avx-vector-shift.h"
#include "avx-pairhmm-template.h"

#include "avx-functions-double.h"
#include "avx-vector-shift.h"
#include "avx-pairhmm-template.h"