//=================================================================================================
/*!
//  \file blaze/util/mpl/Xnor.h
//  \brief Header file for the Xnor class template
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_UTIL_MPL_XNOR_H_
#define _BLAZE_UTIL_MPL_XNOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "../../util/mpl/Bool.h"


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time logical 'not xor' evaluation.
// \ingroup mpl
//
// The Xor alias declaration performs at compile time a logical 'not xor' evaluation of the two
// given compile time conditions:

   \code
   using namespace blaze;

   using Type = int;

   Xnor< IsSigned<Type>  , IsIntegral<Type>      >::value  // Evaluates to 1
   Xnor< IsUnsigned<Type>, IsFloatingPoint<Type> >::value  // Evaluates to 1
   Xnor< IsSigned<Type>  , IsUnsigned<Type>      >::value  // Evaluates to 0
   Xnor< IsIntegral<Type>, IsFloatingPoint<Type> >::value  // Evaluates to
   \endcode
*/
template< typename T1    // Type of the first operand
        , typename T2 >  // Type of the second operand
struct Xnor
   : public Bool< !( T1::value ^ T2::value ) >
{};
//*************************************************************************************************

} // namespace blaze

#endif
