//=================================================================================================
/*!
//  \file blaze/math/functors/Min.h
//  \brief Header file for the Min functor
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

#ifndef _BLAZE_MATH_FUNCTORS_MIN_H_
#define _BLAZE_MATH_FUNCTORS_MIN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "../../math/constraints/SIMDPack.h"
#include "../../math/simd/Min.h"
#include "../../math/typetraits/HasSIMDMin.h"
#include "../../math/typetraits/IsLower.h"
#include "../../math/typetraits/IsStrictlyLower.h"
#include "../../math/typetraits/IsStrictlyUpper.h"
#include "../../math/typetraits/IsSymmetric.h"
#include "../../math/typetraits/IsUniLower.h"
#include "../../math/typetraits/IsUniUpper.h"
#include "../../math/typetraits/IsUpper.h"
#include "../../math/typetraits/YieldsLower.h"
#include "../../math/typetraits/YieldsStrictlyLower.h"
#include "../../math/typetraits/YieldsStrictlyUpper.h"
#include "../../math/typetraits/YieldsSymmetric.h"
#include "../../math/typetraits/YieldsUniLower.h"
#include "../../math/typetraits/YieldsUniUpper.h"
#include "../../math/typetraits/YieldsUpper.h"
#include "../../util/algorithms/Min.h"
#include "../../system/Inline.h"
#include "../../util/IntegralConstant.h"


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for the min() function.
// \ingroup functors
*/
struct Min
{
   //**********************************************************************************************
   /*!\brief Default constructor of the Min functor.
   */
   explicit inline Min()
   {}
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the min() function for the given objects/values.
   //
   // \param a The left-hand side object/value.
   // \param b The right-hand side object/value.
   // \return The result of the min() function for the given objects/values.
   */
   template< typename T1, typename T2 >
   BLAZE_ALWAYS_INLINE decltype(auto) operator()( const T1& a, const T2& b ) const
   {
      return min( a, b );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data types \a T1 and \a T2.
   //
   // \return \a true in case SIMD is enabled for the data types \a T1 and \a T2, \a false if not.
   */
   template< typename T1, typename T2 >
   static constexpr bool simdEnabled() { return HasSIMDMin_v<T1,T2>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the min() function for the given SIMD vectors.
   //
   // \param a The left-hand side SIMD vector.
   // \param b The right-hand side SIMD vector.
   // \return The result of the min() function for the given SIMD vectors.
   */
   template< typename T1, typename T2 >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const T1& a, const T2& b ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T1 );
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T2 );
      return min( a, b );
   }
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsSymmetric<Min,MT1,MT2>
   : public BoolConstant< IsSymmetric_v<MT1> && IsSymmetric_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsLower<Min,MT1,MT2>
   : public BoolConstant< IsLower_v<MT1> && IsLower_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsUniLower<Min,MT1,MT2>
   : public BoolConstant< IsUniLower_v<MT1> && IsUniLower_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsStrictlyLower<Min,MT1,MT2>
   : public BoolConstant< ( IsStrictlyLower_v<MT1> && ( IsStrictlyLower_v<MT2> || IsUniLower_v<MT2> ) ) ||
                          ( IsStrictlyLower_v<MT2> && ( IsStrictlyLower_v<MT1> || IsUniLower_v<MT1> ) ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsUpper<Min,MT1,MT2>
   : public BoolConstant< IsUpper_v<MT1> && IsUpper_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsUniUpper<Min,MT1,MT2>
   : public BoolConstant< IsUniUpper_v<MT1> && IsUniUpper_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsStrictlyUpper<Min,MT1,MT2>
   : public BoolConstant< ( IsStrictlyUpper_v<MT1> && ( IsStrictlyUpper_v<MT2> || IsUniUpper_v<MT2> ) ) ||
                          ( IsStrictlyUpper_v<MT2> && ( IsStrictlyUpper_v<MT1> || IsUniUpper_v<MT1> ) ) >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
